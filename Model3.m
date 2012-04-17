function dY=Model3(t,Y)
dY=zeros(2,1);

%load globals -----------------------
global modelAtm trackOmega trackCov trackCondSink

n = modelAtm.NumBins;
Pop = modelAtm.Pop;
%CStarBasis = modelAtm.CStarBasis;
alpha = modelAtm.SOA.alpha;

if modelAtm.PulseYN==0
    Tau = modelAtm.DecayConstant*3600;
    ROG = modelAtm.Injection/Tau*exp(-t/Tau);
else
    if t>modelAtm.EmitTime
        ROG = 0;
        modelAtm.Sulf.MassConc = eps;
        %xmodelAtm.Sulf.MassConc = 0;
    else
        ROG = modelAtm.ROG;
    end
end

Cpre = Y(n*(1+Pop)+2*Pop+1);
ROG = Cpre*1/(modelAtm.DecayConstant*3600);
dY(n*(1+Pop)+2*Pop+1) = -1*ROG;

%Load current values of integration variables -----------
for i=1:n
    Cv(i)=Y(i);  
    for j=1:Pop
        Cp(j,i)=Y(j*n+i);
    end
end

for j=1:Pop
    Dp(j) = Y((Pop+1)*n+j);
    M_Sulf(j) = Y((Pop+1)*n+Pop+j);
    M_SOA(j) = sum(Cp(j,:));
   % NumConc(j) = Y((Pop+1)*n+2*Pop+j);
end

%
FracSulfSusp = M_Sulf(1)/(M_Sulf(1)+M_SOA(1));

for j=1:Pop
    PopString = int2str(j);
    eval(['modelAtm.Pop' PopString '.Dp = Dp(j);'])
    eval(['NumConc(j) = modelAtm.Pop' PopString '.NumConc0;']);
end


%Load associated properties ----------

for j=1:Pop
    TotalMass(j) = M_Sulf(j)+M_SOA(j); %Total 1=# 2=SurfArea(m2) 3=Volume(m3) per m3
    rho = (M_Sulf(j)+M_SOA(j))/(M_Sulf(j)/modelAtm.Sulf.rho+M_SOA(j)/modelAtm.SOA.rho);
    MW = (M_Sulf(j)+M_SOA(j))/(M_Sulf(j)/modelAtm.Sulf.MW+M_SOA(j)/modelAtm.SOA.MW);
end

for i=1:n
  
   XCp(1,i) = Cp(1,i)/sum(Cp(1,:));
   XCv(i) = Cv(i)/sum(Cv);
   if sum(Cp(1,:)) == 0
       XCp(1,i) = 1;
       XCv(i) = 1;
   end
   
end

t

%Calculate condensation sinks------------

for j=1:Pop
    ParticleKn = 2*modelAtm.SOA.lambda/Dp(j);
    BetaP = FuchsC(ParticleKn);
    CondSinkTot(j) = alpha*NumConc(j)*2*pi*Dp(j)*modelAtm.SOA.Diff*BetaP;
    PopString = int2str(j);
    eval(['ParticleNum = modelAtm.Pop' PopString '.NumConc0;']) %No deposition or particle loss
    %ParticleNum = 1; %#/m3
    CondSinkP(j) = CondSinkTot(j)/ParticleNum;
end
%CondSinkTot(2) = CondSinkTot2(2);

In = length(trackCondSink(:,1))+1;
%trackCondSink(In,1) = t;
%trackCondSink(In,2) = CondSinkTot(1);
%trackCondSink(In,3) = CondSinkTot(2);
%trackCondSink(In,4) = TotalMass(1)/CondSinkTot(1);
%trackCondSink(In,5) = TotalMass(2)/CondSinkTot(2);
%trackCondSink(In,6) = TotalMass(1);
%trackCondSink(In,7) = TotalMass(2);

% %Consider aging of C* = 1 vapors
% ResTime = 2; %hour
CpAging = zeros(1,n);

%k_age = 3e-15;
%R_age = k_age*Cp(1,5)^2*3600;
%modelAtm.R_age = R_age;
%CpAging(1) = R_age;
%CpAging(5) = -R_age;
% CvAging(3) = -1*Cv(3)*1/(3600*ResTime);
% CvAging(1) = -1*CvAging(3);
% CvAging(2) = 0; 
% CvAging(4) = 0; 
% %CvAging = zeros(1,n); %No aging


%Calculates organic flux to bulk particles---------------------------
for j=1:Pop
    trackCov(length(trackCov(:,2)),1) = t;
    [dCv dCp] = ParticleGrowth(Cv,Cp(j,:),CondSinkTot(j),j,modelAtm.KelvinYN,ROG,Dp(j),t); 
    dCvStack(j,:)=dCv;
    for i=1:n
        if NumConc(j)>=1
            dY(j*n+i) = dCp(i);
        else
            dY(j*n+i) = 0;
        end
    end
    J_SOA(j) = sum(dCp);
    J_SOA_Cond(j) = sum(dCp);
end
   
%trackCondSink(In,8) = J_SOA(1)/CondSinkTot(1);
%trackCondSink(In,9) = J_SOA(2)/CondSinkTot(2);
%trackCondSink(In,10) = J_SOA(1);
%trackCondsink(In,11) = J_SOA(2)*1e10;

%itis = J_SOA(2)

for i=1:n
    dY(i) = sum(dCvStack(:,i));
end

     
%Calculate Sulfate growth on single particle and bulk growth-----------------------------------------     
  
for j=1:Pop
    J_Sulf_Cond(j) = CondSinkTot(j)*modelAtm.Sulf.MassConc;
    J_Sulf(j) = J_Sulf_Cond(j); %[=]ug/m3-s S&P eq12.12 p539
 
    dY((Pop+1)*n+Pop+j) = J_Sulf(j);

    rho_T(j) = (J_Sulf_Cond(j)+J_SOA_Cond(j))/(J_Sulf_Cond(j)/modelAtm.Sulf.rho+J_SOA_Cond(j)/modelAtm.SOA.rho);
    rho(j) = (M_Sulf(j)+M_SOA(j))/(M_Sulf(j)/modelAtm.Sulf.rho+M_SOA(j)/modelAtm.SOA.rho);
    rho_T(j) = rho(j);%%%%%%%%%%%%%% assumes density of condensing is same as density of particle itself
    
    J_T(j) = CondSinkP(j)/CondSinkTot(j)*(J_Sulf(j)+J_SOA(j));
    J_T_Cond(j) = CondSinkP(j)/CondSinkTot(j)*(J_Sulf_Cond(j)+J_SOA_Cond(j));
    
    dDp0(j) = (2*J_T_Cond(j))/(rho_T(j)*pi*Dp(j)^2); %growth rate of particles in population at start of timestep
    
    M_All(j) = M_SOA(j) + M_Sulf(j);
    J_All(j) = J_SOA(j) + J_Sulf(j);
    
    dY((Pop+1)*n+j)=dDp0(j);
end
    
%trackCondSink(In,8) = J_All(1)/CondSinkTot(1);
%trackCondSink(In,9) = J_All(2)/CondSinkTot(2);
%trackCondSink(In,10) = J_All(1);
%trackCondsink(In,11) = J_All(2)*1e10;



end    
    
