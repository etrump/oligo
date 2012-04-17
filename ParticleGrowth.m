function [dCv dCp] = ParticleGrowth(Cv, Cp, CondSink, j, KelvinYN,ROG,Dp,t)

global modelAtm

%KelvinYN = 1; %for testing

%j is the population number
n=modelAtm.NumBins;
CStarBasis = modelAtm.CStarBasis;
alpha = modelAtm.SOA.alphaProd;
dimerbin = 8;

Cp;
    
%total_cov = 0;
%for i=2:modelAtm.Pop
%   CovString = int2str(i);
%   eval(['total_cov = total_cov + modelAtm.Pop' CovString '.coverage;']);
%end

%CurPop = int2str(j);

%eval(['Dp = modelAtm.Pop' CurPop '.Dp;'])

%modelAtm.total_cov = total_cov;

%Consider aging of C* = 1 vapors
ResTime = modelAtm.AgingResTime*3600; %1/s
%CvAging(3) = -1*Cv(3)*1/(ResTime);
%CvAging(1) = -1*CvAging(3);
%CvAging(2) = 0; 
%CvAging(4) = 0; 

CpAging = zeros(1,n);
CvAging = zeros(1,n);

%k_age_f = 0.0464;
k_age_f = 0.093*1e10*16/(1*3600)*1e-3;  %THESE ONES ARE WORKING
k_age_r = 5e-5*1e10*16/(1*3600)*1e-3; %/s %THESE ONES ARE WORKING


K = 5000; % from a = 0.005 and b = 0.15.... K = b/a^2
k_age_r = 5e-5; %kr is fixed
k_age_f = K*k_age_r;

k_age_f = 0.1;
k_age_f = 0.5;

modelAtm.k_age_f = k_age_f;

if sum(Cp)>1e-5
    %R_age_f = k_age_f*(Cp(9)/sum(Cp))^2;
    R_age_f = k_age_f*Cp(dimerbin)^2;
    %R_age_r = k_age_r*(Cp(1)/sum(Cp));
    R_age_r = k_age_r*Cp(1);
else
    R_age_f = 0;
    R_age_r = 0;
end



%if Cp(9) < 0.1
%   R_age_f = 0; 
%end

%if Cp(1) < 0.1
%    R_age_r = 0;
%end

modelAtm.R_age_f = R_age_f;

CpAging(1) = 2*R_age_f-R_age_r;
CpAging(1) = R_age_f - R_age_r;
CpAging(dimerbin) = -1*CpAging(1);
CpAging(dimerbin) = 2*R_age_r - 2*R_age_f;
modelAtm.bla = CpAging;
bla = modelAtm.bla;
if (CpAging(dimerbin)*1 + Cp(dimerbin)) < 1e-5
    CpAging(dimerbin) = 0;
    CpAging(1) = 0;
end

if (CpAging(1)*1 + Cp(1)) < 1e-5
    CpAging(dimerbin) = 0;
    CpAging(1) = 0;
end



if modelAtm.AgingYN~=1
    CvAging = zeros(1,n); %No aging
    CpAging = zeros(1,n);
end

if KelvinYN==1
        Kelvin = KelvinTerm(j);
    else
        Kelvin = 1;     %do not consider Kelvin effect       
end
    

for i=1:n
    if j==1 %suspended particles
        P(i) = alpha(i)*ROG;  %Production rate of species i (vapor), [=] ug/m3-s
    else
        P(i) = 0;              % Only produce vapors once
    end

    if sum(Cp)<=0         %If nothing on the particle, set mol frac = 1
        MolFrac(i)=0;
    else
        MolFrac(i)=Cp(i)/(sum(Cp));  %Calculate mole fraction of organic i in particle phase
    end
    
%    if KelvinYN==1
 %       Kelvin(i) = KelvinTerm(j);
 %   else
 %       Kelvin(i) = 1;     %do not consider Kelvin effect       
 %   end
    
% 
%     Ceq(i) = Kelvin(i)*MolFrac(i)*CStarBasis(i);
%    
%     %jString = int2str(j);
%     %eval(['modelAtm.Pop' jString '.Ceq = Ceq(i);']);
     QPhi(i)=modelAtm.UptakeFact*CondSink*(Cv(i)-Kelvin*MolFrac(i)*CStarBasis(i)); % calculte flux assuming no other limitations [=]ug/m3-s
%     DriveForce(i) = Cv(i)-Ceq(i);
%     %eval(['modelAtm.Pop' jString '.DriveForce = DriveForce(i);']);
%     
%     tol = 1e-16;
%     tol = 0;
%     if(QPhi(i)>=0)  % Always allow organic flux toward bulk particles
%         Phi(i)=QPhi(i);
%         elseif (Cp(i)>abs(QPhi(i)))  % Allow flux away from particle if suffient mass on particle
%             Phi(i)=QPhi(i);
%             elseif (Cp(i)<=tol)  % When no particle-phase mass, nothing can flux away
%                 Phi(i)=0;
%                 else
%                     Phi(i)=-1*Cp(i);  % But can allow all of particle-phase mass to flux away
%     end
%        
%     %CvAging(i) = 0;
%     
%    if QPhi(i)>0 && QPhi(i)>Cv(i)
%         QPhi(i) = Cv(i);
%    end
    
%    tol = eps;
%    if (QPhi(i)<0) && ((Cp(i)+(QPhi(i))*10)<tol)
%       QPhi(i) = 0; 
%    end
  %%%  if (QPhi(i)<0) && ((Cp(i)+(QPhi(i))*1)<tol)
  %%%     QPhi(i) = 0; 
  %%%  end
    
   Phi(i) = QPhi(i);
   
   Dp_cut = 100;
  % if i == modelAtm.EmitBin
  % if ((Phi(i) < 0) && (Dp*1e9 < Dp_cut))
  %     Phi(i) = 0;
  % end
%   if ((Phi(i) < 0) && (Cp(i) < tol))
%       Phi(i) = 0;
%   end
  % end
   
   
  %  if (Phi(i)+Cp(i))<tol %%%%replaced 0
  %      dCv(i) = P(i)+CvAging(i);
  %      dCp(i) = 0+CvAging(i);
  %      dCp(i) = 0;  
  %  else
        dCv(i)=P(i)-Phi(i);
        dCp(i)=Phi(i)+CpAging(i);
  %      dCp(i) = Phi(i);
  %  end
    
      %Change in gas-phase concentration, or production minus flux to particles, [=] ug/m3-s
      %Change in particle-phase concentration, or flux to particle, [=] ug/m3-s
end
    


CpAging;
