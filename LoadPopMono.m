function [TotalCS BGCp TotalSulfMass] = LoadPopMono(Tag, TotalN, ModeM, FracOrg, Bin,CurrPop)

Mode = ModeM*1e9; %convert mode diameter from meters to nanometers
global modelAtm

rho_SOA = modelAtm.SOA.rho;
rho_Sulf = modelAtm.Sulf.rho;
rho = rho_SOA*FracOrg + rho_Sulf*(1-FracOrg);

lambda = modelAtm.SOA.lambda; %m
Diff = modelAtm.SOA.Diff; %m2/s

Dp_0 = 0.01; %[=]nm, start at 0.010 um
Dp_f = 1e5;   %end at 100 um

if CurrPop==2
    logSigma = 0.05;  % choose a standard deviation for the distribution (see S&P fig 8.3 pg371 for typical values)
else
    logSigma = 0.3;
    logSigma = 0.23;
    %logSigma = 2;
end
%logSigma = 0.12;
Sigma = 10^(logSigma);

% Load particle size range 
logDp_0 = log10(Dp_0);
logDp_f = log10(Dp_f);
logDp = logDp_0:0.001:logDp_f;
 for i = 1:length(logDp)
     Dp(i) = 10^(logDp(i));
 end

 Kn = zeros(length(Dp),1);
 logMode = log10(Mode);

 
for i = 1:length(Dp)-1
   % if Dp(i)==Mode
    d(i) = TotalN/((2*pi)^(1/2)*logSigma)*exp(-(logDp(i)-logMode)^2/(2*(logSigma)^2)); 
    %d(i) = Total; %S&P pg 370, eq 8.54
    Kn(i) = 2*lambda*1e9/Dp(i);
    Beta(i) = FuchsC(Kn(i));
    %Beta(i) = 1;
    if Tag==1  %Number distribution given    
        dCS(i) = d(i)*2*pi*Dp(i)*Diff*Beta(i)*1e-9; %[=] 1/s
        n_N(i) = d(i);  %#/m3 
        n_S(i) = d(i)*pi*Dp(i)^2*1e-18; %[=] m2/m3
        n_V(i) = d(i)*pi/6*Dp(i)^3*1e-27; %m3/m3 
    elseif Tag==2    
        dCS(i) = d(i)*2/Dp(i)*Diff*Beta(i)*1e9; %[=] 1/s
        n_N(i) = d(i)/(pi*Dp(i)^2)*1e18; %#/m3
        n_S(i) = d(i);  %m2/m3
        n_V(i) = d(i)*Dp(i)/6*1e-9;
    elseif Tag==3
        dCS(i) = d(i)*12/Dp(i)^2*Diff*Beta(i)*1e18;
        n_N(i) = d(i)*6/pi*1/Dp(i)^3*1e27;
        n_S(i) = d(i)*6/Dp(i)*1e9;
        n_V(i) = d(i);
    elseif Tag==4
        dCS(i) = d(i)/rho*12/Dp(i)^2*Diff*Beta(i)*1e18;
        n_N(i) = d(i)/rho*6/pi*1/Dp(i)^3*1e27;
        n_S(i) = d(i)/rho*6/Dp(i)*1e9;
        n_V(i) = d(i)/rho;
    end
        
    CS(i) = dCS(i)*(logDp(i+1)-logDp(i));
    D(i) = d(i)*(logDp(i+1)-logDp(i));
    V(i) = n_V(i)*(logDp(i+1)-logDp(i));
    N(i) = n_N(i)*(logDp(i+1)-logDp(i));
    S(i) = n_S(i)*(logDp(i+1)-logDp(i));

end  


  dCS(length(Dp)) = 0;
  d(length(Dp)) = 0;
  dN(length(Dp)) = 0;
  dV(length(Dp)) = 0;
  n_N(length(Dp))=0;

  figure
%semilogy(Dp,dCS)
semilogx(Dp,n_N)
%axis([0 1e3 0 1e9])
  
  n_N(length(Dp))=n_N(length(Dp)-1);
  n_S(length(Dp))=n_S(length(Dp)-1);
  n_V(length(Dp))=n_V(length(Dp)-1);
  
TotalCS = sum(CS);
TotalD = sum(D);  % check if program is working correctly
TotalV = sum(V);
TotalS = sum(S);
TotalN = sum(N);

TotalOrgMass = FracOrg*TotalV*rho_SOA;  % [=] ug/m3
TotalOrgMass = FracOrg*TotalV*rho;
TotalSulfMass = (1-FracOrg)*TotalV*rho_Sulf; % [=] ug/m3
TotalSulfMass = (1-FracOrg)*TotalV*rho;

CStarBasis = modelAtm.CStarBasis;

BGCp = zeros(1,4);   %BGCp [=] ug/m3


BGCp(Bin) = TotalOrgMass;

%if CurrPop==1 && modelAtm.AmmonSeedYN==0 && modelAtm.AlphaPYN==1
if CurrPop==1 && modelAtm.AlphaPYN==1
   BGCp = TotalOrgMass*modelAtm.SOA.alphaProdAlphaP;
end

%if CurrPop==2 && modelAtm.AmmonSeedYNPop2==0 && modelAtm.AlphaPYNPop2==1
if CurrPop==2 && modelAtm.AlphaPYNPop2==1
   BGCp = TotalOrgMass*modelAtm.SOA.alphaProdAlphaP; 
end




