function [TotalCS BGCp TotalSulfMass] = LoadPopMono2Modes(Tag, TotalN, ModeM, TotalN2, ModeM2, TotalN3, ModeM3, FracOrg, Bin,CurrPop)

CurrPop = 1; %TESTING
Tag = 1; %TESTING

Mode = ModeM*1e9; %convert meters to nanometers
Mode2 = ModeM2*1e9;
Mode3 = ModeM3*1e9;
global modelAtm

rho_SOA = modelAtm.SOA.rho;
rho_Sulf = modelAtm.Sulf.rho;
rho = rho_SOA*FracOrg + rho_Sulf*(1-FracOrg);

lambda = modelAtm.SOA.lambda; %m
Diff = modelAtm.SOA.Diff; %m2/s

Dp_0 = 1; %[=]nm, start at 0.001 um
Dp_f = 1e5;   %end at 100 um

if CurrPop==2
    logSigma = 0.05;  % choose a standard deviation for the distribution (see S&P fig 8.3 pg371 for typical values)
    logSigma2= 0.05;
    logSigma3 = 0.05;
else
    logSigma = 0.245;
    %logSigma = 0.225;
    logSigma2 = 0.666;
    %logSigma2 = 0.557;
    logSigma3 = 0.337;
    %logSigma3 = 0.266;
end

Sigma = 10^(logSigma);
Sigma2 = 10^(logSigma2);
Sigma3 = 10^(logSigma3);

% Load particle size range 
logDp_0 = log10(Dp_0);
logDp_f = log10(Dp_f);
logDp = logDp_0:0.01:logDp_f;
logDp2 = logDp;
logDp3 = logDp;
 for i = 1:length(logDp)
     Dp(i) = 10^(logDp(i));
     Dp2(i) = Dp(i);
     Dp3(i) = Dp(i);
 end

 Kn = zeros(length(Dp),1);
 Kn2 = Kn;
 Kn3 = Kn;
 logMode = log10(Mode);
 logMode2 = log10(Mode2);
 logMode3= log10(Mode3);

 
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



for i = 1:length(Dp)-1
   % if Dp(i)==Mode
    d2(i) = TotalN2/((2*pi)^(1/2)*logSigma2)*exp(-(logDp2(i)-logMode2)^2/(2*(logSigma2)^2)); 
    %d(i) = Total; %S&P pg 370, eq 8.54
    Kn2(i) = 2*lambda*1e9/Dp2(i);
    Beta2(i) = FuchsC(Kn2(i));
    %Beta(i) = 1;
    if Tag==1  %Number distribution given    
        dCS2(i) = d2(i)*2*pi*Dp2(i)*Diff*Beta2(i)*1e-9; %[=] 1/s
        n_N2(i) = d2(i);  %#/m3 
        n_S2(i) = d2(i)*pi*Dp2(i)^2*1e-18; %[=] m2/m3
        n_V2(i) = d2(i)*pi/6*Dp2(i)^3*1e-27; %m3/m3 
    elseif Tag==2    
        dCS2(i) = d2(i)*2/Dp2(i)*Diff*Beta2(i)*1e9; %[=] 1/s
        n_N2(i) = d2(i)/(pi*Dp2(i)^2)*1e18; %#/m3
        n_S2(i) = d2(i);  %m2/m3
        n_V2(i) = d2(i)*Dp2(i)/6*1e-9;
    elseif Tag==3
        dCS2(i) = d2(i)*12/Dp2(i)^2*Diff*Beta2(i)*1e18;
        n_N2(i) = d2(i)*6/pi*1/Dp2(i)^3*1e27;
        n_S2(i) = d2(i)*6/Dp2(i)*1e9;
        n_V2(i) = d2(i);
    elseif Tag==4
        dCS2(i) = d2(i)/rho*12/Dp2(i)^2*Diff*Beta2(i)*1e18;
        n_N2(i) = d2(i)/rho*6/pi*1/Dp2(i)^3*1e27;
        n_S2(i) = d2(i)/rho*6/Dp2(i)*1e9;
        n_V2(i) = d2(i)/rho;
    end
        
    CS2(i) = dCS2(i)*(logDp2(i+1)-logDp2(i));
    D2(i) = d2(i)*(logDp2(i+1)-logDp2(i));
    V2(i) = n_V2(i)*(logDp2(i+1)-logDp2(i));
    N2(i) = n_N2(i)*(logDp2(i+1)-logDp2(i));
    S2(i) = n_S2(i)*(logDp2(i+1)-logDp2(i));

end  



for i = 1:length(Dp)-1
   % if Dp(i)==Mode
    d3(i) = TotalN3/((2*pi)^(1/2)*logSigma3)*exp(-(logDp3(i)-logMode3)^2/(2*(logSigma3)^2)); 
    %d(i) = Total; %S&P pg 370, eq 8.54
    Kn3(i) = 2*lambda*1e9/Dp3(i);
    Beta3(i) = FuchsC(Kn3(i));
    %Beta(i) = 1;
    if Tag==1  %Number distribution given    
        dCS3(i) = d3(i)*2*pi*Dp3(i)*Diff*Beta3(i)*1e-9; %[=] 1/s
        n_N3(i) = d3(i);  %#/m3 
        n_S3(i) = d3(i)*pi*Dp3(i)^2*1e-18; %[=] m2/m3
        n_V3(i) = d3(i)*pi/6*Dp3(i)^3*1e-27; %m3/m3 
    elseif Tag==2    
        dCS3(i) = d3(i)*2/Dp3(i)*Diff*Beta3(i)*1e9; %[=] 1/s
        n_N3(i) = d3(i)/(pi*Dp3(i)^2)*1e18; %#/m3
        n_S3(i) = d3(i);  %m2/m3
        n_V3(i) = d3(i)*Dp3(i)/6*1e-9;
    elseif Tag==3
        dCS3(i) = d3(i)*12/Dp3(i)^2*Diff*Beta3(i)*1e18;
        n_N3(i) = d3(i)*6/pi*1/Dp3(i)^3*1e27;
        n_S3(i) = d3(i)*6/Dp3(i)*1e9;
        n_V3(i) = d3(i);
    elseif Tag==4
        dCS3(i) = d3(i)/rho*12/Dp3(i)^2*Diff*Beta3(i)*1e18;
        n_N3(i) = d3(i)/rho*6/pi*1/Dp3(i)^3*1e27;
        n_S3(i) = d3(i)/rho*6/Dp3(i)*1e9;
        n_V3(i) = d3(i)/rho;
    end
        
    CS3(i) = dCS3(i)*(logDp3(i+1)-logDp3(i));
    D3(i) = d3(i)*(logDp3(i+1)-logDp3(i));
    V3(i) = n_V3(i)*(logDp3(i+1)-logDp3(i));
    N3(i) = n_N3(i)*(logDp3(i+1)-logDp3(i));
    S3(i) = n_S3(i)*(logDp3(i+1)-logDp3(i));

end  

  dCS(length(Dp)) = 0;
  d(length(Dp)) = 0;
  dN(length(Dp)) = 0;
  dV(length(Dp)) = 0;
  n_N(length(Dp))=0; 
  n_S(length(Dp))=0;
  
  dCS2(length(Dp2)) = 0;
  d2(length(Dp2)) = 0;
  dN2(length(Dp2)) = 0;
  dV2(length(Dp2)) = 0;
  n_N2(length(Dp2))=0;
  n_S2(length(Dp))=0;
  
   dCS3(length(Dp3)) = 0;
  d3(length(Dp3)) = 0;
  dN3(length(Dp3)) = 0;
  dV3(length(Dp3)) = 0;
  n_N3(length(Dp3))=0;
  n_S3(length(Dp))=0;

  figure(555)
semilogx(Dp,n_N+n_N2+n_N3)
%plot(Dp,n_N+n_N2+n_N3)
%axis([1 1e5 0 2e5])

figure(556)
semilogx(Dp,n_S+n_S2+n_S3)
%plot(Dp,n_S+n_S2+n_S3)

  n_N(length(Dp))=n_N(length(Dp)-1);
  n_S(length(Dp))=n_S(length(Dp)-1);
  n_V(length(Dp))=n_V(length(Dp)-1);
  
  n_N2(length(Dp2))=n_N2(length(Dp2)-1);
  n_S2(length(Dp2))=n_S2(length(Dp2)-1);
  n_V2(length(Dp2))=n_V2(length(Dp2)-1);
  
  
  n_N3(length(Dp3))=n_N3(length(Dp3)-1);
  n_S3(length(Dp3))=n_S3(length(Dp3)-1);
  n_V3(length(Dp3))=n_V3(length(Dp3)-1);
  
TotalCS1 = sum(CS);
TotalD1 = sum(D);  % check if program is working correctly
TotalV1 = sum(V);
TotalS1 = sum(S);
TotalN1 = sum(N);

  
TotalCS2 = sum(CS2);
TotalD2 = sum(D2);  % check if program is working correctly
TotalV2 = sum(V2);
TotalS2 = sum(S2);
TotalN2 = sum(N2);


TotalCS3 = sum(CS3);
TotalD3 = sum(D3);  % check if program is working correctly
TotalV3 = sum(V3);
TotalS3 = sum(S3);
TotalN3 = sum(N3);

TotalCS = TotalCS1 + TotalCS2+ TotalCS3;
TotalD = TotalD1+TotalD2+TotalD3;
TotalV = TotalV1+TotalV2+TotalV3;
TotalS = TotalS1+TotalS2+TotalS3;
TotalN = TotalN1+TotalN2+TotalN3;

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




