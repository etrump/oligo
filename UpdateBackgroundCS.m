function [TotalCS] = UpdateBackgroundCS(TotalMass, Dpg0, rho_T, CurrPop)
%TotalMass in ug/m3



global modelAtm

rho_SOA = modelAtm.SOA.rho;
rho_Sulf = modelAtm.Sulf.rho; %ug/m3

rho = rho_SOA;

lambda = modelAtm.SOA.lambda; %m
Diff = modelAtm.SOA.Diff; %m2/s


if Dpg0<0 || Dpg0~=Dpg0
   display('Erorr!!!!!!!!!')
   modelAtm.OMG = Dpg0;
else
    modelAtm.OMG = 0;
end

Dpg = 1e9*Dpg0; %[=] nm
%rho = rho_T*1e-18; %[=] ug/um3

% mode is in nm
%Diff = 7e-2; %cm2/s
%lambda = 1.8e-8; % [=]m
%Integrate over entire particle distribution
Dp_0 = 1; %[=]nm, start at 0.010 um
Dp_f = 1e4;   %end at 100 um
%Dp_0 = 10;
%Dp_f = Dpg+300; %NOTE THIS COULD AUSE ERROR IF IT IS CUTTING OFF PART OF BACKGROUND DISTRIBUTION - ONLY FOR SMALL DEVIATIONS

logSigma = 0.23;  % choose a standard deviation for the distribution (see S&P fig 8.3 pg371 for typical values)
%logSigma = 0.01;
%logSigma = 0.12;

CurrPop = 1; %FOR TESTING ONLY, and to make similar to 8/24

if CurrPop==2
    logSigma = 0.05;  % choose a standard deviation for the distribution (see S&P fig 8.3 pg371 for typical values)
else
    logSigma = 0.3;
end

Sigma = exp(logSigma);


% Load particle size range 
logDp_0 = log(Dp_0);
logDp_f = log(Dp_f);
logDp = logDp_0:0.1:logDp_f;
 for i = 1:length(logDp)
     Dp(i) = exp(logDp(i));
 end

 %Dp;
 
 Kn = zeros(length(Dp),1);
 logDpg = log(Dpg);
 DpgS = exp(logDpg+2*logSigma^2);
 DpgV = exp(logDpg+3*logSigma^2);
 logDpgV = log(DpgV);
 
  for i = 1:length(Dp)-1
    d(i) = TotalMass/((2*pi)^(1/2)*logSigma)*exp(-(logDp(i)-logDpgV)^2/(2*(logSigma)^2));
    Kn(i) = 2*lambda*1e9/Dp(i);
    Beta(i) = FuchsC(Kn(i));

    dCS(i) = d(i)/rho*12/Dp(i)^2*Diff*Beta(i)*1e18;
    n_N(i) = d(i)/rho*6/pi*1/Dp(i)^3*1e27;
    n_S(i) = d(i)/rho*6/Dp(i)*1e9;
    n_V(i) = d(i)/rho;
        
    CS(i) = dCS(i)*(logDp(i+1)-logDp(i));
    D(i) = d(i)*(logDp(i+1)-logDp(i));
    V(i) = n_V(i)*(logDp(i+1)-logDp(i));
  end  

  dCS(length(Dp)) = 0;
  d(length(Dp)) = 0;
  
  %n_N;
  n_N(length(Dp))=n_N(length(Dp)-1);
  n_S(length(Dp))=n_S(length(Dp)-1);
  n_V(length(Dp))=n_V(length(Dp)-1);
  
TotalCS = sum(CS);
TotalD = sum(D);  % check if program is working correctly
TotalV = sum(V);

CheckMe = (TotalD-TotalMass)/TotalD;
modelAtm.CheckMe = CheckMe;