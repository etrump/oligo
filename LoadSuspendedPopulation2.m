function [CondSink0 Cp0 MSulf0 Dp0 NumConc0] = LoadSuspendedPopulation2(Total, Tag,DiamMode,Sigma)
%Total 1=# 2=SurfArea(m2) 3=Volume(m3) per m3
%Tag =1, specify total number concentration, 4 = mass

%Tag = 1;
%Total = 5e8;
%Tag = 4;
%Total = 10;
DiamMode = 50;
Sigma = 50;

global modelAtm

n = modelAtm.NumBins;
DpMode0 = DiamMode*1e-9;

%Define population parameters

if modelAtm.AmmonSeedYN == 1
    FracOrg = 0; 
    rho = modelAtm.Sulf.rho;
else
    FracOrg = 1; %Particles are entirely organic
    rho = modelAtm.SOA.rho;
    FracOrg = 0.5;
    rho = modelAtm.SOA.rho*FracOrg + modelAtm.Sulf.rho*(1-FracOrg);
end

if Tag==1
    N_tot = Total;
elseif Tag ==4
    N_tot = Total/rho*6/pi*1/DpMode0^3
end
N_tot

logSigma = log10(Sigma);
logMode = log10(DiamMode); %based on nm

Dp_0m = 10; %[=]nm, start at 0.010 um
Dp_fm = 1000;   % match SMPS

logDp_0m = log10(Dp_0m);
logDp_fm = log10(Dp_fm);
logDpm = logDp_0m:0.001:logDp_fm;
for i = 1:length(logDpm)
     Dpm(i) = 10^(logDpm(i));
end
 
for i = 1:length(Dpm)
    d(i) = N_tot/((2*pi)^(1/2)*logSigma)*exp(-(logDpm(i)-logMode)^2/(2*(logSigma)^2)); %dNdlogDp, fitted
    %d(i) = Total; %S&P pg 370, eq 8.54
    %Beta(i) = 1;
    n_N(i) = d(i);  %#/m3 
end

for i = 1:length(Dpm)-1
   N(i) =  (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2;
end
N(length(Dpm)) = 0;

figure(55)
semilogx(Dpm,n_N)

TotalPopBG = 5;
modelAtm.TotalPop = TotalPopBG+1;
modelAtm.TotalPopBG = TotalPopBG;
xn = Dp_fm;
x0 = Dp_0m;
delh = (log10(xn)-log10(x0))/TotalPopBG;
crap = 1:TotalPopBG;
Bin_split = 10.^(delh*crap + log10(x0));
Dp_bin_vect(1) = (x0+(Bin_split(1)))/2;
for i = 2:TotalPopBG
    Dp_bin_vect(i) = (Bin_split(i)+Bin_split(i-1))/2;
end

N_bin = zeros(TotalPopBG,length(Dpm)-1);

for i = 1:length(Dpm)-1
   N_check(i) =  (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2;
   
   dump = 0;
   for j=1:TotalPopBG
       if Dpm(i) < Bin_split(j) && dump~=1
           N_bin(j,i) = (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2;
           dump = 1;
       end
   end
   if dump == 0
       N_bin(1,i) = (log10(Dpm(i+1))-log10(Dpm(i)))*(n_N(i)+n_N(i+1))*1/2;
       dump = 1;
   end
   if N_bin(j,i) == 0
       N_bin(j,i) = 1;
   end
           
end

for j=1:TotalPopBG
    N_tot_bin(j) = sum(N_bin(j,:));
end

for j = 1:TotalPopBG-1
    nN_bin(j) = 1/2*(N_tot_bin(j+1)+N_tot_bin(j))*(log10(Dp_bin_vect(j+1))-log10(Dp_bin_vect(j))); 
end
nN_bin(TotalPopBG) = 0;

N_tot_bin;
N_tot;
sum(N_check);

%Binned_tot = sum(N_tot_bin)
Unbinned_tot = sum(N_check)
Binned_tot = sum(N_tot_bin)

figure(2)
%plot(Dp_bin_vect,nN_bin)
DiamSusp_m = 1e-9*Dp_bin_vect;
TotalSusp = N_tot_bin;
semilogx(Dp_bin_vect,N_tot_bin/max(N_tot_bin),Dpm(1:length(N)-1),N(1:length(N)-1)/max(N))

CondSink0 = 0;
for j = 1:TotalPopBG
%j = 1;
    %[BGCp_1bin MSulf_1bin] = LoadPopulation(TotalSusp(j),DiamSusp_m(j),j);
    [CondSink BGCp_1bin MSulf_1bin] = LoadPopBin(TotalSusp(j),DiamSusp_m(j),FracOrg);
    
    CondSink0 = CondSink0 + CondSink;
    for i = 1:n
     Cp0((j-1)*n+i) = BGCp_1bin(i); % only for one-bin organics!
    end
    
     MSulf0(j) = MSulf_1bin;
     Dp0(j) = DiamSusp_m(j);
     NumConc0(j) = TotalSusp(j);
       
    PopString = int2str(j);
    eval(['modelAtm.Pop' PopString '.Dp0 = DiamSusp_m(j);']) 
    eval(['modelAtm.Pop' PopString '.NumConc0 = TotalSusp(j);'])
    eval(['modelAtm.Pop' PopString '.CondSink0 = CondSink;'])
    eval(['modelAtm.Pop' PopString '.Cp0 = BGCp_1bin;'])
    
    %Cv0_contrib(j,:) = LoadEqmVapors2(j);
    LoadPopulationProps(j);
    
end

TotalSusp
FracOrg
Cp0

for i = 1:5
    ugh(i) = TotalSusp(i)*pi/6*(DiamSusp_m(i))^3*modelAtm.SOA.rho;
end
ugh
sum(ugh)
   