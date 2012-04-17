function [Cp0 MSulf0] = LoadSecondPopulation2(j,Total,Dp)
%Tag specifies units of "Total" 1=# 2=SurfArea(m2) 3=Volume(m3) all per m3
Tag = 1;
NumConc = Total;
global modelAtm
Diam = Dp*1e-9;
%j = 2
%Total = 2
%Tag = 2
AmmonSeedYN = modelAtm.AmmonSeedYNPop2;

if AmmonSeedYN == 1
    FracOrg = 0;
    FracOrg = 0;
    rho = modelAtm.Sulf.rho;
else
    FracOrg = 1; %Particles are entirely organic
    rho = modelAtm.SOA.rho;
end
 
PopString = int2str(j);

eval(['modelAtm.Pop' PopString '.NumConc0 = NumConc;'])
eval(['modelAtm.Pop' PopString '.Dp0 = Diam;'])

%[CondSinkBG BGCp MSulf] = LoadPopMono(1,NumConc,Diam,FracOrg,modelAtm.BGBin,2);
[CondSinkBG BGCp MSulf] = LoadPopBin(NumConc,Diam,FracOrg);

Cp0 = BGCp;
MSulf0 = MSulf;


eval(['modelAtm.Pop' PopString '.CondSink0=CondSinkBG;'])  
eval(['modelAtm.Pop' PopString '.Cp0 = BGCp;'])
eval(['modelAtm.Pop' PopString '.MSulf0 = MSulf;'])


