function LoadSecondPopulation(j,Total,Tag,Dp)
%Tag specifies units of "Total" 1=# 2=SurfArea(m2) 3=Volume(m3) all per m3

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

eval(['modelAtm.Pop' PopString '.NumConc0 = 1;'])
eval(['modelAtm.Pop' PopString '.Dp0 = Diam;'])
if j==2
   if Tag==1
       modelAtm.Pop2.NumConc0 = Total;
   elseif Tag==2
       modelAtm.Pop2.NumConc0 = Total/rho*6/pi*1/Diam^3; 
   else
       modelAtm.Pop2.NumConc0 = 1; 
   end
end


eval(['NumConc = modelAtm.Pop' PopString '.NumConc0;'])
eval(['modelAtm.Pop' PopString '.NumConc = modelAtm.Pop' PopString '.NumConc0;'])

[CondSinkBG BGCp MSulf] = LoadPopMono(1,NumConc,Diam,FracOrg,modelAtm.BGBin,2);

eval(['modelAtm.Pop' PopString '.CondSink0=CondSinkBG;'])  
eval(['modelAtm.Pop' PopString '.Cp0 = BGCp;'])
eval(['modelAtm.Pop' PopString '.MSulf0 = MSulf;'])


