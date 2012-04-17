function LoadSecondPopulation(j,Total,Tag)

global modelAtm

%n = modelAtm.NumBins;
%Tag = 1;  %Tag =1, specify total number concentration
AmmonSeedYN = modelAtm.AmmonSeedYN;
AmmonSeedYN = 1;

if AmmonSeedYN == 1
    FracOrg = 0;
    rho = modelAtm.Sulf.rho;
else
    FracOrg = 1; %Particles are entirely organic
    rho = modelAtm.SOA.rho;
end
 %Total 1=# 2=SurfArea(m2) 3=Volume(m3) per m3
 
PopString = int2str(j);

eval(['modelAtm.Pop' PopString '.NumConc0 = 1;'])
eval(['Diam = modelAtm.Pop' PopString '.Dp0;'])
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

%[CondSinkBG BGCp MSulf fakevariable] = LoadPopMonoWall(1,NumConc,Diam,FracOrg,modelAtm.BGBin);
[CondSinkBG BGCp MSulf] = LoadPopMono(1,NumConc,Diam,FracOrg,modelAtm.BGBin);
%if j==2
%   BGCp = eps*ones(1,4); 
%end
eval(['modelAtm.Pop' PopString '.CondSink0=CondSinkBG;'])  
eval(['modelAtm.Pop' PopString '.Cp0 = BGCp;'])
eval(['modelAtm.Pop' PopString '.MSulf0 = MSulf;'])


