function LoadInitWallPopulation(Total, Tag, Diam)

global modelAtm

if modelAtm.AmmonSeedYN == 1
    FracOrg = 0;
    rho = modelAtm.Sulf.rho;
else
    FracOrg = 1; %Particles are entirely organic
    rho = modelAtm.SOA.rho;
end


modelAtm.Pop2.Dp0 = Diam;

if Tag==1
    modelAtm.Pop2.NumConc0 = Total;
elseif Tag==2
    modelAtm.Pop2.NumConc0 = Total/rho*6/pi*1/Diam^3; 
else
    modelAtm.Pop2.NumConc0 = 1;
end

NumConc = modelAtm.Pop2.NumConc0;

[CondSinkBG BGCp MSulf] = LoadPopMono(1,NumConc,Diam,FracOrg,modelAtm.BGBin);

modelAtm.Pop2.CondSink0=CondSinkBG; 
modelAtm.Pop2.Cp0 = BGCp;
modelAtm.Pop2.MSulf0 = MSulf;
