function LoadPopulations(Pop)

global modelAtm

%Define population parameters
modelAtm.Pop1.NumConc0 = 5*1e9;
modelAtm.Pop1.Dp0 = 50; %Mode of # distribution, (in nm)
% %modelAtm.Pop2.NumConc0 = 2*1e9;
% modelAtm.Pop2.NumConc0 = eps;
% modelAtm.Pop2.Dp0 = 50;
% modelAtm.Pop3.NumConc0 = eps;
% modelAtm.Pop3.Dp0 = 50;

n = modelAtm.NumBins;
Tag = 1;  %Tag =1, specify total number concentration
FracOrg = 1; %Particles are entirely organic
 %Total 1=# 2=SurfArea(m2) 3=Volume(m3) per m3
 
for j = 1:Pop
    PopString = int2str(j);
%     if j>1  %Wall populations intially empty
%         eval(['modelAtm.Pop' PopString '.NumConc0 = eps;']);
%         eval(['modelAtm.Pop' PopString '.Dp0 = 50;']);
%     end
    eval(['NumConc = modelAtm.Pop' PopString '.NumConc0;'])
    eval(['Diam = modelAtm.Pop' PopString '.Dp0;'])
    [CondSinkBG BGCp MSulf] = LoadPopMono3(Tag,NumConc,Diam,FracOrg);
    eval(['modelAtm.Pop' PopString '.CondSink0=CondSinkBG;'])  
    eval(['modelAtm.Pop' PopString '.Cp0 = BGCp;'])
    eval(['modelAtm.Pop' PopString '.MSulf0 = MSulf;'])
    eval(['modelAtm.Pop' PopString '.coverage = 0;'])
end

CpTot = zeros(1,n);
for i = 1:n  
    CpTot(i) = modelAtm.Pop1.Cp0(i)+modelAtm.Pop2.Cp0(i);  
end

modelAtm.Cv0 = LoadEqmVapors(CpTot);