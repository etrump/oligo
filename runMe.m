clear
clc
close all
global modelAtm trackCov MassforPlot MassforPlotAlpha Struct
global Tagg CpT DpT  CvT Yagg_true
global trackCondSink

trackCondSink = zeros(1,11);

modelAtm.NumBins=10;
n = modelAtm.NumBins;

%Set simulation options: 1=yes, 0 = no ----------------------------
modelAtm.AgingYN = 1; 

modelAtm.AmmonSeedYN = 1;
modelAtm.AmmonSeedYNPop2 = 1;
modelAtm.PulseYN = 0;
modelAtm.KelvinYN = 1; %!!!!!!!
modelAtm.AlphaPYN = 0; %Alpha pinene organic = 1; single-bin organic = 0;
modelAtm.AlphaPYNPop2 = 0;
modelAtm.AlphaPYNVap = 1;
modelAtm.DilutionYN = 1;
modelAtm.DilutionTime = 5*3600;
modelAtm.DF = 150;
modelAtm.DF = 10; %120409

%Set simulation parameters ----------------------------
TimeVector = [0 10]*3600;
%TimeVector = [0 497.93]

EmitTime = 3*3600;
modelAtm.EmitTime = EmitTime;
modelAtm.DecayConstant = 0.171;     %hr-1

modelAtm.Injection = 1287;
modelAtm.Injection = 429; %120326
modelAtm.AgingResTime = 5;      %hr  

modelAtm.BGBin = 1; %Specify C* bin of initial particles; 1 is C*=1e-2, 2 is C* =1e-1, 3 is C*=1e0, 4 is C*=1\  /X\( .. )/X\ 
modelAtm.EmitBin = 3;

SulfMassConc = eps;

modelAtm.ROG = 1e-4;
DiamSusp = 10;  % in nanomaters 
DiamSuspM = DiamSusp*1e-9;

TotalSusp = 3e9; %exper on 120202
TotalSusp = 1e10*100;
TotalSusp = 4.5e11;
TotalSusp = 3.79e11;%120326
TotalSusp = 3.7943e+11; %120409
TotalSusp = 3.7943e+10; %120409 TESTING
TagSusp = 1; %1 = specify #/m3    2=specify ug/m3

TotalOnWall = eps;
TagOnWall = 1; %1 = specify #/m3    2=specify ug/m3

LoadAtmos;
LoadSOAProps;
modelAtm.SOA.lambda = 100e-9;
LoadSulfProps(SulfMassConc);
Sigma = 50;
modelAtm.TotalPopBG = 1;

if modelAtm.AmmonSeedYN == 1
    FracOrg = 0;
else
    FragOrg = 1;
end

modelAtm.Pop1.Dp0 = DiamSuspM; 
modelAtm.Pop1.NumConc0 = TotalSusp;

[CS0 Cp0BG MSulf0BG] = LoadPopBin(TotalSusp,DiamSuspM,FracOrg);
modelAtm.Pop1.Cp0 = Cp0BG;
modelAtm.Pop1.MSulf0 = MSulf0BG;
Dp0BG = DiamSusp;

for i = 1:modelAtm.TotalPopBG
    LoadPopulationProps(i); % what is this for?
end
%modelAtm.Cv01=LoadEqmVapors2(ceil(modelAtm.TotalPopBG/2+1));
%modelAtm.Cv0 = modelAtm.Cv01;
modelAtm.Cv0 = zeros(1,n);
Dp0 = Dp0BG;
Cp0 = Cp0BG;
MSulf0 = MSulf0BG;
NumConc0BG = TotalSusp;
modelAtm.CS0 = CS0;

modelAtm.V_small = 2;

Cpre0 = modelAtm.Injection/modelAtm.V_small; % m3;

modelAtm.Pop = 1+modelAtm.TotalPopBG;
modelAtm.Pop = 1;
modelAtm.total_cov = 0;
trackCov = zeros(1,4);    

Cp0_new = eps*ones(1,n);
Yagg = [];
Tagg = [];
Yagg_true = [];
Cv0 = modelAtm.Cv0;    

Dp0 = DiamSusp*1e-9;
    
    T_break = TimeVector(2);
    if modelAtm.DilutionYN == 1
        T_break = modelAtm.DilutionTime;
    end
    %tspan = [TimeVector(1) TimeVector(2)];
    tspan = [TimeVector(1) T_break];
    y0=[Cv0 Cp0 Dp0 MSulf0 Cpre0];
    
    [T1,Y1] = Main3(tspan,y0);
    
    if modelAtm.DilutionYN == 1
    tspan2 = [modelAtm.DilutionTime TimeVector(2)];
    y02 = Y1(length(T1),:);
    modelAtm.Pop1.NumConc0 = 1/modelAtm.DF*modelAtm.Pop1.NumConc0;
    
    for i = 1:n*(modelAtm.Pop+1)
        y02(i) = y02(i)*1/modelAtm.DF;
    end
    
    for i = n*(modelAtm.Pop+1)+modelAtm.Pop+1:n*(modelAtm.Pop+1)+modelAtm.Pop*2
        y02(i) = y02(i)*1/modelAtm.DF;
    end
    
    
    
    [T2,Y2] = Main3(tspan2,y02);
    Y2_true = Y2;
  for j = 1:length(T2)  
    for i = 1:n*(modelAtm.Pop+1)
        Y2(j,i) = Y2(j,i)*modelAtm.DF;
    end
    
    for i = n*(modelAtm.Pop+1)+modelAtm.Pop+1:n*(modelAtm.Pop+1)+modelAtm.Pop*2
        Y2(j,i) = Y2(j,i)*modelAtm.DF;
    end
  end
    
    for i = 1:(length(T1)+length(T2))
        if i <= length(T1)
            Tagg(i) = T1(i);
            Yagg(i,:) = Y1(i,:);
            Yagg_true(i,:) = Y1(i,:);
        else
            Tagg(i) = T2(i-length(T1));
            Yagg(i,:) = Y2((i-length(T1)),:);
            Yagg_true(i,:) = Y2_true((i-length(T1)),:);
        end
    end
    else
    Yagg = Y1;
    Yagg_true = Y1;
  
    Tagg = T1;
    end

    
Struct.Yagg = Yagg;
Struct.Tagg = Tagg;

ProcessData;

