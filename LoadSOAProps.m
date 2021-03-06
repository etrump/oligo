function LoadSOAProps()

global modelAtm

n = modelAtm.NumBins;

%Assign SOA Properties-----------------------------------------
modelAtm.SOA.MW = 152; %[=]g/mol
modelAtm.SOA.MW = 215;
modelAtm.SOA.MW = 200; % PAPER BASECASE D_K
%modelAtm.SOA.MW = 250; % PAPER MAX D_K
%modelAtm.SOA.MW = 100; % PAPER MIN D_K

modelAtm.SOA.rho = 1e12;              %[=]ug/m3, liquid
modelAtm.SOA.rho = 1.4e12; % PAPER BASECASE D_K
%modelAtm.SOA.rho = 1e12; %PAPER MAX D_K
%modelAtm.SOA.rho = 2e12; % PAPER MIN D_K

%%modelAtm.SOA.MassConc = 0.015;       %[=]ug/m3
modelAtm.SOA.CollDiam = 5.34e-8;     %[=]cm, with air
modelAtm.SOA.alpha = 1;              %sticking coeff
%%modelAtm.SOA.NumConc = MassConc_SOA*1e-12*1/(modelAtm.SOA.MW/6.023e23); %[=]#/cm3
z_SOA=modelAtm.SOA.MW/modelAtm.air.MW        %ratio of molecular masses, S&P p401
%modelAtm.SOA.lambda2=(1/(pi*(1+z_SOA)^(0.5)*modelAtm.air.NumConc*modelAtm.SOA.CollDiam^2))*(1/100);  %[=]m  S&P eq9.11 p401
modelAtm.SOA.Diff =(1e-4)*3/(8*pi)*(pi*(1.38065e-23)^3*modelAtm.Temp^3*(1+z_SOA)/...
    (2*(modelAtm.SOA.MW*1e-3/6.022e23)))^(1/2)/(modelAtm.Press*(modelAtm.SOA.CollDiam*1e-2)^2);
        %Diffusivity   [=]m2/s
modelAtm.SOA.Diff =3/(8*pi)*(pi*(1.38065e-23)^3*modelAtm.Temp^3*(1+z_SOA)/...
    (2*(modelAtm.SOA.MW*1e-3/6.022e23)))^(1/2)/(modelAtm.Press*(modelAtm.SOA.CollDiam*1e-2)^2);
modelAtm.SOA.va = (8*1.38e-23*298/(pi*modelAtm.SOA.MW*1e-3/6.022e23))^(1/2); %m/s
%modelAtm.SOA.va = 171

modelAtm.SOA.lambda = 3*modelAtm.SOA.Diff/modelAtm.SOA.va; %m
%modelAtm.SOA.lambda = 1e-7;
%  Use Fuchs & Sutugin version

modelAtm.SOA.sigma_guess = 0.05; % N/m &PAPERBASE CASE D_K
%%%%%modelAtm.SOA.sigma_guess = 0.03;
%modelAtm.SOA.sigma_guess = 0.1; % PAPER MAX D_K
%modelAtm.SOA.sigma_guess = 0.01; % PAPER MIN D_K



modelAtm.SOA.alphaProd = zeros(1,n);
modelAtm.SOA.alphaProd(modelAtm.EmitBin) = 1;
modelAtm.SOA.InitialComp = zeros(1,n);
modelAtm.SOA.InitialComp(modelAtm.BGBin) = 1;
%modelAtm.SOA.alphaProdAlphaP = [0.0284 0 0.3617 0.6099]; 

modelAtm.UptakeFact = 0.1;
%modelAtm.UptakeFact = 1;

if modelAtm.AlphaPYNVap == 1
    modelAtm.SOA.alphaProd = [0.0284 0 0.3617 0.6099]; %This is the alpha-pinene distribution (considering only first 4 bins)
    %modelAtm.SOA.alphaProd = [0 0 0.8 0.1]; %For testing!
    modelAtm.SOA.alphaProd = [0 0.004 0.000 0.051 0.086 0.120 0.183 0.400];
    modelAtm.SOA.alphaProd = [0 0.004 0.000 0.051 0.086 0.120 0.183 0.400 0.350 0.200];
    %modelAtm.SOA.alphaProd = [0 0.004 0.000 0.051 0.086-0.043 0.120 0.183 0.400 0.350 0.200+0.043];
     modelAtm.SOA.alphaProd = [0 0.004 0.000 0.051 0.086 0.120-0.06 0.183 0.400 0.350 0.200+0.06];
     modelAtm.SOA.alphaProd = [0 0.004 0.000 0.051 0.086-0.029 0.120-0.09 0.183 0.400 0.350 0.200+0.09+0.029];
    %modelAtm.SOA.alphaProd = [0 0.004 0.000 0.051 0.086 0.120 0.183 0.400 0.350 0.200];
    modelAtm.SOA.alphaProd = [0 0.004 0.000 0.051 0.086-0.029 0.120-0.12 0.183-0.12 0.400 0.350 0.200+0.12+0.12+0.029];
    modelAtm.SOA.alphaProd = [0 0.004 0.000 0.051 0.086-0.029 0.120-0.12 0.183-0.12 0.400+0.12+0.12+0.029 0.350 0.200];
    modelAtm.SOA.alphaProdAlphaP = modelAtm.SOA.alphaProd;
end


end

