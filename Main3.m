function [T, Y] = Main3(tspan, y0)

%*System Parameters*******************************


global modelAtm


% Set simulation parameters -----

%modelAtm.Pop = Pop;
n = modelAtm.NumBins;



fun=@Model3;

options=odeset('RelTol',1e-14,'AbsTol',1e-16,'NonNegative',[1:16]);
%options=odeset('RelTol',1e-12,'AbsTol',1e-14);
%options=odeset('RelTol',1e-11,'AbsTol',1e-13);
%options=odeset('RelTol',1e-9,'AbsTol',1e-12);
options=odeset('RelTol',1e-9,'AbsTol',1e-11,'NonNegative',[1:length(y0)],'MaxStep',10);
%options=odeset('RelTol',1e-7,'AbsTol',1e-10);
%options=odeset('RelTol',1e-6,'AbsTol',1e-9);
%options=odeset('RelTol',1e-5,'AbsTol',1e-8);
%options=odeset('RelTol',1e-1,'AbsTol',1e-1);
[T Y]=ode15s(fun,tspan,y0,options);





