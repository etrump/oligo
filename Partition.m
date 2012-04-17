function [Caer Cvap] = Partition(Mass_reacted)
%this function calculates partitioning.  C is the total gas + aerosol conc
%of species i in ug/m^3, Cstar is the saturation concentrations, Coa is a
%guess for the total organic aerosol concentration, maxiter is the maximum
%number of iterations, and tol is the solution tolerance in Coa.
%Mass_reacted = 1.66e3
global modelAtm


Cstar = modelAtm.CStarBasis
C = Mass_reacted*modelAtm.SOA.alphaProd*1/modelAtm.Vol
maxiter = 1e5;
tol = 1e-15;
Coa = 1; %initial guess


for n = 1:maxiter
    xi = (1 + Cstar/Coa ) .^ (-1);
    CoaNEW = sum(C .* xi);
    if abs(CoaNEW - Coa) < tol
        break
    else
        Coa = CoaNEW;
    end
end

if n == 1000
    error('Did not Converge')
end

Caer = C * CoaNEW ./ Cstar .* (1 + CoaNEW ./ Cstar) .^ (-1)
Cvap = C - Caer
