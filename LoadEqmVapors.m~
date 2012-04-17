function LoadEqmVapors()

global modelAtm

BGCp = modelAtm.Pop1.Cp0;

CStarBasis = modelAtm.CStarBasis;
n = modelAtm.NumBins;

for i=1:n
   if modelAtm.KelvinYN==1
        Kelvin(i) = exp(4*modelAtm.SOA.MW*modelAtm.SOA.sigma_guess/(8.314*298*modelAtm.SOA.rho*1e-6*modelAtm.Pop1.Dp0));
    else
        Kelvin(i) = 1;            
   end
    
   x(i) = BGCp(i)/sum(BGCp);
   BGCv(i) = x(i)*Kelvin(i)*CStarBasis(i);
  % BGCv(i) = 0;
  
   modelAtm.Cv0 = BGCv; 
end