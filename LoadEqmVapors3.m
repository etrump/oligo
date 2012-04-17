function [BGCv] = LoadEqmVapors(BGCp)

global modelAtm

CStarBasis = modelAtm.CStarBasis;
n = modelAtm.NumBins;

for j = 1:modelAtm.Pop
    
    CurPop = int2str(j);
    eval(['Dp = modelAtm.Pop' CurPop '.Dp'])
    if modelAtm.KelvinYN==1
        if CurPop==1
            ShapeFactor = 1;
        else
           %eval(['ShapeFactor = modelAtm.Pop' CurPop '.diamFact']);  
           ShapeFactor = 1.29;
        end
        Kelvin(i) = exp(4*modelAtm.SOA.MW*modelAtm.SOA.sigma_guess/(8.314*298*ShapeFactor*modelAtm.SOA.rho*1e-6*Dp));
    else
        Kelvin(i) = 1;            
    end

for i=1:n
   x(i) = x(i)+BGCp(i)/sum(BGCp);
   BGCv(i) = x(i)*CStarBasis(i);
  % BGCv(i) = 0;
end

end