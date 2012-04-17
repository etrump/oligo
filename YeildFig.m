
%plot(Tagg,DpT(:,1))
%Pop = modelAtm.Pop;
%index_pre = n*(1+Pop)+2*Pop+1;
%length(Tagg)

%MSOA = zeros(1,length(Tagg));
%for i = 1:length(Tagg)
%    for j = 1:modelAtm.NumBins
%        MSOA(i) = MSOA(i) + CpT(i,j);
%    end
%    
%    Yeild_SOA(i) = MSOA(i)/(Yagg(1,index_pre)-Yagg(i,index_pre));
%end

%figure(56)
%plot(Tagg,MSOA)

%figure(57)
%plot(Tagg,Yeild_SOA)

%figure(58)
%plot(MSOA,Yeild_SOA)
%
%alpha_P_added = modelAtm.Injection
%Ult_yield = Yeild_SOA(length(Tagg))
%Ult_loading = MSOA(length(Tagg))

Cstar = modelAtm.CStarBasis
C_tot = 0.5:0.1:4000;

C_tot = 1:2800;
for i = 1:length(C_tot)
   
    C = C_tot(i) * modelAtm.SOA.alphaProd;
    
    

    maxiter = 1e15;
    tol = 1e-25;
    Coa = 1; %initial guess

 for n = 1:maxiter
    xi = (1 + Cstar/Coa ) .^ (-1);
    CoaNEW = sum(C .* xi);
    if abs(CoaNEW - Coa) < tol
     %   finalDiff = CoaNEW-Coa
        break
    else
        diffNOW = CoaNEW-Coa;
        Coa = CoaNEW;
    end
 end
 C_SOA(i) = Coa;
 Coa;
 sum(C);
 AMF(i) = Coa/sum(C);
 Yeild_SOA(i) = Coa/C_tot(i);    
end

figure(99)
%plot(C_SOA,AMF)
%plot(C_SOA,Yeild_SOA)
semilogx(C_SOA,Yeild_SOA)
hold on
axis([0.01 1e4 0 0.4]);

C_dyn = [10.5037 274.5 471.38 933.4]
Y_dyn = [0.105 0.2745 0.3143 0.3734]

%figure(8)
semilogx(C_dyn,Y_dyn)