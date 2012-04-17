plot(Tagg,CpT(:,1))

plot(Tagg,DpT(:,1))

plot(Tagg,(CpT(:,1))./(CpT(:,1)+CpT(:,9)))

for i = 1:length(Tagg)
    if CpT(i,9) > 0
        X_frac(i,9) = CpT(i,9)./(CpT(i,1)+CpT(i,2)+CpT(i,3)+CpT(i,4)+CpT(i,5)+CpT(i,6)+CpT(i,7)+CpT(i,8)+CpT(i,9)+CpT(i,10));
    else
        X_frac(i,9) = 0;
    end
end
plot(Tagg, CvT(:,9)-X(:,9)*modelAtm.CStarBasis(9))

Yagg(Indy,length(y0))+sum(CpT(Indy,:))+sum(CvT(Indy,:))
Yagg(1,length(y0))*sum(modelAtm.SOA.alphaProd)+sum(CpT(1,:))+sum(CvT(1,:))

CpTot(1,:) = CpT(:,1)+CpT(:,2)+CpT(:,3)+CpT(:,4)+CpT(:,5)+CpT(:,6)+CpT(:,7)+CpT(:,8)+CpT(:,9)+CpT(:,10);

plot(Tagg,(CpT(:,1))./(CvT(:,9)+CpT(:,1)))