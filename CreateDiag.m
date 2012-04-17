
global Tagg DpT modelAtm CpT_true CvT_true
dimerbin = 8;

CpTot = zeros(length(Tagg));
CvTot = zeros(length(Tagg));
X_frac_part = zeros(length(Tagg),modelAtm.NumBins);
X_frac_vap = zeros(length(Tagg),modelAtm.NumBins);


for i = 1:length(Tagg)
    
%CpTot(i) = CpT(i,1)+CpT(i,2)+CpT(i,3)+CpT(i,4)+CpT(i,5)+CpT(i,6)+CpT(i,7)+CpT(i,8)+CpT(i,9)+CpT(i,10);
CpTot(i) = sum(CpT_true(i,:));
CpTot_Mass(i) = sum(CpT(i,:))+CpT(i,1);
CvTot(i) = CvT(i,1)+CvT(i,2)+CvT(i,3)+CvT(i,4)+CvT(i,5)+CvT(i,6)+CvT(i,7)+CvT(i,8)+CvT(i,9)+CvT(i,10);

if CpTot(i) > eps
X_frac_part(i,1) = CpT_true(i,1)./CpTot(i);
X_frac_part(i,9) = CpT_true(i,9)./CpTot(i);
X_frac_part(i,dimerbin) = CpT_true(i,dimerbin)./CpTot(i);
else
X_frac_part(i,1) = 0;
X_frac_part(i,9) = 0;
X_frac_part(i,8) = 0;    
end

%X_frac_vap(i,1) = CvT(i,1)./CvTot(i);
%X_frac_vap(i,8) = CvT(i,8)./CvTot(i);
%X_frac_vap(i,9) = CvT(i,9)./CvTot(i);

end

DistEqm(:,1) = CvT_true(:,1)-X_frac_part(:,1)*modelAtm.CStarBasis(1);
DistEqm(:,8) = CvT_true(:,8)-X_frac_part(:,8)*modelAtm.CStarBasis(8);
DistEqm(:,dimerbin) = CvT_true(:,dimerbin)-X_frac_part(:,dimerbin)*modelAtm.CStarBasis(dimerbin);

figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure1,'LineWidth',3,'FontSize',16);
box('on');
hold('all');

DistEqm(length(Tagg),dimerbin);

plot1 = plot(Tagg/3600, DistEqm(:,dimerbin),'LineWidth',3)

xlabel({'Time (hr)'},'FontSize',16);
ylabel({'Cv-XC*(8)'},'FontSize',16);
%title({'k_f = 0.75'},'FontSize',16)

Olig_Frac = CpT(:,1)./(CpTot_Mass(:));

figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure1,'LineWidth',3,'FontSize',16);
box('on');
hold('all');

plot2 = plot(Tagg/3600,Olig_Frac,'LineWidth',3)

