function CreateFigS()

figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure1,'LineWidth',3,'FontSize',16);
box('on');
hold('all');

global Tagg CvT CpT modelAtm
Tagg_hr = Tagg/3600;

MolFrac(:,1) = CvT(:,1)./(CvT(:,1)+CvT(:,2)+CvT(:,3)+CvT(:,4)); 
MolFrac(:,3) = CvT(:,3)./(CvT(:,1)+CvT(:,2)+CvT(:,3)+CvT(:,4)); 

MolFrac(:,1) = CpT(:,1)./(CpT(:,1)+CpT(:,2)+CpT(:,3)+CpT(:,4))
MolFrac(:,3) = CpT(:,3)./(CpT(:,1)+CpT(:,2)+CpT(:,3)+CpT(:,4))

%plot1 = semilogy(Tagg_hr,CvT(:,1)./(modelAtm.CStarBasis(1).*MolFrac(:,1)),Tagg_hr,CvT(:,3)./(modelAtm.CStarBasis(3).*MolFrac(:,3)));

plot1 = semilogy(Tagg_hr,CvT(:,1)./(modelAtm.CStarBasis(1)),Tagg_hr,CvT(:,3)./(modelAtm.CStarBasis(3)));

set(plot1(1),'Color',[0 0 0],'LineWidth',3,'DisplayName','S, lC*=-2');
set(plot1(2),'Color',[0 0 0],'LineWidth',3,'DisplayName','S, lC*=0','LineStyle','--');

xlabel({'Time (hr)'},'FontSize',16);
ylabel({'Saturation Ratio (C_{v,i}/C*_{i})'},'FontSize',16);




paperSize = [5 4];
set(gcf,'Units','inches');
pos = get(gcf,'Position');
set(gcf,'Position',[[pos(1) pos(2)]  paperSize]);
set(gcf,'PaperPosition',[[0 0] paperSize]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', paperSize);
set(gcf,'Color','none');

legend1 = legend(axes1,'show');
set(legend1,'Location','NorthEast');


axis([0 5 0 2]);
% Make sure there is a place for figures and save as a pdf
if ~exist('./figs','dir'); mkdir('./figs'); end;
saveas(gcf,['./figs/Sat.pdf'],'pdf'); 