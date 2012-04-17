figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure1,'LineWidth',3,'FontSize',16);
box('on');
hold('all');

global Tagg DpT modelAtm Tagg_No DpT_No

Tagg_hr = Tagg/3600;
Tagg_hr_No = Tagg_No/3600;
DpT_nm = DpT*1e9;
DpT_nm_No = DpT_No*1e9;

a = length(DpT_nm(:,1))
b = length(DpT_nm(:,2))
c = length(Tagg_hr)


plot1 = plot(Tagg_hr,DpT_nm(:,2),Tagg_hr_No,DpT_nm_No(:,2),'LineWidth',3)
%hold on
%plot1 = plot(Tagg_hr_No,DpT_nm_No(:,2),'Color',[0 1 0],'LineWidth',3)
%hold off

%plot1 = plot(X1,YMatrix1,'Parent',axes1,'LineWidth',3);
%set(plot1(1),'DisplayName','Background','Color',[0 0.6 0]);
%set(plot1(2),'DisplayName','Nanoparticle','Color',[0 0 1]);
set(plot1(1),'DisplayName','Case 1')
set(plot1(1),'Color',[0 0.4 0])
set(plot1(2),'DisplayName','Case 2')
set(plot1(2),'Color',[0 1 0])

xlabel({'Time (hr)'},'FontSize',16);

% Create ylabel
ylabel({'D_{p,ave} (nm)'},'FontSize',16);

% Create legend
paperSize = [5 4];
set(gcf,'Units','inches');
pos = get(gcf,'Position');
set(gcf,'Position',[[pos(1) pos(2)]  paperSize]);
set(gcf,'PaperPosition',[[0 0] paperSize]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', paperSize);
set(gcf,'Color','none');

legend1 = legend(axes1,'show');
axis([0 5  50 90]);
%set(legend1,'Position',[0.5034 0.1988 0.37 0.1601]);
set(legend1,'Location','NorthWest')

% Make sure there is a place for figures and save as a pdf
if ~exist('./figs','dir'); mkdir('./figs'); end;
saveas(gcf,['./figs/DiamCompare.pdf'],'pdf'); 

