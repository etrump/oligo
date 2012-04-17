function CreateFigMass()

figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure1,'LineWidth',3,'FontSize',16);
box('on');
hold('all');

global Tagg MassforPlot NumConcT
global modelAtm
Tagg_hr = Tagg/3600;


%y = [MassforPlot(5,:)' MassforPlot(1,:)' MassforPlot(2,:)'];
y = [MassforPlot(5,:)' MassforPlot(2,:)' MassforPlot(1,:)'];
area1 = area(Tagg_hr,y,'LineWidth',3);

%set(area1(2),'FaceColor',[0 0.4 0],'DisplayName','lC*=-2');
%set(area1(3),'FaceColor',[0 1 0],'DisplayName','lC*=  0');
set(area1(1),'FaceColor',[1 0 0],'DisplayName','Sulfate');
set(area1(3),'FaceColor',[0 0.4 0],'DisplayName','lC*=-2');
set(area1(2),'FaceColor',[0 1 0],'DisplayName','lC*=  0');

xlabel({'Time (hr)'},'FontSize',16);
ylabel({'Mass ( \mug  m^{-3})'},'FontSize',16);
%title('Nanoparticle')
%title('Aged Background') 

%axis([0 10 -1 1]);

paperSize = [5 4];
set(gcf,'Units','inches');
pos = get(gcf,'Position');
set(gcf,'Position',[[pos(1) pos(2)]  paperSize]);
set(gcf,'PaperPosition',[[0 0] paperSize]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', paperSize);
set(gcf,'Color','none');
legend1 = legend(axes1,'show');


set(legend1,'Location','NorthWest')
axis([0 5 0 5e-10]);
axis([0 5 0 0.6e-3]); %
axis([0 5 0 1e-3]);
% Make sure there is a place for figures and save as a pdf
if ~exist('./figs','dir'); mkdir('./figs'); end;
saveas(gcf,['./figs/MassNanoP.pdf'],'pdf'); 







% %%%%

figure2 = figure('InvertHardcopy','off','Color',[1 1 1]);
 axes2 = axes('Parent',figure2,'LineWidth',3,'FontSize',16);
 box('on');
 hold('all'); 
 %global Tagg MassforPlot
 %Tagg_hr = Tagg/3600; 
% y = [MassforPlot(6,:)' MassforPlot(3,:)' MassforPlot(4,:)']
 y = [MassforPlot(6,:)' MassforPlot(4,:)' MassforPlot(3,:)']
 plot2 = area(Tagg_hr,y,'LineWidth',3);
% 
 set(plot2(1),'FaceColor',[1 0 0],'DisplayName','Sulfate');
 %set(plot2(3),'FaceColor',[0 1 0],'DisplayName','lC*=  0');
 %set(plot2(2),'FaceColor',[0 0.4 0],'DisplayName','lC*=-2');
 set(plot2(2),'FaceColor',[0 1 0],'DisplayName','lC*=  0');
 set(plot2(3),'FaceColor',[0 0.4 0],'DisplayName','lC*=-2');
 xlabel({'Time (hr)'},'FontSize',16);
 ylabel({'Mass ( \mug  m^{-3})'},'FontSize',16);
 %title('Background')

 
paperSize = [5 4];
set(gcf,'Units','inches');
pos = get(gcf,'Position');
set(gcf,'Position',[[pos(1) pos(2)]  paperSize]);
set(gcf,'PaperPosition',[[0 0] paperSize]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', paperSize);
set(gcf,'Color','none');
axis([0 5 0 1.5]);
%axis([0 5 0 1.5e3]);

legend2 = legend(axes2,'show');
%set(legend2,'Position',[0.1687 0.6124 0.3856 0.264]);
%set(legend2,'Position',[0.1687 0.6124 0.2611 0.264]);
set(legend2,'Location','NorthWest','Orientation','horizontal')

% Make sure there is a place for figures and save as a pdf
if ~exist('./figs','dir'); mkdir('./figs'); end;
saveas(gcf,['./figs/MassBG.pdf'],'pdf'); 
