%function createfigure1(X1, YMatrix1)
%CREATEFIGURE1(X1,YMATRIX1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 12-Mar-2010 13:40:11

global Tagg Struct modelAtm
% Create figure
figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',figure1,'LineWidth',2,'FontSize',16);
% Uncomment the following line to preserve the X-limits of the axes
% xlim([0 3]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim([0.5 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim([-1 1]);
box('on');
hold('all');

for m = 1:9
modelAtm.Pop1.XCpT_1(m) = NaN;
modelAtm.Pop2.XCpT_1(m) = NaN;
end


plot1 = plot(Tagg/3600,modelAtm.Pop1.XCpT_1,Tagg/3600,modelAtm.Pop2.XCpT_1,'Parent',axes1,'LineWidth',2);
set(plot1(1),'DisplayName','Background');
set(plot1(2),'DisplayName','Nanopaticle');

xlabel({'Time (hr)'},'FontSize',16);
ylabel({'X(C*=1)'},'FontSize',16);


%set(legend1,'Position',[0.5206 0.1837 0.3589 0.1601])

% Set the size of the figure on the screen and for printing
% Rediculously complicated but such is life
paperSize = [5 4];
set(gcf,'Units','inches');
pos = get(gcf,'Position');
set(gcf,'Position',[[pos(1) pos(2)]  paperSize]);
set(gcf,'PaperPosition',[[0 0] paperSize]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', paperSize);
set(gcf,'Color','none');
legend1 = legend(axes1,'show');
if modelAtm.BGBin==1
set(legend1,'Position',[0.5139 0.1753 0.3589 0.1601]);
set(legend1,'Location','NorthWest')
else
set(legend1,'Position',[0.5206 0.7364 0.3589 0.1601]);
end
axis([0 5 0 1]);
% Make sure there is a place for figures and save as a pdf
if ~exist('./figs','dir'); mkdir('./figs'); end;
saveas(gcf,['./figs/MolFrac.pdf'],'pdf'); 