function CreateFigMass()

figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure1,'LineWidth',3,'FontSize',16);
box('on');
hold('all');

global Tagg MassforPlot
Tagg_hr = Tagg/3600;

y = [MassforPlot(1,:)' MassforPlot(2,:)' MassforPlot(3,:)' MassforPlot(4,:)'];
area1 = area(Tagg_hr,y,'LineWidth',3);

set(area1(1),'FaceColor',[0 0 0.4],'DisplayName','Wall C*=0.01');
set(area1(2),'FaceColor',[0 0 0.8],'DisplayName','Wall C*=1');
set(area1(3),'FaceColor',[0 0.4 0],'DisplayName','Suspended C*=0.01');
set(area1(4),'FaceColor',[0 0.8 0],'DisplayName','Suspended C*=1');
xlabel({'Time (hr)'},'FontSize',16);
ylabel({'Mass ( \mug  m^{-3})'},'FontSize',16);
%axis([0 10 -1 1]);
legend1 = legend(axes1,'show');
set(legend1,'Position',[0.4585 0.2485 0.4187 0.264]);


%%%%
figure2 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes2 = axes('Parent',figure2,'LineWidth',3,'FontSize',16);
box('on');
hold('all');

global Tagg MassforPlot
Tagg_hr = Tagg/3600;

y = [MassforPlot(1,:)' MassforPlot(3,:)' MassforPlot(2,:)' MassforPlot(4,:)'];
plot2 = area(Tagg_hr,y,'LineWidth',3);

set(plot2(1),'FaceColor',[0 0 0.4],'DisplayName','Wall C*=0.01');
set(plot2(3),'FaceColor',[0 0 0.8],'DisplayName','Wall C*=1');
set(plot2(2),'FaceColor',[0 0.4 0],'DisplayName','Suspended C*=0.01');
set(plot2(4),'FaceColor',[0 0.8 0],'DisplayName','Suspended C*=1');
xlabel({'Time (hr)'},'FontSize',16);
ylabel({'Mass ( \mug  m^{-3})'},'FontSize',16);
%axis([0 10 -1 1]);
legend2 = legend(axes2,'show');
set(legend2,'Position',[0.4585 0.2485 0.4187 0.264]);
