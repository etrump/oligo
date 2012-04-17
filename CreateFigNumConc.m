figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure1,'LineWidth',3,'FontSize',16);
box('on');
hold('all');

global Tagg NumConcT modelAtm 

Tagg_hr = Tagg/3600;
DpT_nm = DpT*1e9;
if modelAtm.Pop>3
    plot1 = plot(Tagg_hr,NumConcT(:,1),Tagg_hr,NumConcT(:,2),Tagg_hr,NumConcT(:,3),Tagg_hr,NumConcT(:,4),'LineWidth',3)
else
    plot1 = plot(Tagg_hr,NumConcT(:,1),Tagg_hr,NumConcT(:,2),'LineWidth',3)
end
%plot1 = plot(X1,YMatrix1,'Parent',axes1,'LineWidth',3);
set(plot1(1),'DisplayName','Suspension');
set(plot1(2),'DisplayName','Wall');
if modelAtm.Pop>3
set(plot1(3),'DisplayName','Pop 3');
set(plot1(4),'DisplayName','Pop 4');
end

% Create xlabel
xlabel({'Time (hr)'},'FontSize',16);

% Create ylabel
ylabel({'Number Concentration (m^{-3})'},'FontSize',16);

% Create legend
legend1 = legend(axes1,'show');
%set(legend1,'Position',[0.5638 0.2009 0.2955 0.264]);