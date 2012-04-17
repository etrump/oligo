function CreateFigOmega()

figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure1,'LineWidth',3,'FontSize',16);
box('on');
hold('all');

global Tagg omegaNum modelAtm omegaOverall

Tagg_hr = Tagg/3600;
%plot1 = plot(Tagg_hr,omegaNum(:,2),Tagg_hr,omegaNum(:,3),Tagg_hr,omegaNum(:,4),Tagg_hr,omegaNum(:,5),'LineWidth',3)
plot1 = plot(Tagg_hr,omegaNum(:,2),'LineWidth',3)

%plot1 = plot(X1,YMatrix1,'Parent',axes1,'LineWidth',3);


% Create xlabel
xlabel({'Time (hr)'},'FontSize',16);

% Create ylabel
ylabel({'Omega,numeric'},'FontSize',16);

axis([0 10 -1 1]);

% Create legend
%legend1 = legend(axes1,'show');
%set(legend1,'Position',[0.5638 0.2009 0.2955 0.264]);