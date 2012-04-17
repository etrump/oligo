function CreateFigOmega()

figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure1,'LineWidth',3,'FontSize',16);
box('on');
hold('all');

global Tagg omegaNum modelAtm omegaOverall trackOmega omegaCheck omegaOverall_inst omegaNum_inst omegaW

Tagg_hr = Tagg/3600;
%plot1 = plot(Tagg_hr,omegaNum(:,2),Tagg_hr,omegaNum(:,3),Tagg_hr,omegaNum(:,4),Tagg_hr,omegaNum(:,5),'LineWidth',3)
plot1 = plot(Tagg_hr,omegaNum(:,2),trackOmega(:,1)/3600,trackOmega(:,2),Tagg_hr,omegaOverall',...
    Tagg_hr,omegaOverall_inst',Tagg_hr,omegaNum_inst(:,2),Tagg_hr,omegaW,'LineWidth',3)
%plot1 = plot(Tagg_hr,omegaCheck(:,2))


set(plot1(1),'DisplayName','Numeric, dDp3/dDp3susp');
set(plot1(2),'DisplayName','Diffusive limit');
set(plot1(3),'DisplayName','Numeric, dM/dMsusp');
set(plot1(4),'DisplayName','dM/dM_susp inst');
set(plot1(5),'DisplayName','dDp3/dDp3susp inst');
set(plot1(6),'DisplayName','OMEGA');

% Create xlabel
xlabel({'Time (hr)'},'FontSize',16);

% Create ylabel
ylabel({'Omega'},'FontSize',16);

axis([0 10 -1 2]);

% Create legend
%legend1 = legend(axes1,'show');
%set(legend1,'Position',[0.5638 0.2009 0.2955 0.264]);
legend1 = legend(axes1,'show');