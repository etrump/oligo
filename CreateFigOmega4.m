function CreateFigOmega()

figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
%hold on
axes1 = axes('Parent',figure1,'LineWidth',3,'FontSize',16);
box('on');
hold('all');

global Tagg omegaNum modelAtm omegaOverall trackOmega omegaCheck omegaOverall_inst omegaNum_inst omegaW

trackOmega(1,2) = NaN;
omegaW(1) = NaN;
omegaW(2) = NaN;
omegaW(3) = NaN;

for jj = 1:25
   omegaW(jj) = NaN; 
end
Tagg_hr = Tagg/3600;
%plot1 = plot(Tagg_hr,omegaNum(:,2),Tagg_hr,omegaNum(:,3),Tagg_hr,omegaNum(:,4),Tagg_hr,omegaNum(:,5),'LineWidth',3)
checkzero = zeros(length(Tagg_hr),1);
plot1 = plot(trackOmega(:,1)/3600,trackOmega(:,2),'--',Tagg_hr,omegaW,Tagg_hr,checkzero,'LineWidth',3)
%plot1 = plot(Tagg_hr,omegaCheck(:,2))
%plot2 =plot(Tagg/3600,0)

set(plot1(1),'DisplayName','\omega_{Diffusive}','Color',[0 0 1]);
set(plot1(2),'DisplayName','\omega_{Experimental}','Color',[0 0 1]);
set(plot1(3),'LineWidth',1,'LineStyle','--','Color',[0 0 0]);
legend(plot1,'\omega_{Diffusive}','\omega_{Experimental}')

% Create xlabel
xlabel({'Time (hr)'},'FontSize',16);

% Create ylabel
ylabel({'Omega'},'FontSize',16);

axis([0 10 -1 1]);
%axis([0 3 -1 1]);

% Create legend
%legend1 = legend(axes1,'show');
%set(legend1,'Position',[0.5638 0.2009 0.2955 0.264]);
legend1 = legend(axes1,'show');


if ~exist('./figs','dir'); mkdir('./figs'); end;
saveas(gcf,['./figs/Omega.png'],'png'); 