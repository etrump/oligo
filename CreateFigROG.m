global modelAtm Tagg

Tau = modelAtm.DecayConstant*3600;
for i = 1:length(Tagg)
   
ROG(i) = modelAtm.Injection/Tau*exp(-Tagg(i)/Tau);
  
end

Tagg_hr = Tagg/3600;

figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure1,'LineWidth',3,'FontSize',16);
box('on');
hold('all');

plot1 = plot(Tagg_hr,ROG,'LineWidth',3);

%plot1 = plot(X1,YMatrix1,'Parent',axes1,'LineWidth',3);
%set(plot1(1),'DisplayName','Suspension');
%set(plot1(2),'DisplayName','Wall');


% Create xlabel
xlabel({'Time (hr)'},'FontSize',16);

% Create ylabe
ylabel({'R_{COC} (\mug m^{-3}s^{-1})'},'FontSize',16);

% Create legend
%legend1 = legend(axes1,'show');