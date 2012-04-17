function ProcessData()

global modelAtm trackCov MassforPlot MassforPlotAlpha Struct
global Tagg CpT DpT NumConcT CvT Yagg_true CvT_true CpT_true
global trackCondSink

Yagg = Struct.Yagg;
Tagg = Struct.Tagg;
%EndIndex = Struct.EndIndex;
EmitTime = modelAtm.EmitTime;



n = modelAtm.NumBins;
TotalPop = modelAtm.Pop;

%trackCondSink = zeros(1,11);

%modelAtm.Dpf = Dpf;
EndIndex(1) = 0;

CpT = zeros(0,n); DpT = zeros(0,1); MSulfT = zeros(0,1); %NumConcT = zeros(0,1);

CvT = Yagg(:,(1:n));
CvT_true = Yagg_true(:,(1:n));
p = modelAtm.Pop;

%for p = 2:TotalPop
   CpT = Yagg(:,n+1:(p+1)*n);
   CpT_true = Yagg_true(:,n+1:(p+1)*n);
   DpT = Yagg(:, n*(p+1)+1:n*(p+1)+p);
   MSulfT = Yagg(:, n*(p+1)+p+1:n*(p+1)+2*p);
   %NumConc_new = Yagg(EndIndex(p-1)+1:EndIndex(p), n*(p+1)+2*p+1:n*(p+1)+3*p);
   
   %CpT = [Cp_new];
   %DpT = [Dp_new];
   %MSulfT = [MSulf_new];
   %NumConcT = [NumConcT zeros(EndIndex(p-1),1); NumConc_new];
%end

% %Check mass balance---------------------------------------------------------------------------
% TotalSuspMSOA = 0;
% TotalPop2MSOA = 0;
% TotalVap = 0;
% 
% for j=1:n
%     TotalPop2MSOA = TotalPop2MSOA + CpT(length(CpT),n+j);
%     TotalSuspMSOA = TotalSuspMSOA + CpT(length(CpT),j);
%     TotalVap = TotalVap + CvT(length(CpT),j);
% end
% 
% TotalMSOA0 = sum(modelAtm.Pop1.Cp0)+sum(modelAtm.Pop2.Cp0);
% Tau1 = modelAtm.DecayConstant*3600;
% 
% if modelAtm.PulseYN==1
%     CheckMassBal = sum(modelAtm.Cv0)+modelAtm.ROG*modelAtm.EmitTime*sum(modelAtm.SOA.alphaProd)+TotalMSOA0-TotalPop2MSOA-TotalSuspMSOA-TotalVap;
% else
%     CheckMassBal = sum(modelAtm.Cv0)+modelAtm.Injection*(1-exp(-FinalTime/Tau1))*sum(modelAtm.SOA.alphaProd)+TotalMSOA0-TotalPop2MSOA-TotalSuspMSOA-TotalVap;
% end
% 
% toler = 1e-7;
% if abs(CheckMassBal) > 1e-7
%     Error = 'Mass balance is not correct!';
% else
%     Error = 'NO mass balance error!';
% end
% 
% CheckEqm(1) = (CpT(length(CpT),1)+CpT(length(CpT),n+1))/(CvT(length(CpT),1)+CpT(length(CpT),1)+CpT(length(CpT),1*n+1))-...
%     1/(1+modelAtm.CStarBasis(1)/(TotalPop2MSOA+TotalSuspMSOA));
% CheckEqm(2) = NaN;
% CheckEqm(3) = (CpT(length(CpT),3)+CpT(length(CpT),n+3))/(CvT(length(CpT),3)+CpT(length(CpT),3)+CpT(length(CpT),1*n+3))-...
%     1/(1+modelAtm.CStarBasis(3)/(TotalPop2MSOA+TotalSuspMSOA));
% 
% FracMassToSuspended = TotalSuspMSOA/(TotalSuspMSOA+TotalPop2MSOA);
% 
% %Calculate mole fractions of organics on particles--------------------------
% for p = 1:TotalPop
%    PopString = int2str(p);
%    XCpT_1 = CpT(:,(p-1)*n+3)./(CpT(:,(p-1)*n+1)+CpT(:,(p-1)*n+3));
%    XCpT_01 = CpT(:,(p-1)*n+1)./(CpT(:,(p-1)*n+1)+CpT(:,(p-1)*n+3));
%   % eval(['Struct.P' PopString '.XCpT_1 = XCpT_1;']);
%    eval(['modelAtm.Pop' PopString '.XCpT_1 = XCpT_1;']);
%   % eval(['Struct.P' PopString '.XCpT_1 = XCpT_1;']);
%    eval(['modelAtm.Pop' PopString '.XCpT_1 = XCpT_1;']);
% end
% 
%  for i=1:length(Tagg)
%      Cp1Wall(i) = 0;
%      Cp01Wall(i) = 0;
%      Cp_1Wall(i) = 0;
%      Cp1_Wall(i) = 0;
%      for j = 2:TotalPop
%         Cp1Wall(i) = Cp1Wall(i) + CpT(i,4*(j-1)+3);
%         Cp_1Wall(i) = Cp_1Wall(i) + CpT(i,4*(j-1)+2);
%         Cp01Wall(i) = Cp01Wall(i) + CpT(i,4*(j-1)+1);
%         Cp1_Wall(i) = Cp1_Wall(i) + CpT(i,4*(j-1)+4);
%      end
%      MSulfWall(i) = MSulfT(i,2);
%      Cp1Sus(i) = CpT(i,3);
%      Cp1_Sus(i) = CpT(i,4);
%      Cp_1Sus(i) = CpT(i,2);
%      Cp01Sus(i) = CpT(i,1);
%      MSulfSus(i) = MSulfT(i,1);
%  end
% 
% MassforPlot = [Cp01Wall; Cp1Wall; Cp01Sus; Cp1Sus; MSulfWall; MSulfSus];
% MassforPlotAlpha = [Cp01Wall; Cp_1Wall; Cp1Wall; Cp1_Wall; Cp01Sus; Cp_1Sus; Cp1Sus; Cp1_Sus; MSulfWall; MSulfSus];
% 
% 
% 
% for j = 1:length(Tagg)
%    TotalOrgMass = sum(CpT(j,1:modelAtm.NumBins));
%    TotalMass(j) = TotalOrgMass+MSulfT(j,1);
%    rho = (TotalOrgMass+MSulfT(j,1))/(MSulfT(j,1)/modelAtm.Sulf.rho+TotalOrgMass/modelAtm.SOA.rho);
%    CondSinkBG(j) = UpdateBackgroundCS(TotalMass(j),DpT(j,1),rho);
%    BG_MtoCS(j) = TotalMass(j)/CondSinkBG(j);
%    
%    
%    TotalOrgMass = sum(CpT(j,5:2*modelAtm.NumBins));
%    TotalMass(j) = TotalOrgMass + MSulfT(j,2);
%    rho = (TotalOrgMass+MSulfT(j,2))/(MSulfT(j,2)/modelAtm.Sulf.rho+TotalOrgMass/modelAtm.SOA.rho);
%    
%    ParticleKn = 2*modelAtm.SOA.lambda/DpT(j,2);
%    BetaP = FuchsC(ParticleKn);
%    CondSinkNanoP(j) = modelAtm.SOA.alpha*modelAtm.Pop2.NumConc0*2*pi*DpT(j,2)*modelAtm.SOA.Diff*BetaP;
%    %CondSinkNanoP(j) = UpdateBackgroundCS(TotalMass(j),DpT(j,2),rho);
%    NanoP_MtoCS(j) = TotalMass(j)/CondSinkNanoP(j);
%    
%    
%    
%    
%    
% end
% 
% for j = 1:length(trackCondSink(:,1))
% totalCondSink = trackCondSink(j,2)+trackCondSink(j,3);
%    BG_PctCondSink(j) = trackCondSink(j,2)/totalCondSink;
%    NP_PctCondSink(j) = trackCondSink(j,3)/totalCondSink;
%    
%    totalFlux = trackCondSink(j,10)+trackCondSink(j,11);
%    BG_PctFlux(j) = trackCondSink(j,10)/totalFlux;
%    NP_PctFlux(j) = trackCondSink(j,11)/totalFlux;
% end
% 
% 
% 
% figure(78)
% plot(trackCondSink(:,1),trackCondSink(:,8),trackCondSink(:,1),trackCondSink(:,9))
% %plot(Tagg/3600,BG_MtoCS,Tagg/3600,NanoP_MtoCS)
% %plot(
% %DISPLAY FIGURES----------------------------------------------------
% %CreateFigNumConc; 
% 
% CreateFig;
% if modelAtm.AlphaPYN == 1 || modelAtm.AlphaPYNPop2 == 1 || modelAtm.AlphaPYNVap == 1
%     CreateFigMassAlpha2;
%     CreateMolFracFigAlpha;
% else
%     CreateFigMass;
%     CreateMolFracFig;
% end
% 
% if modelAtm.EmitBin ==1
%    DispEmitBin = 1e-2;
% elseif modelAtm.EmitBin==3
%    DispEmitBin = 1e0;
% end
% 
% if modelAtm.BGBin ==1
%    DispBGBin = 1e-2;
% elseif modelAtm.BGBin==3
%    DispBGBin = 1e0;
% end


%table.ROG = modelAtm.ROG;
%table.EmitTime = EmitTime/3600;
%table.EmitBin = DispEmitBin;
%table.BGBin = DispBGBin;
%table.TotalSusp = TotalSusp;
%table.TotalOnWall = TotalOnWall;
%table.kwall = modelAtm.kwall;
%table;