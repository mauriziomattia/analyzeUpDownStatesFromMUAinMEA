% analyzeUpDownStateVsRecordingSet.m
%
%   Copyright 2015 Maurizio Mattia @ Ist. Super. Sanità, Rome - Italy
%   Version: 1.0 - Feb. 18, 2015
%


%% Load analysis params and information about recording dataset.
%
setParamsAndOptions


%% Open data file...
%
hSMR = openSONandLoadInfo(DataFile);
if isempty(hSMR)
   return
end


%% Analysis loop on recording set...
%
try
for crsRecSet = ChannelSet

   %% Makes the working directory and change the current dir to it...
   %
   tempdir = [AnalysisDir RecordingSet(crsRecSet).label];
   if exist(tempdir,'dir') == 0
      mkdir(tempdir);
   end
   cd(tempdir);
   disp(['Recording set: ' RecordingSet(crsRecSet).label]);

   %% Loads dataset...
   %
   RawSig = getSONAnalogData(hSMR, RecordingSet(crsRecSet).port, Options.PeriodToAnalyze);
   RawSig.dt = diff(RawSig.time(1:2));
   fprintf('Analyzed period: [%g,%g] s\n', RawSig.time(1), RawSig.time(end));

   clear ndx
   close all
   
   
   %%  Compute baseline for spectral MUA estimate...
   %
   MovingWindowSize = 1/Options.LogMUA.FreqBand(1);
   DetrendingOrder = 0;
   [PyyFreq, PyyBaseline] = plotMedianPSDofLFP(RawSig, MovingWindowSize, DetrendingOrder);
%    ndx = find(PyyFreq>=Options.LogMUA.FreqBand(1)-0.1 & PyyFreq<=Options.LogMUA.FreqBand(2)+0.1);
   ndx = find(PyyFreq>=Options.LogMUA.FreqBand(1)*0.99 & PyyFreq<=Options.LogMUA.FreqBand(2)*1.01);
   PyyFreq = PyyFreq(ndx);
   PyyBaseline = PyyBaseline(ndx);
   
   
   %% Compute spectral estimate of MUA...
   %
   [MUA.time, MUA.value] = computeSpectralEstimateOfMUA(RawSig.time, RawSig.value, ...
         Options.LogMUA.FreqBand, 0, PyyBaseline);
   MUA.dt = diff(MUA.time(1:2));

   
   %% Compute smoothed log(MUA): rsMUA...
   %
   logMUA.value = log(MUA.value);
   logMUA.time = MUA.time;
   logMUA.dt = MUA.dt;
   
   if Options.LogMUA.SmoothingWindow > 0
      [rsMUA.time, rsMUA.value] = computeMovingAverage(logMUA, ...
         round(Options.LogMUA.SmoothingWindow / logMUA.dt), 1);
   else
      rsMUA = logMUA;
   end
   rsMUA.dt = diff(rsMUA.time(1:2));
   
   
   %% Subsample LFP...
   %
   WindowSize = round(1/RawSig.dt/Options.LFP.SubsamplingRate);
   [ssLFP.time,ssLFP.value] = computeMovingAverage(RawSig, WindowSize, WindowSize);
   ssLFP.dt = diff(ssLFP.time(1:2));
   

   %% High-pass filtering of raw signal (LFP), if required...
   %
   if isfield(RecordingSet(crsRecSet),'HPFCutOffFreq')
      if ~isempty(RecordingSet(crsRecSet).HPFCutOffFreq)
         disp('LFP high-pass filtering...')
         [lpLFP.time,lpLFP.value] = computeMovingAverage(ssLFP, ...
            round(1/RecordingSet(crsRecSet).HPFCutOffFreq/ssLFP.dt), ...
            round(1/RecordingSet(crsRecSet).HPFCutOffFreq/ssLFP.dt/10));
%          ssLFP.value = ssLFP.value - interp1(lpLFP.time,lpLFP.value,ssLFP.time,'cubic');
         ssLFP.value = ssLFP.value - interp1(lpLFP.time,lpLFP.value,ssLFP.time,'pchip');
         clear lpLFP
      end
   end
   
   %% Low-pass filtering of LFP...
   %
   if Options.LFP.SmoothingWindow > 0
      [ssLFP.time,ssLFP.value] = computeMovingAverage(ssLFP, ...
         round(Options.LFP.SmoothingWindow/ssLFP.dt), 1);
   end

   
   %% Remove bad periods...
   %
   if isfield(RecordingSet(crsRecSet), 'PeriodToRemove')
      if numel(RecordingSet(crsRecSet).PeriodToRemove) == 2
         ndx = find(rsMUA.time < RecordingSet(crsRecSet).PeriodToRemove(1) | ...
            rsMUA.time > RecordingSet(crsRecSet).PeriodToRemove(2));
         rsMUAbackup = rsMUA;
         rsMUA.time = rsMUA.time(ndx);
         rsMUA.value = rsMUA.value(ndx);
      end
   end
   
   
   %% Computes the bimodal distribution of resampled MUA...
   %
   ModeParams = plotBimodalHistogram(rsMUA.value);
   close(gcf);
   
   
   %% Restore bad periods...
   %
   if isfield(RecordingSet(crsRecSet), 'PeriodToRemove')
      if numel(RecordingSet(crsRecSet).PeriodToRemove) == 2
         rsMUA = rsMUAbackup;
         clear rsMUAbackup
      end
   end
   
   
   %% Shift logMUA to have the first peak in 0...
   %  NOTE: this is not the final reference value for log(MUA)
   %  corresponding to no firing activity... 
   %
   if isfield(ModeParams,'Mu1')
      LogMUAReference = ModeParams.Mu1;
   else
      LogMUAReference = ModeParams.Mu2;
   end
   RecordingSet(crsRecSet).LogMUAReference = LogMUAReference;
   rsMUA.value = rsMUA.value - RecordingSet(crsRecSet).LogMUAReference;
   ModeParams = plotBimodalHistogram(rsMUA.value);
   
   xlabel('log(MUA)');
   set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.5 3.5 5 4]);
   print('-deps2c', ['rsLogMUA.bihist.SW_' ...
      num2str(Options.LogMUA.SmoothingWindow*1000) '.eps']);


   %% PCA to extract the signal to single out Up-Down states...
   %
   PCA_BASED_UPDOWNDETECTION = 0;
   if isfield(RecordingSet(crsRecSet),'PCAforUpDownDetection')
      if ~isempty(RecordingSet(crsRecSet).PCAforUpDownDetection)
         PCA_BASED_UPDOWNDETECTION = RecordingSet(crsRecSet).PCAforUpDownDetection;
      end
   end
   if PCA_BASED_UPDOWNDETECTION
      Y = [rsMUA.value' ssLFP.value'];
      CovMat = cov(Y);
      [pc,variances,explained] = pcacov(CovMat);
      
      MUALFP.value = (Y * pc(:,1))';
      if pc(1,1) < 0
         MUALFP.value = -1.*MUALFP.value;
      end
      MUALFP.time = rsMUA.time;
      MUALFP.dt = rsMUA.dt;
      
      ModeParams = plotBimodalHistogram(MUALFP.value);
      xlabel('PC_1 (MUA,LFP)');
      set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.5 3.5 5 4]);
      print('-deps2c', ['MUALFP.bihist.SW_' ...
         num2str(Options.LogMUA.SmoothingWindow*1000) '.eps']);
      
      MUAShift = mean(rsMUA.value(abs(MUALFP.value-ModeParams.Mu1)<=std(MUALFP.value)/100));
      rsMUA.value = rsMUA.value - MUAShift;
      LogMUAReference = LogMUAReference + MUAShift;
   end

   
   %%  Set parameters for thresholding of resampled MUA...
   %
   if isfield(RecordingSet(crsRecSet), 'AbsoluteThreshold')
      if isempty(RecordingSet(crsRecSet).AbsoluteThreshold)
         if isfield(ModeParams,'Mu1') && isfield(ModeParams,'Mu2')
            UD_THRESHOLD = Options.UD.ThresholdModulation * ...
               (ModeParams.Mu2-ModeParams.Mu1)+ ModeParams.Mu1;
         else
            UD_THRESHOLD = NaN;
         end
      else
         UD_THRESHOLD = RecordingSet(crsRecSet).AbsoluteThreshold;
      end
   else
      if isfield(ModeParams,'Mu1') && isfield(ModeParams,'Mu2')
         UD_THRESHOLD = Options.UD.ThresholdModulation * ...
            (ModeParams.Mu2-ModeParams.Mu1)+ ModeParams.Mu1;
      else
         UD_THRESHOLD = NaN;
      end
   end
   if isfield(RecordingSet(crsRecSet), 'MinStateDuration')
      if ~isempty(RecordingSet(crsRecSet).MinStateDuration)
         MIN_STATE_LEN = RecordingSet(crsRecSet).MinStateDuration;
      end
   else
      if isfield(RecordingSet(crsRecSet), 'UDcycleFraction')
         if ~isempty(RecordingSet(crsRecSet).UDcycleFraction)
            MIN_STATE_LEN = UDcycle * RecordingSet(crsRecSet).UDcycleFraction;
         else
            MIN_STATE_LEN = 0.0;
         end
      else
         MIN_STATE_LEN = 0.0;
      end
   end

   
   %% The following part of is perfomed only if a threshold for MUA Up and
   %  Down state detection has been set..
   %
   LogMUAbaseline = NaN;
   UpState = [];
   DownState = [];
   
   if ~isnan(UD_THRESHOLD)
      disp('Up and Down states detectable...')

      %%  Set parameters for thresholding of resampled MUA...
      %
      if PCA_BASED_UPDOWNDETECTION
         UD = computeBinaryUDState(MUALFP, UD_THRESHOLD, MIN_STATE_LEN);
      else
         UD = computeBinaryUDState(rsMUA, UD_THRESHOLD, MIN_STATE_LEN);
      end
      
      
      %% Single out and adjust transitions times...
      %
      if PCA_BASED_UPDOWNDETECTION
         [Trans, UD] = findAndAdjustUDTransition(MUALFP, UD, UD_THRESHOLD);
      else
         [Trans, UD] = findAndAdjustUDTransition(rsMUA, UD, UD_THRESHOLD);
      end
      
      
      %% Save if needed
      %
      if isfield(Options, 'SaveMUA')
         if Options.SaveMUA
%             save('UD.mat', 'UD', 'Trans');
            save('MUA.mat', 'rsMUA', 'Trans');
         end
      end
      
      %% Plots together a sample of logMUA and RawSig...
      %
      PERIOD_TO_PLOT = [-20 0] + RawSig.time(end);
      
      plotLFPlogMUAandUD(RawSig, rsMUA, UD, PERIOD_TO_PLOT);
      
      
      %%  Histograms of Up and Down state duration...
      %
      if isfield(RecordingSet(crsRecSet), 'PeriodToRemove')
         if numel(RecordingSet(crsRecSet).PeriodToRemove) == 2
            [UpStateLen, DownStateLen] = plotStateDurationHist(Trans, ...
               RecordingSet(crsRecSet).PeriodToRemove);
         else
            [UpStateLen, DownStateLen] = plotStateDurationHist(Trans);
         end
      else
         [UpStateLen, DownStateLen] = plotStateDurationHist(Trans);
      end
      
      DownState.mean_duration = mean(DownStateLen);
      DownState.std_duration = std(DownStateLen);
      DownState.numel_duration = length(DownStateLen);
      UpState.mean_duration = mean(UpStateLen);
      UpState.std_duration = std(UpStateLen);
      UpState.numel_duration = length(UpStateLen);
      
      save('UpDownStateDuration', 'UpStateLen', 'DownStateLen');
      
      disp(['<UD cycle> = ' num2str(mean(UpStateLen)+mean(DownStateLen),3) ...
         's (st.err.: ' num2str(std(UpStateLen)/sqrt(numel(UpStateLen)) ...
         + std(DownStateLen)/sqrt(numel(DownStateLen)),2) 's)']);
      
      
      %%  Rasters and PSTH of log(MUA) and UD around UPWARD transitions...
      %
      WF_TIME_RANGE = [-0.5 1.5]; % Time range around the trigger event.
      
      % Extract useful upward transitions times, the trigger for raster and PSTH...
      ndx = Trans.val == 1;
      Triggers = Trans.time(ndx);
      SortStateLen = UpStateLen;
      Triggers = Triggers(1:numel(SortStateLen));
      ndx = find(Triggers - UD.time(1) > -WF_TIME_RANGE(1) & ...
         UD.time(end) - Triggers > WF_TIME_RANGE(2));
      Triggers = Triggers(ndx);
      SortStateLen = SortStateLen(ndx);
      
      % Remove bad waveforms...
      if isfield(RecordingSet(crsRecSet), 'WFsToRemove')
         if numel(RecordingSet(crsRecSet).WFsToRemove) > 0
            ndx = 1:numel(Triggers);
            ndx = setdiff(ndx, RecordingSet(crsRecSet).WFsToRemove);
            Triggers = Triggers(ndx);
         end
      end
      save('UpwardTransitions.mat', 'Triggers');
      
      %  Rasters of UD...
      WFs = extractEventCenteredWaveForms(UD.time, UD.value, ...
         Triggers, WF_TIME_RANGE);
      
      plotRasterOfWFs(WFs, [0.0 1.0], 'Up/Down');
      delete('RasterOfWFs.eps')
      
      colormap((1-gray/2));
      YLim = get(gca, 'YLim');
      hold on
      plot([0 0], YLim, 'r-');
      xlabel('Time (s)');
      ylabel('Down-Up Transition');
      
      set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]);
      print('-deps2c', 'UD.Raster.eps');
      
      
      %%  PSTH of log(MUA)...
      %
      WFs = extractEventCenteredWaveForms(rsMUA.time, rsMUA.value, ...
         Triggers, WF_TIME_RANGE);
      
      [tWF, MeanWFs, StdWFs] = plotAverageWFs(WFs, 5);
      
      set(gca, 'TickDir', 'out', 'Box', 'on', 'Layer', 'top', ...
         'XLim', WF_TIME_RANGE);
      
      xlabel('Time (s)');
      ylabel('log(MUA)');
      
      set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [0.25 3.5 7.5 4]);
      print('-deps2c', 'LogMUA.Average.eps')
      
      SampleSize = length(WFs);
      save('LogMUA.Average.mat', 'tWF', 'MeanWFs', 'StdWFs', 'SampleSize');
      
      
      %% Average Log(MUA) around upward transitions of pure sequencies DOWN-UP...
      %
      [t, MeanY, StdY] = ...
         plotAverageUpwardTransition(rsMUA, UD, Triggers, WF_TIME_RANGE);
      
      save('LogMUAUpTrans.Average.mat', 't', 'MeanY', 'StdY');
      
      
      %% Data for correlation plot between Up state duration and activity level...
      %
      AVERAGE_WINDOW = [0.0 0.160];
      X = UpStateLen;
      if isfield(RecordingSet(crsRecSet), 'WFsToRemove')
         if numel(RecordingSet(crsRecSet).WFsToRemove) > 0
            ndx = 1:numel(X);
            ndx = setdiff(ndx, RecordingSet(crsRecSet).WFsToRemove);
            X = X(ndx);
         end
      end
      X = X(1:numel(WFs));
      Y = zeros(size(X));
      for k = 1:length(WFs)
         ndx = find(WFs(k).time>AVERAGE_WINDOW(1) & ...
            WFs(k).time<=AVERAGE_WINDOW(2));
         Y(k) = max(WFs(k).value(ndx));
      end
      
      
      %%  Rasters of log(MUA)...
      %
      plotRasterOfWFs(WFs, prctile(rsMUA.value,[2.5 100-2.5]), 'log(MUA)');
      [RasterUpward.t, RasterUpward.Ys] = computeRasterOfWFs(WFs);
      delete('RasterOfWFs.eps')
      
      xlabel('Time (s)');
      ylabel('Down-Up Transition');
      set(gca,'TickDir','out');
      
      set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]);
      print('-deps2c', 'LogMUA.Raster.eps');
      
      
      %%  Rasters of log(MUA) sorted by Up state length...
      %
      [val,ndx] = sort(SortStateLen);
      WFs = extractEventCenteredWaveForms(rsMUA.time, rsMUA.value, ...
         Triggers(ndx), WF_TIME_RANGE);
      plotRasterOfWFs(WFs, prctile(rsMUA.value,[2.5 100-2.5]), 'log(MUA)');
      delete('RasterOfWFs.eps')
      
      xlabel('Time (s)');
      ylabel('Down-Up Transition');
      
      set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]);
      print('-deps2c', 'LogMUA_Sorted.Raster.eps');
      
      
      %%  Rasters and PSTH of log(MUA) around DOWNWARD transitions...
      %
      WF_TIME_RANGE = [-0.5 1.5]; % Time range around the trigger event.
      
      % Extract useful downward transitions times, the trigger for raster and PSTH...
      Triggers = Trans.time(Trans.val == 0);
      SortStateLen = DownStateLen;
      Triggers = Triggers(1:numel(SortStateLen));
      ndx = find(Triggers - UD.time(1) > -WF_TIME_RANGE(1) & ...
         UD.time(end) - Triggers > WF_TIME_RANGE(2));
      Triggers = Triggers(ndx);
      SortStateLen = SortStateLen(ndx);
      
      % Remove bad waveforms...
      if isfield(RecordingSet(crsRecSet), 'WFsToRemove')
         if numel(RecordingSet(crsRecSet).WFsToRemove) > 0
            ndx = 1:numel(Triggers);
            ndx = setdiff(ndx, RecordingSet(crsRecSet).WFsToRemove);
            Triggers = Triggers(ndx);
         end
      end
      save('DownwardTransitions.mat', 'Triggers');
      
      %  PSTH of log(MUA)...
      WFs = extractEventCenteredWaveForms(rsMUA.time, rsMUA.value, ...
         Triggers, WF_TIME_RANGE);
      
      [tWF, MeanWFs, StdWFs] = plotAverageWFs(WFs, 5);
      
      set(gca, 'TickDir', 'out', 'Box', 'off', 'Layer', 'top', ...
         'XLim', WF_TIME_RANGE);
      
      xlabel('Time (s)');
      ylabel('log(MUA)');
      
      set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [0.25 3.5 7.5 4]);
      print('-deps2c', 'LogMUA.Downward.Average.eps')
      
      SampleSize = length(WFs);
      save('LogMUA.Downward.Average.mat', 'tWF', 'MeanWFs', 'StdWFs', 'SampleSize');
      
      LogMUAbaseline = min(MeanWFs(tWF > 0));
      
      %  Rasters of log(MUA)...
      %
      plotRasterOfWFs(WFs, prctile(rsMUA.value,[2.5 100-2.5]), 'log(MUA)');
      %       plotRasterOfWFs(WFs, [-0.75 3.75], 'log(MUA)');
      [RasterDownward.t, RasterDownward.Ys] = computeRasterOfWFs(WFs);
      delete('RasterOfWFs.eps')
      
      xlabel('Time (s)');
      ylabel('Up-Down Transition');
      
      set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]);
      print('-deps2c', 'LogMUA.Downward.Raster.eps');
      
      
      %%  Rasters of log(MUA) sorted by Down state duration...
      %
      [val,ndx] = sort(SortStateLen);
      WFs = extractEventCenteredWaveForms(rsMUA.time, rsMUA.value, ...
         Triggers(ndx), WF_TIME_RANGE);
      %    T = Triggers(ndx);
      %    WFs = extractEventCenteredWaveForms(rsMUA.time, rsMUA.value, ...
      %             T([1:37 39:end]), [-1 10]);
      plotRasterOfWFs(WFs, prctile(rsMUA.value,[2.5 100-2.5]), 'log(MUA)');
      delete('RasterOfWFs.eps')
      
      xlabel('Time (s)');
      ylabel('Up-Down Transition');
      
      set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [1.0 3.5 7.0 4]);
      print('-deps2c', 'LogMUA_Sorted.Downward.Raster.eps');
      
      
      %% Average Log(MUA) around downward transitions...
      %
      [t, MeanY, StdY] = ...
         plotAverageDownwardTransition(rsMUA, UD, Triggers, WF_TIME_RANGE);
      
      save('LogMUADownTrans.Average.mat', 't', 'MeanY', 'StdY');
      
      
      %% Correlation plot between Up state duration and activity level...
      %
      Y = Y - LogMUAbaseline;
      [RCC, pCC] = corrcoef(X, Y);
      
      figure
      hold on
      plot(X*1000., Y, 'bo');
      
      % Linear fit...
      XRange = get(gca, 'XLim')/1000.;
      pl = polyfit(X, Y, 1);
      plX = linspace(XRange(1), XRange(2), 20);
      plY = polyval(pl, plX);
      
      plot(plX*1000, plY, 'r--');
      title(sprintf('R = %.2f, p = %.3f, n = %d', RCC(1,2), pCC(1,2), length(X)));
      xlabel('Up state duration (ms)');
      ylabel('Up state activity');
      set(gca,'Layer', 'top', 'TickDir', 'out', 'Box', 'off');
      set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [2 3.5 4 4]);
      print('-deps2c', 'UpDurationVsActivity.eps');
      
      
   end % if ~isnan(UD_THRESHOLD)
   
   
   %% Plot the flow and density of LFP-vs-log(MUA)...
   %
   MinDensity = 0.025/100; % Threshold value where to plot the flow.
   GRID_SIZE = [30 30];
   for FIXED_PLOT_RANGE = 0:1
      
      X = rsMUA.value;
      if FIXED_PLOT_RANGE
         Y = ssLFP.value;
         X_RANGE = [-0.75 2.25]; % log(MUA)
         if isfield(Options.LFPMUAplot,'MUArange')
            X_RANGE = Options.LFPMUAplot.MUArange;
         end
         Y_RANGE = [-2500 3500]; % LFP
         if isfield(Options.LFPMUAplot,'LFPrange')
            Y_RANGE = Options.LFPMUAplot.LFPrange;
         end
      else
         Y = ssLFP.value/std(ssLFP.value); % Here we use a normalized measure.
%          X_RANGE = [min(rsMUA.value) max(rsMUA.value)];
         X_RANGE = prctile(rsMUA.value,[0.005 100-0.02]);
         w = 2.1;
         q = prctile(Y,[25 50 75]);
         Y_RANGE = [q(1)-w*(q(3)-q(1)) q(3)+w*(q(3)-q(1))];
         Y_RANGE = [-1 1]*max(abs(Y_RANGE));
%          Y_RANGE = prctile(LFP,[0.05 100-0.05]);
      end
%       Z_RANGE = [0 1.]; % Density
      
      Density = zeros(GRID_SIZE);
      XFlow = zeros(GRID_SIZE);
      YFlow = zeros(GRID_SIZE);
      
      X = X(Y>=Y_RANGE(1) & Y<=Y_RANGE(2));
      Y = Y(Y>=Y_RANGE(1) & Y<=Y_RANGE(2));
      Y = Y(X>=X_RANGE(1) & X<=X_RANGE(2));
      X = X(X>=X_RANGE(1) & X<=X_RANGE(2));
      
      nx = floor((X - X_RANGE(1))/diff(X_RANGE)*(GRID_SIZE(1)-1))+1;
      ny = floor((Y - Y_RANGE(1))/diff(Y_RANGE)*(GRID_SIZE(2)-1))+1;
      for k = 1:numel(nx)-1
         Density(ny(k),nx(k)) = Density(ny(k),nx(k)) + 1;
         XFlow(ny(k),nx(k)) = XFlow(ny(k),nx(k)) + X(k+1) - X(k);
         YFlow(ny(k),nx(k)) = YFlow(ny(k),nx(k)) + Y(k+1) - Y(k);
      end
      ndx = find(Density(:)>0);
      XFlow(ndx) = XFlow(ndx)./Density(ndx)/rsMUA.dt;
      YFlow(ndx) = YFlow(ndx)./Density(ndx)/rsMUA.dt;
      dx = diff(X_RANGE)/GRID_SIZE(1);
      dy = diff(Y_RANGE)/GRID_SIZE(2);
      Density = Density/(numel(rsMUA.value)-1)/(dx*dy);
      
      % Compute contour levels containing different percentiles of samples...
%       Prcs = 5:5:50;
      Prcs = 10:10:90;
      sD = sort(Density(:));
      csD = cumsum(sD)*(dx*dy);
      csD = csD + 1-csD(end);
      ContourLevels = zeros(size(Prcs));
      for k = 1:numel(Prcs)
         ContourLevels(k) = sD(find(csD>=Prcs(k)/100,1,'first'));
      end
      Z_RANGE = [0 ContourLevels(end)];
      
      figure
      hold on
      xx = linspace(X_RANGE(1),X_RANGE(2),GRID_SIZE(1)+1);
      xx = xx(1:end-1) + dx/2;
      yy = linspace(Y_RANGE(1),Y_RANGE(2),GRID_SIZE(2)+1);
      yy = yy(1:end-1) + dy/2;
      %    imagesc(xx,yy,Density)
      contourf(xx,yy,Density,linspace(Z_RANGE(1),Z_RANGE(2),numel(Prcs)*4))
      contour(xx,yy,Density,ContourLevels,'w')
      shading flat
      colormap(1-gray())
      brighten(-0.75)
      colorbar
      caxis(Z_RANGE)
      set(gca,'YDir','norm')
      
      grid('on')
      
      hold on
      %    NormX =  std(XFlow(Density(:)>MinDensity))*sqrt(dy/dx);
      %    NormY =  std(YFlow(Density(:)>MinDensity))*sqrt(dx/dy);
%       NormX =  std(XFlow(Density(:)>MinDensity))/1.25/dx;
%       NormY =  std(YFlow(Density(:)>MinDensity))/1.25/dy;
      NormX =  std(XFlow(Density(:)>MinDensity/(dx*dy)))/1/dx;
      NormY =  std(YFlow(Density(:)>MinDensity/(dx*dy)))/1/dy;
      for r = 1:GRID_SIZE(2)
         for c = 1:GRID_SIZE(1)
%             if Density(r,c) > MinDensity
            if Density(r,c) > MinDensity/(dx*dy)
               plot(xx(c),yy(r),'r.');
               plot(xx(c)+[0 XFlow(r,c)/NormX],yy(r)+[0 YFlow(r,c)/NormY],'r')
            end
         end
      end
      
      xlabel('log(MUA)')
      if FIXED_PLOT_RANGE
         ylabel('LFP (\muV)')
      else
         ylabel('LFP/st.dev.')
      end
      title('Density of samples')
      
      set(gca,'XLim',X_RANGE,'YLim',Y_RANGE);
      set(gca,'Layer','top','TickDir','out','Box','on');
      
      set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [0 0 5 4]);
      if FIXED_PLOT_RANGE
         print('-deps2c', 'LFPvsLogMUA_DensityAndFlow.eps');
      else
         print('-deps2c', 'LFPvsLogMUA_DensityAndFlow_ZoomAll.eps');
      end
   end % for FIXED_PLOT_RANGE = ...



   %% Saves some analysis results...
   %
%    save('AnalysisSummary.mat', 'ModeParams', 'UDcycle', 'UpDuration', ...
%       'UpState', 'DownState', 'LogMUAbaseline','PyyBaseline','LogMUAReference');
   save('AnalysisSummary.mat', 'ModeParams', ...
      'UpState', 'DownState', 'LogMUAbaseline','PyyBaseline','LogMUAReference');


   %% Returns to the root directory...
   %
   cd('..');
end


%% Analysis loop on recording set...
%
SDs = NaN(1,numel(RecordingSet));
for crsRecSet = 1:length(RecordingSet)
% for crsRecSet = length(RecordingSet)

   workingdir = [AnalysisDir RecordingSet(crsRecSet).label];
   load([workingdir '/AnalysisSummary.mat'], 'ModeParams');

   PeakMean = Inf;
   if isfield(ModeParams,'Mu1')
      PeakMean = ModeParams.Mu1;
      PeakStd = ModeParams.Sigma1;
   end
   if isfield(ModeParams,'Mu2')
      if PeakMean > ModeParams.Mu2
         PeakMean = ModeParams.Mu2;
         PeakStd = ModeParams.Sigma2;
      end
   end
   SDs(crsRecSet) = PeakStd;
end


%% To check the stability of a MUA threshold across channels.
%  In other words if the same threshold can be used for all the channels.
%  If yes, this should provide more reliable profiles of traveling
%  wavefronts.
%
figure

stem(SDs)
set(gca,'YLim',[0 0.35])
xlabel('Channels')
ylabel('Down MUA SD')
title(['Median SD = ' num2str(median(SDs),3)])
set(gcf, 'PaperUnits', 'inch', 'PaperPosition', [0 0 5 4]);
print('-deps2c', 'DownMUASDsVsChannels.eps');


catch
   msg = lasterror;
   fprintf('%s', msg.message)
   cd('..');
   if ~isempty(hSMR)
      closeSON(hSMR);
   end
   rethrow(lasterror);
end

if ~isempty(hSMR)
   closeSON(hSMR);
end
