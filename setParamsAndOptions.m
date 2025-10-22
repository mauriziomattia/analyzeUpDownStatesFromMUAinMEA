% setParamsAndOptions.m
%
%   Load analysis params and information about recording dataset.
%
%   Copyright 2015 Maurizio Mattia @ Ist. Super. Sanita', Rome - Italy
%   Version: 1.0 - Feb. 18, 2015
%


%% Load analysis params and information about recording dataset.
Options.PeriodToAnalyze = [0 500];

Options.SaveMUA = 0;

Options.LogMUA.FreqBand = [200 1500];
Options.LogMUA.SmoothingWindow = 0.040; % It depends on the sampling rate.
Options.UD.ThresholdModulation = 0.6;   % Between 0 and 1. 0.5 is the mean of the Up and Down levels.

Options.LFP.SmoothingWindow = Options.LogMUA.SmoothingWindow; % s...
Options.LFP.SubsamplingRate = Options.LogMUA.FreqBand(1);   % Hz...

Options.LFPMUAplot.MUArange = [-0.75 3.25];
Options.LFPMUAplot.LFPrange = [-2000 2000];

% N.B. struct UpTimeLags depends on the electrode array
Options.UpTimeLags.MaxAbsTimeLag = 0.800;       % Maximum reasonable time lag between electrodes...
% Options.UpTimeLags.ReferenceChannel = 11;        % Channel whose triggers are used as reference. If not defined use channel close to the middle is automatically selected.
Options.UpTimeLags.SortingChannels = [05 06 07]; % A guess to sorting timelags array. Usually the channels belonging to the same cluster.
Options.UpTimeLags.PCNumOfSTD = 3;              % Num of st.dev. in the PC space to avoid outliers in the PCA.
% % Options.UpTimeLags.InvertElectrodePos = 1;      % If exist and is 1 invert the order of the electrodes position.
Options.UpTimeLags.NumOfClusters = 10;          % Number of clusters of wavefronts.
Options.UpTimeLags.WavefrontTimeGap = 20;       % Time gap between wavefront in ms.

BaseDir = 'C:/Users/2016.AnesthesiaFadingOutInECoGrecs/';
AnalysisDir = [BaseDir '161203/Light/'];
DataDir = 'C:/Users/antonio.pazienti/DATA/ISS/2016.AnesthesiaFadingOutInECoGrecs/';
DataFile = [DataDir '161203/161203_rec13_Propagation_000_Light.smr']; 

hemisphere = 'LH'; % or 'RH'
SelectedArea = 'A'; % A = all; or 'V', or 'R', or 'P' ... (see above)

% -------------------------------------------------------------------------
% ELECTRODE ARRAY - 32 channels
% -------------------------------------------------------------------------
% --> distinguish between RIGHT (RH) and LEFT (LH) Hemisphere
% --> define a reference system
% --> enable the opportunity to select a specific Cortical Area:
%   V = Visual Cortex
%   R = Retrosplenial Cortex 
%   P = Parietal Association Area (PtA)
%   S = Somatosensory Cortex
%   M = Motor Cortex



if strcmpi(hemisphere,'RH')
    
    n = 1;
    RecordingSet(n).label ='01.Ch27';
    RecordingSet(n).port = 'E27';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -0*0.550;
    RecordingSet(n).YPos = 0*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;

    n = n + 1;
    RecordingSet(n).label ='02.Ch18';
    RecordingSet(n).port = 'E18';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -1*0.550;
    RecordingSet(n).YPos = 0*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;

    n = n + 1;
    RecordingSet(n).label ='03.Ch28';
    RecordingSet(n).port = 'E28';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -0*0.550;
    RecordingSet(n).YPos = 1*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;

    n = n + 1;
    RecordingSet(n).label ='04.Ch19';
    RecordingSet(n).port = 'E19';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;    
    RecordingSet(n).XPos = -1*0.550;
    RecordingSet(n).YPos = 1*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;

    n = n + 1; 
    RecordingSet(n).label ='05.Ch10';
    RecordingSet(n).port = 'E10';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -2*0.550;
    RecordingSet(n).YPos = 1*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;

    n = n + 1; 
    RecordingSet(n).label ='06.Ch29';
    RecordingSet(n).port = 'E29';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -0*0.550;
    RecordingSet(n).YPos = 2*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;

    n = n + 1; 
    RecordingSet(n).label ='07.Ch20';
    RecordingSet(n).port = 'E20';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -1*0.550;
    RecordingSet(n).YPos = 2*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;

    n = n + 1; 
    RecordingSet(n).label ='08.Ch11';
    RecordingSet(n).port = 'E11';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -2*0.550;
    RecordingSet(n).YPos = 2*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;

    n = n + 1;
    RecordingSet(n).label ='09.Ch30';
    RecordingSet(n).port = 'E30';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -0*0.550;
    RecordingSet(n).YPos = 3*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
 
    n = n + 1;
    RecordingSet(n).label ='10.Ch21';
    RecordingSet(n).port = 'E21';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -1*0.550;
    RecordingSet(n).YPos = 3*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;

    n = n + 1;
    RecordingSet(n).label ='11.Ch12';
    RecordingSet(n).port = 'E12';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -2*0.550;
    RecordingSet(n).YPos = 3*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='12.Ch4';
    RecordingSet(n).port = 'E04';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -3*0.550;
    RecordingSet(n).YPos = 3*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='13.Ch31';
    RecordingSet(n).port = 'E31';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -0*0.550;
    RecordingSet(n).YPos = 4*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='14.Ch22';
    RecordingSet(n).port = 'E22';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration = 0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -1*0.550;
    RecordingSet(n).YPos = 4*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='15.Ch13';
    RecordingSet(n).port = 'E13';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -2*0.550;
    RecordingSet(n).YPos = 4*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='16.Ch5';
    RecordingSet(n).port = 'E05';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -3*0.550;
    RecordingSet(n).YPos = 4*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='17.Ch32';
    RecordingSet(n).port = 'E32';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -0*0.550;
    RecordingSet(n).YPos = 5*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='18.Ch23';
    RecordingSet(n).port = 'E23';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration = 0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -1*0.550;
    RecordingSet(n).YPos = 5*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='19.Ch14';
    RecordingSet(n).port = 'E14';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -2*0.550;
    RecordingSet(n).YPos = 5*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='20.Ch6';
    RecordingSet(n).port = 'E06';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -3*0.550;
    RecordingSet(n).YPos = 5*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='21.Ch33';
    RecordingSet(n).port = 'E33';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -0*0.550;
    RecordingSet(n).YPos = 6*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='22.Ch24';
    RecordingSet(n).port = 'E24';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -1*0.550;
    RecordingSet(n).YPos = 6*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='23.Ch15';
    RecordingSet(n).port = 'E15';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration = 0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -2*0.550;
    RecordingSet(n).YPos = 6*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='24.Ch7';
    RecordingSet(n).port = 'E07';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -3*0.550;
    RecordingSet(n).YPos = 6*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='25.Ch3';
    RecordingSet(n).port = 'E03';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -4*0.550;
    RecordingSet(n).YPos = 6*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='26.Ch34';
    RecordingSet(n).port = 'E34';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -0*0.550;
    RecordingSet(n).YPos = 7*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='27.Ch25';
    RecordingSet(n).port = 'E25';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -1*0.550;
    RecordingSet(n).YPos = 7*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='28.Ch16';
    RecordingSet(n).port = 'E16';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -2*0.550;
    RecordingSet(n).YPos = 7*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='29.Ch8';
    RecordingSet(n).port = 'E08';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration = 0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -3*0.550;
    RecordingSet(n).YPos = 7*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='30.Ch26';
    RecordingSet(n).port = 'E26';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -1*0.550;
    RecordingSet(n).YPos = 8*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='31.Ch17';
    RecordingSet(n).port = 'E17';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -2*0.550;
    RecordingSet(n).YPos = 8*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='32.Ch9';
    RecordingSet(n).port = 'E09';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = -3*0.550;
    RecordingSet(n).YPos = 8*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
        
elseif strcmpi(hemisphere,'LH')
    
    n = 1;
    RecordingSet(n).label ='01.Ch10';
    RecordingSet(n).port = 'E10';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 0*0.550;
    RecordingSet(n).YPos = 0*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;

    n = n + 1;
    RecordingSet(n).label ='02.Ch19';
    RecordingSet(n).port = 'E19';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 0*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;

    n = n + 1;
    RecordingSet(n).label ='03.Ch9';
    RecordingSet(n).port = 'E09';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 0*0.550;
    RecordingSet(n).YPos = 1*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;

    n = n + 1;
    RecordingSet(n).label ='04.Ch18';
    RecordingSet(n).port = 'E18';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;    
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 1*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;

    n = n + 1; 
    RecordingSet(n).label ='05.Ch27';
    RecordingSet(n).port = 'E27';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 2*0.550;
    RecordingSet(n).YPos = 1*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;

    n = n + 1; 
    RecordingSet(n).label ='06.Ch8';
    RecordingSet(n).port = 'E08';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 0*0.550;
    RecordingSet(n).YPos = 2*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;

    n = n + 1; 
    RecordingSet(n).label ='07.Ch17';
    RecordingSet(n).port = 'E17';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 2*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;

    n = n + 1; 
    RecordingSet(n).label ='08.Ch26';
    RecordingSet(n).port = 'E26';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 2*0.550;
    RecordingSet(n).YPos = 2*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;

    n = n + 1;
    RecordingSet(n).label ='09.Ch7';
    RecordingSet(n).port = 'E07';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 0*0.550;
    RecordingSet(n).YPos = 3*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
 
    n = n + 1;
    RecordingSet(n).label ='10.Ch16';
    RecordingSet(n).port = 'E16';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 3*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;

    n = n + 1;
    RecordingSet(n).label ='11.Ch25';
    RecordingSet(n).port = 'E25';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 2*0.550;
    RecordingSet(n).YPos = 3*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='12.Ch33';
    RecordingSet(n).port = 'E33';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 3*0.550;
    RecordingSet(n).YPos = 3*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='13.Ch6';
    RecordingSet(n).port = 'E06';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 0*0.550;
    RecordingSet(n).YPos = 4*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='14.Ch15';
    RecordingSet(n).port = 'E15';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration = 0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 4*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='15.Ch24';
    RecordingSet(n).port = 'E24';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 2*0.550;
    RecordingSet(n).YPos = 4*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='16.Ch32';
    RecordingSet(n).port = 'E32';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 3*0.550;
    RecordingSet(n).YPos = 4*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='17.Ch5';
    RecordingSet(n).port = 'E05';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 0*0.550;
    RecordingSet(n).YPos = 5*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='18.Ch14';
    RecordingSet(n).port = 'E14';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration = 0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 5*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='19.Ch23';
    RecordingSet(n).port = 'E23';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 2*0.550;
    RecordingSet(n).YPos = 5*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='20.Ch31';
    RecordingSet(n).port = 'E31';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 3*0.550;
    RecordingSet(n).YPos = 5*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='21.Ch4';
    RecordingSet(n).port = 'E04';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 0*0.550;
    RecordingSet(n).YPos = 6*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='22.Ch13';
    RecordingSet(n).port = 'E13';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 6*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='23.Ch22';
    RecordingSet(n).port = 'E22';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration = 0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 2*0.550;
    RecordingSet(n).YPos = 6*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='24.Ch30';
    RecordingSet(n).port = 'E30';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 3*0.550;
    RecordingSet(n).YPos = 6*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='25.Ch34';
    RecordingSet(n).port = 'E34';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 4*0.550;
    RecordingSet(n).YPos = 6*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='26.Ch3';
    RecordingSet(n).port = 'E03';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 0*0.550;
    RecordingSet(n).YPos = 7*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='27.Ch12';
    RecordingSet(n).port = 'E12';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 7*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='28.Ch21';
    RecordingSet(n).port = 'E21';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 2*0.550;
    RecordingSet(n).YPos = 7*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='29.Ch29';
    RecordingSet(n).port = 'E29';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration = 0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 3*0.550;
    RecordingSet(n).YPos = 7*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='30.Ch11';
    RecordingSet(n).port = 'E11';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 1*0.550;
    RecordingSet(n).YPos = 8*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='31.Ch20';
    RecordingSet(n).port = 'E20';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 2*0.550;
    RecordingSet(n).YPos = 8*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
    
    n = n + 1;
    RecordingSet(n).label ='32.Ch28';
    RecordingSet(n).port = 'E28';
    RecordingSet(n).AbsoluteThreshold = 0.40;
    RecordingSet(n).MinStateDuration =  0.080;
    RecordingSet(n).HPFCutOffFreq = 0.1;
    RecordingSet(n).XPos = 3*0.550;
    RecordingSet(n).YPos = 8*0.550;
    % RecordingSet(n).PCAforUpDownDetection = 1;
   
else 
    disp('Please, specify if Left Hemisphere or Right Hemisphere')
end
    

% --- SelectedArea ---
if SelectedArea == 'V'     % V = Visual Cortex
    ChannelSet = [18 19 20 22 23 24 25 27 28 29 30 31 32];
elseif SelectedArea == 'R' % R = Retrosplenial Cortex 
    ChannelSet = [13 17 21 26];
elseif SelectedArea == 'P' % P = Parietal Association Area (PtA)
    ChannelSet = [14 15 16];
elseif SelectedArea == 'S' % S = Somatosensory Cortex
    ChannelSet = [4 5 7 8 10 11 12];
elseif SelectedArea == 'M' % M = Motor Cortex
    ChannelSet = [1 2 3 6 9];
elseif SelectedArea == 'A' % A = All the electrodes
    ChannelSet = [1:32];
else
    disp('Please, specify the Selected Area of the Cortex')
end