
%% 3 lck long stim
clearvars

AllData= [];
All_traces= [];


%% Information about your images
Settings.MainDir = 'E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f';

Settings.AnimalNames = {
    'RG14',...
    'RG16',...
    'RG17',...
    'RG18',...
    };
Settings.ScoreSheetNames = {
    'RG14_Scoresheet_LongStim.xls',...
    'RG16_Scoresheet_LongStim.xls',...
    'RG17_Scoresheet_LongStim.xls',...
    'RG18_Scoresheet_LongStim.xls',...
    };
Settings.NameConditions = {'Stim'};

channel = struct('Ca_Memb_Astro',1,'Ca_Neuron',2);
plotMotion =0;

% final data file name
SaveFiles{1,1} = fullfile(Settings.MainDir, 'Results', 'LckGC&RC_2D_longstim_28_04_2017.csv');%'Control_Peaks_3Conds.csv'); % all data
SaveFiles{1,2} = fullfile(Settings.MainDir, 'Results', 'LckGC&RC_2D_longstim_28_04_2017.mat');%'Control_Peaks_3Conds.csv'); % all data
SaveFiles{1,3}= fullfile(Settings.MainDir, 'Results','LckGC&RC_traces_2D_longstim_28_04_2017.mat'); %'Control_TraceAUC_20sWindow_3Conds.csv'); % neuronal hand click cell scan

%% Load calibration file
calibration ='E:\matlab\2p-img-analysis\tests\res\calibration_20x.mat';
CalFile = Calibration_PixelSize.load(calibration);

%% make a mixing matrix for spectral unmixing of RCaMP and GCaMP

if ~exist(fullfile(Settings.MainDir, 'Results','RCaMP_mGCaMP_Matrix.mat'),'file')
    % load example high res image
    unmixFile = 'E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\RG14\Awake\test\highres_spot1_long021.tif';
    unmixImg = SCIM_Tif(unmixFile, channel, CalFile);
    
    [unmixImg, RCaMP_mGCaMP_Matrix] = unmix_chs(unmixImg); %returns the mixing matrix to be used for all imaging
    
    % save the mixing matrix for loading later
    cd(fullfile(Settings.MainDir, 'Results'));
    % write matrix to created file
    save('RCaMP_mGCaMP_Matrix.mat', 'RCaMP_mGCaMP_Matrix');
else
    load(fullfile(Settings.MainDir, 'Results','RCaMP_mGCaMP_Matrix.mat'));
end


%% load scoresheet and loop through animal, spot, etc.

Settings.ScoreSheetPath = fullfile(Settings.MainDir,'Scoresheets',Settings.ScoreSheetNames);
CurrentSheet = Settings.ScoreSheetPath{1,1};
Settings = readScoresheet2(CurrentSheet, Settings);

BL_frames=59;
% Create an array of ScanImage Tiffs
ImgArray =  SCIM_Tif([], channel, CalFile);

% examples:
% RG14, 03_08_2016, spot1, trial 10, ROI3, Neuron 8,9
% RG14, 02_26_2016, spot1, trial 3, process 15 vs dendrites
% RG17, 03_08_2016, spot1, trial 3, process 15 vs dendrites

% Spectral Unmixing of GCaMP and RCaMP
ImgArray= ImgArray.unmix_chs(false, [], cell2mat(RCaMP_mGCaMP_Matrix));

HighRes = SCIM_Tif([],channel, CalFile);
% Extract a reference image
refImg = mean(HighRes.rawdata(:,:,2,:),4);
%figure, imagesc(refImg), axis image, axis off, colormap(gray)

ImgArray = ImgArray.motion_correct('refImg',refImg, 'ch',2,'maxShift', 10,'minCorr', 0.4);%,'doPlot',true);


%% Configs for Finding ROIs
%ASTROCYTES
%AWAKE astrocyte membrane calcium
% 2D automated selection for peaks
AC_findConf{1} = ConfigFindROIsFLIKA_2D.from_preset('ca_memb_astro', 'baselineFrames',...
    BL_frames,'freqPassBand',1,'sigmaXY', 2,...
    'sigmaT', 0.1,'threshold_std', 7, 'threshold_2D', 0.2,...
    'min_rise_time',0.0845, 'max_rise_time', 1,'minPuffArea', 10,...
    'dilateXY', 5, 'dilateT', 0.3,'erodeXY', 1, 'erodeT', 0.1,...
    'discardBorderROIs',true);

AC_findConf{2} = ConfigFindROIsFLIKA_2p5D.from_preset('ca_memb_astro', 'baselineFrames',...
    BL_frames,'freqPassBand',1,'sigmaXY', 2,...
    'sigmaT', 0.1,'threshold_std', 7, 'threshold_2D', 0.2,...
    'min_rise_time',0.0845, 'max_rise_time', 1,'minPuffArea', 10,...
    'dilateXY', 5, 'dilateT', 0.3,'erodeXY', 1, 'erodeT', 0.1,...
    'discardBorderROIs',true);

% 3D automated selection for time and space estimations
AC_findConf{3} = ConfigFindROIsFLIKA_3D.from_preset('ca_memb_astro', 'baselineFrames',...
    BL_frames,'freqPassBand',1,'sigmaXY', 2,...
    'sigmaT', 0.1,'threshold_std', 7, 'threshold_2D', 0.2,...
    'min_rise_time',0.0845, 'max_rise_time', 1,'minPuffArea', 10,...
    'dilateXY', 5, 'dilateT', 0.3,'erodeXY', 1, 'erodeT', 0.1,...
    'discardBorderROIs',true);

% NEURONS
% 2D FLIKA selected for peaks from "dendrites"
Neur_findConf{1} = ConfigFindROIsFLIKA_2D.from_preset('ca_neuron', 'baselineFrames',...
    BL_frames,'freqPassBand',1,'sigmaXY', 2,...
    'sigmaT', 0.1,'threshold_std', 7, 'threshold_2D', 0.2,...
    'min_rise_time',0.0845, 'max_rise_time', 2,'minPuffArea', 10,...
    'dilateXY', 3, 'dilateT', 0.5,'erodeXY', 2, 'erodeT', 0.3,...
    'discardBorderROIs',true);

% hand selected for peaks from somata
x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
Neur_findConf{2} = ConfigFindROIsDummy.from_ImageJ([], x_pix, y_pix,1);

% 3D FLIKA selected for time and space estimations
Neur_findConf{3} = ConfigFindROIsFLIKA_3D.from_preset('ca_neuron', 'baselineFrames',...
    BL_frames,'freqPassBand',1,'sigmaXY', 2,...
    'sigmaT', 0.1,'threshold_std', 7, 'threshold_2D', 0.2,...
    'min_rise_time',0.0845, 'max_rise_time', 2,'minPuffArea', 10,...
    'dilateXY', 3, 'dilateT', 0.5,'erodeXY', 2, 'erodeT', 0.3,...
    'discardBorderROIs',true);

%% Configuration for measuring ROIs


detectConf{1} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...'normMethod', 'z-score','zIters', 100,...
    'propagateNaNs', false, 'excludeNaNs', false, 'lpWindowTime', 1.5, 'spFilterOrder', 2,...
    'spPassBandMin',0.05, 'spPassBandMax', 0.5, 'thresholdSD_low', 3,'thresholdSD_band', 5);

% AWAKE neuron calcium
detectConf{2} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,... 'normMethod','z-score',... 'zIters', 10000,...
    'propagateNaNs', false,'excludeNaNs', false, 'lpWindowTime', 2, 'spFilterOrder', 2,...
    'spPassBandMin',0.1, 'spPassBandMax', 1, 'thresholdSD_low', 3,'thresholdSD_band', 5);

% for 3D FLIKA
detectConf{3} = ConfigDetectSigsDummy();

% for calculating AUC for each trace
measureConf = ConfigMeasureROIsDummy('baselineFrames', BL_frames);

%%
% Combine the configs into a CellScan config
configCS{1,1} = ConfigCellScan(AC_findConf{1,1}, measureConf,detectConf{1,1}); % astrocyte 2D FLIKA
configCS{1,2} = ConfigCellScan(AC_findConf{1,2}, measureConf,detectConf{1,3}); % astrocyte 2.5D FLIKA
configCS{1,3} = ConfigCellScan(AC_findConf{1,3}, measureConf,detectConf{1,3}); % astrocyte 3D FLIKA
configCS{1,4} = ConfigCellScan(Neur_findConf{1,1}, measureConf,detectConf{1,2}); % neuronal 2D FLIKA
configCS{1,5} = ConfigCellScan(Neur_findConf{1,2}, measureConf,detectConf{1,2}); % neuronal hand peaks
configCS{1,6} = ConfigCellScan(Neur_findConf{1,3}, measureConf,detectConf{1,3}); % neuronal 3D FLIKA

%% Create CellScan objects
CSArray_Ch1_FLIKA_2D = CellScan([], ImgArray, configCS{1,1}, 1);
CSArray_Ch1_FLIKA_2p5D = CellScan([], ImgArray, configCS{1,2}, 1);
CSArray_Ch1_FLIKA_3D = CellScan([], ImgArray, configCS{1,3}, 1);

CSArray_Ch2_Hand = CellScan([], ImgArray, configCS{1,5}, 2);
CSArray_Ch2_FLIKA_2D= CellScan([], ImgArray, configCS{1,4}, 2);
CSArray_Ch2_FLIKA_3D= CellScan([], ImgArray, configCS{1,6}, 2);




%% Process the images
CSArray_Ch1_FLIKA_3D =CSArray_Ch1_FLIKA_3D.process();
CSArray_Ch1_FLIKA_2D =CSArray_Ch1_FLIKA_2D.process();
CSArray_Ch1_FLIKA_2p5D =CSArray_Ch1_FLIKA_2p5D.process();
CSArray_Ch2_Hand =CSArray_Ch2_Hand.process();
CSArray_Ch2_FLIKA_3D =CSArray_Ch2_FLIKA_3D.process();
CSArray_Ch2_FLIKA_2D =CSArray_Ch2_FLIKA_2D.process();


%% View Movies
CSArray_Ch1_FLIKA_2D.plot();
CSArray_Ch1_FLIKA_2D.plot('video');
%CSArray_Ch1_FLIKA_2D.plot('video','plotROIs',3);
CSArray_Ch1_FLIKA_2D.plot('plotROIs',15);
CSArray_Ch1_FLIKA_2D.plot('video','plotROIs',15);

CSArray_Ch1_FLIKA_2p5D.plot();
CSArray_Ch1_FLIKA_2p5D.plot('video');

CSArray_Ch1_FLIKA_3D.plot('video');


CSArray_Ch2_FLIKA_2D.plot();
CSArray_Ch2_FLIKA_2D.plot('video');
CSArray_Ch2_FLIKA_2D.plot('video','plotROIs',5);

CSArray_Ch2_FLIKA_3D.plot('video');

CSArray_Ch2_Hand.plot('video','plotROIs',[8,9]);
CSArray_Ch2_Hand.plot();


%%
filename='RG17_2016_08_03_trial8';
write_tiff_stacks(CSArray_Ch1_FLIKA_2D, filename);

% for working out find peaks parameters
%CSArray_Ch2_FLIKA(1).opt_config()



