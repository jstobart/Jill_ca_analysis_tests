close all; clc;

%% Example script for CellScan application
% Add path to 2p-img-analysis
addpath('D:\Code\Matlab\2p-img-analysis');
% Add path to experiment directory
addpath(genpath('D:\Experiments\Tendon Biomechanics'));

cd(utils.CHIPS_rootdir)

%% First, establish variables
% Image file and channels
fNameImg = 'D:\Experiments\Tendon Biomechanics\resources\Sample 1 RTF deconvolved.tif';
info = imfinfo(fNameImg);
num_im = length(info);

channels = struct('Ca_Cyto_Astro', 1, 'blood plasma', 2);

% Calibration 
fNameCal = fullfile(utils.CHIPS_rootdir, 'tests/res', 'calibration_20x.mat');
calibration = Calibration_PixelSize.load(fNameCal);

%% Instatiate raw image
% for TF8 images
rawImg = SCIM_Tif(fNameImg, channels, calibration, false);

%% Configuration for finding ROIs (hand-selected)
zipPath = 'D:/....';
xPix = 128;
yPix = 127;
scaleF = 1;
findConf = ConfigFindROIsDummy.from_imageJ(zipPath, xPix, yPix, scaleF);

findConf = ConfigFindROIsFLIKA_3D.from_preset('ca_memb_astro', 'freqPassBand',1,'sigmaXY', 2,...
    'sigmaT', 0.1,'threshold_std', 7, 'threshold2D', 0.2,...
    'min_rise_time',0.0845, 'max_rise_time', 1,'minPuffArea', 10,...
    'dilateXY', 5, 'dilateT', 0.3,'erodeXY', 1, 'erodeT', 0.1,...
    'discardBorderROIs',true);

%% Configuration for measuring ROIs
% Change these parameters in case you're not happy with the peaks it
% identified or missed
bgLevel = 15; % background in percentile

measureConf = ConfigMeasureROIsZScore(...
    'backgroundLevel', bgLevel, ...
    'zIters', 10, ...
    'zSDs', 2, ...
    'propagateNaNs', true);

detectConf = ConfigDetectSigsClsfy(); % kann verändert werden


% Combine the three configs into a CellScan config
configCS = ConfigCellScan(findConf, measureConf, detectConf);

%% Run preprocessing steps
% Spectral unmixing
rawImg.unmix_ch

% Motion correction
channelToUseMC = 2; % which channel to use
refImg = squeeze(mean(rawImg.rawdata(:,:,channelToUseMC, 1:20),4));
rawImg.motion_correct('ch', channelToUseMC, 'fillBadData', 'inpaint', ...
    'refImg', refImg);

%% Instantiate CellScan 
testCS = CellScan(fNameImg, rawImg, configCS);

%% Process CellScan
testCS.process;

%% Plot results
testCS.plot;
testCS.opt_config;
%% Output results
testCS.output_data;