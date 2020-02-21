% Fangying


% Specify the channels relevant for this raw image
channels003 = struct('Ca_Cyto_Astro', 1);

% load the made up calibration object (created from known pixel size from
% microscope
calibration='G:\Data\Fangying\FangYingCalibration.mat';
cal=load(calibration);

% Specify some data about the image acquisition
nRows003 = 1024;
nCols003 = 1024;
acq003 = struct('isBiDi', false, 'lineTime', 5000/nRows003, 'zoom', 1, ...
    'nLinesPerFrameOrig', nRows003, 'nPixelsPerLineOrig', nCols003 );

% load the tif stack saved in ImageJ from confocal czi file
rid003 = RawImgDummy([], channels003, cal.ans, acq003);

% View the RawImgDummy object
rid003.plot()



%% Define calcium analysis parameters

BL_frames= 7;

AC_findConf1 = ConfigFindROIsFLIKA_3D.from_preset('ca_memb_astro', 'baselineFrames',...
    BL_frames,'freqPassBand',0.001,'sigmaXY', 2,...
    'sigmaT', 5,'thresholdPuff', 5,'minROITime', 5,...
    'minRiseTime',5, 'maxRiseTime', 50,'minROIArea', 500,...
    'maxROIArea', 10000, 'dilateXY', 0, 'dilateT', 0,'erodeXY', 10, 'erodeT', 5);

AC_findConf2 = ConfigFindROIsFLIKA_2D.from_preset('ca_memb_astro', 'baselineFrames',...
    BL_frames,'freqPassBand',0.001,'sigmaXY', 2,...
    'sigmaT', 5,'thresholdPuff', 5,'threshold2D', 0,...
    'minRiseTime',5, 'maxRiseTime', 50,'minROIArea', 500,...
    'maxROIArea', 10000,'dilateXY', 0, 'dilateT', 0,'erodeXY', 10, 'erodeT', 5);

% measure ROIs (extract the traces)
AC_measureConf = ConfigMeasureROIsDummy('baselineFrames', BL_frames);

% filter the traces to detect the peaks and get info about them
AC_detectConf = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...
    'propagateNaNs', false, 'excludeNaNs', false, 'lpWindowTime', 1, 'spFilterOrder', 10,...
    'spPassBandMin',0.005, 'spPassBandMax', 0.05, 'thresholdLP', 5,'thresholdSP', 7);

% Combine the configs into a CellScan config for membrane tagged GCaMP
AC_configCS1= ConfigCellScan(AC_findConf1, AC_measureConf, AC_detectConf); % 3D region of interest (ROI) finding (unique regions in time)
AC_configCS2= ConfigCellScan(AC_findConf2, AC_measureConf, AC_detectConf); % 2D ROI finding (same region over time)

%% start the analysis

% create a cell scan object with the movie we entered and the parameters
% above

% 3D ROIs
cs001=CellScan([],rid003, AC_configCS1, 1);

cs001.process();

cs001.plot();
cs001.plot('video');
cs001.opt_config();

%% 2D ROIs

% 3D ROIs
cs002=CellScan([],rid003, AC_configCS2, 1);

cs002.process();

cs002.plot();

cs002.opt_config();
