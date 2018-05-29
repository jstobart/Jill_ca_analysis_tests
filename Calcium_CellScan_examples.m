
%% Example- basic CellScan

% Create the CellScan object interactively
cs001 = CellScan();  % load image, select channels, etc.

% Process the CellScan object
cs001.process();

% Produce a plot
cs001.plot();

cs001.opt_config();


%% FLIKA Parameter examples

BL_frames= 59;  %number of frames before stimulation

% For cytosolic GCaMP- 128x128 pixels, 12Hz

findConf{1} = ConfigFindROIsFLIKA_2D.from_preset('ca_cyto_astro', 'baselineFrames',...
    BL_frames,'sigmaXY', 2,...
    'sigmaT', 0.5,'thresholdPuff', 5,...
    'minRiseTime',0.1689, 'maxRiseTime', 8,...
    'dilateXY', 4,...
    'dilateT', 0.5,'erodeXY', 2);

% For membrane GCaMP- 128x128 pixels, 12Hz

findConf{2} = ConfigFindROIsFLIKA_2D.from_preset('ca_memb_astro', 'baselineFrames',...
    BL_frames,'freqPassBand',1,'sigmaXY', 2,...
    'sigmaT', 0.1,'thresholdPuff', 7, 'threshold2D', 0.2,...
    'minRiseTime',0.0845, 'maxRiseTime', 1,'minPuffArea', 10,...
    'dilateXY', 5, 'dilateT', 0.3,'erodeXY', 1, 'erodeT', 0.1);


%% measure ROIs and signal detection
                
% for calculating AUC for each trace
measureConf = ConfigMeasureROIsDummy('baselineFrames', BL_frames);

% most peaks are "single peaks"
detectConf{1} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...
    'propagateNaNs', false, 'excludeNaNs', false, 'lpWindowTime', 3, 'spFilterOrder', 2,...
    'spPassBandMin',0.05, 'spPassBandMax', 0.2, 'thresholdLP', 3,'thresholdSP', 5);


% mix of single, multi and plateaus
detectConf{2} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...'normMethod', 'z-score','zIters', 100,...
    'propagateNaNs', false, 'excludeNaNs', false, 'lpWindowTime', 5, 'spFilterOrder', 6,...
    'spPassBandMin',0.025, 'spPassBandMax', 0.0667, 'thresholdLP', 7,'thresholdSP', 7);
                

% put all the configuration parameters together
configCS{1,1} = ConfigCellScan(findConf{1,1}, measureConf,detectConf{1,1}); 
configCS{1,2} = ConfigCellScan(findConf{1,2}, measureConf,detectConf{1,2}); 
% the variable above can also be saved and loaded each time (instead of
% writing out the above settings)

%% CellScan with predefined parameters

FLIKA_Example1 = CellScan([], [], configCS{1,1}, 1);  % last parameter is the channel to use (change to 2 if astrocytes are the second channel)

FLIKA_Example2 = CellScan([], [], configCS{1,2}, 1);

%channel = struct('Ca_Cyto_Astro',1,'blood_plasma',2);
