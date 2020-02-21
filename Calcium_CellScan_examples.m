% Fangying

% load in the data from the confocal
Images = RawImgDummy();  % make an artifical image stack


% fnRID003 = fullfile(utils.CHIPS_rootdir, 'tests', 'res', ...
%     'ArtificialData_SNR_0.1_FOV_250.tif');
%%
% Specify the channels relevant for this raw image
channels003 = struct('cellular_signal', 1);
calibration='C:\JillsFiles\matlab\Jill_ca_analysis_tests\CalibrationFiles\calibration_dummy.mat';
% Specify some data about the image acquisition
nRows003 = 1024;
nCols003 = 1024;
acq003 = struct('isBiDi', false, 'lineTime', 300000/nRows003, 'zoom', 1, ...
    'nLinesPerFrameOrig', nRows003, 'nPixelsPerLineOrig', nCols003, 'pixelSize', 0.4151);

% Create the RawImgDummy object without any interaction
rid003 = RawImgDummy([], channels003, calibration, acq003);
%%
% View the RawImgDummy object
rid003.plot()

Alice1.plot();

% unmix the two fluorophores
% an unmixing matrix must be created for each combination of sensors
% (GCaMP/RCaMP, RCaMP/FITC Dextran, GCaMP/Texas Red Dextran)

Alice1b= Alice1.unmix_chs();
Alice1b.plot();

% motion correction in XY

% open or create a reference image

refImg = BioFormats(); % load high resolution image (but it must be the same zoom factor!)
refImg = mean(Alice1b.rawdata(:,:,1,:),4);  % make a reference image by taking the average of the first 10 frames from the T-series


Alice1c= Alice1.motion_correct('refImg',refImg, 'ch',1,'doPlot',true);
           
Alice1.plot();


%% Example- basic CellScan

% Create a CellScan object with the images that have been imported and
% processed
cs001 = CellScan([], Alice1, [], []);

% Create the CellScan object interactively
cs001 = CellScan();  % load image, select channels, etc.

% Process the CellScan object
cs001.process();

% Produce a plot
cs001.plot();

% change the parameters
cs001.opt_config();

% Output the data
cs001.output_data();

% TODO: show the data file(s)

% TODO: show where to retrieve the data from inside the object


%% Use ImageJ Rois

% TODO: Show how to create the ImageJ ROIs

% Create the CellScan object interactively
cs002 = CellScan([], Alice1);

% Process and plot the CellScan object
cs002.process();
cs002.plot();

% Look at some alternative plots
cs002.plot('rois');
cs002.plot('traces', 'doWholeFrame', false);

%% Use FLIKA

% Create the CellScan object interactively
cs003 = CellScan([], Alice1);

% Process and plot the CellScan object
cs003.process();
cs003.plot();

% Look at some alternative plots
cs003.plot('rois');
cs003.plot('traces', 'doWholeFrame', false);

%% Use FLIKA with alternate parameters


BL_frames= 20;  %number of frames before stimulation


% For membrane GCaMP- 128x128 pixels, 12Hz

findConf{1} = ConfigFindROIsFLIKA_2D.from_preset('ca_memb_astro', 'baselineFrames',...
    BL_frames,'freqPassBand',1,'sigmaXY', 2,...
    'sigmaT', 0.1,'thresholdPuff', 3, 'threshold2D', 0.2,...
    'minRiseTime',0.07, 'maxRiseTime', 1,'minROIArea', 10,...
    'dilateXY', 5, 'dilateT', 0.3,'erodeXY', 1, 'erodeT', 0.1);


%% measure ROIs and signal detection
                
% for calculating AUC for each trace
measureConf = ConfigMeasureROIsDummy('baselineFrames', BL_frames);

% most peaks are "single peaks"
detectConf{1} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...'normMethod', 'z-score','zIters', 100,...
                'propagateNaNs', false, 'excludeNaNs', false, 'lpWindowTime', 1.5, 'spFilterOrder', 2,...
                'spPassBandMin',0.05, 'spPassBandMax', 0.5, 'thresholdLP', 3,'thresholdSP', 5);
             

% put all the configuration parameters together
configCS{1,1} = ConfigCellScan(findConf{1}, measureConf,detectConf{1,1}); 

% the variable above can also be saved and loaded each time (instead of
% writing out the above settings)

%% CellScan with predefined parameters


% last parameter is the channel to use (change to 2 if astrocytes are the second channel)


%channel = struct('Ca_Cyto_Astro',1,'blood_plasma',2);


%% Create a simple CellScan object without any interaction

% TODO: Point out that basically everything can be scripted


