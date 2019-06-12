%% ASTROCYTE-NEURON CALCIUM,  RCaMP T-SERIES for CALCIUM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Overview of CHIPS processing

% You will need CHIPS folder, Bfmatlab folder, load_prairie_line script, xml2struct script on your path

% - Create RawImg
% - Do things to it (mix up, combine, calculate, motion correct)

%   ch_calc         - Perform mathematical calculations on image channels
%   ch_calc_ratio   - Calculate the ratio of two image channels
%   ch_calc_sum2    - Calculate the sum of two image channels
%   check_ch        - Check that the appropriate channels are present
%   copy            - Copy MATLAB array of handle objects
%   denoise         - Denoise the images
%   downsample      - Downsample the images in space and/or time
%   exclude_frames  - Excludes the specified frame(s) from the image
%   get_ch_name     - Get channel names from channel numbers
%   get_mc          - Get the motion correction information
%   has_ch          - Determine if particular channels are present
%   motion_correct  - Motion correct the images
%   plot            - Plot a figure
%   split1          - Split the image data along a given dimension
%   to_long         - Convert the images to long format
%   unmix_chs       - Unmix image channels
%
% RawImgDummy static methods:
%   cat_data        - Concatenate the data from RawImgDummy objects
%   from_files      - Create a RawImgDummy object from a list of files


% - Create Config
% - Create ProcessedImg
% - Process
% - Plot
% - Output and/or analyse results in more detail


% TODO: discuss parallel processing

%% Example- load some data and adjust it

Alice1 = BioFormats();  % BioFormats plugin

% Alice 2019.06.03 data

Alice1.plot();

% unmix the two fluorophores
% an unmixing matrix must be created for each combination of sensors
% (GCaMP/RCaMP, RCaMP/FITC Dextran, GCaMP/Texas Red Dextran)

Alice1= Alice1.unmix_chs();
Alice1.plot();

% motion correction in XY

% open or create a reference image

refImg = BioFormats(); % load high resolution image (but it must be the same zoom factor!)
refImg = mean(Alice1.rawdata(:,:,2,1:10),4);  % make a reference image by taking the average of the first 10 frames from the T-series


Alice1= Alice1.motion_correct('refImg',refImg, 'ch',2,'maxShift', 10,'minCorr', 0.4);%,'doPlot',true);
           
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

findConf{2} = ConfigFindROIsFLIKA_2D.from_preset('ca_memb_astro', 'baselineFrames',...
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
configCS{1,1} = ConfigCellScan(findConf{1,1}, measureConf,detectConf{1,1}); 

% the variable above can also be saved and loaded each time (instead of
% writing out the above settings)

%% CellScan with predefined parameters

cs004 = CellScan([], [], configCS{1,1}, 1);  % last parameter is the channel to use (change to 2 if astrocytes are the second channel)


%channel = struct('Ca_Cyto_Astro',1,'blood_plasma',2);


%% Create a simple CellScan object without any interaction

% TODO: Point out that basically everything can be scripted


