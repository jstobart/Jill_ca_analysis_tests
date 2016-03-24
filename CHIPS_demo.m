%% Introduction

% TODO: Thank Jill and Kim for help.
%
%   Jill for initial setting up of of Ca2+ imaging + analysis, which
%   started the whole CellScan thing, and for lots of initial debugging and
%   development work.
%
%   Kim for implementing most of the CellScan functionality, especially
%   with the FLIKA algorithm, as well as lots of other utility functions,
%   debugging, ideas etc.

%% What CHIPS is and is not

% IS NOT: A magical program to do everything you could possibly ever
% imagine, taking your images and producing your final paper with
% stats, publication quality images etc

% IS: A tool to help you extract initial data from your images

%% Status, path forward

% 

%% Overview of CHIPS processing

% - Create RawImg
% - Do things to it (mix up, combine, calculate, motion correct)
% - Create Config
% - Create ProcessedImg
% - Process
% - Plot
% - Output and/or analyse results in more detail

%% Create a simple CellScan object interactively

% TODO: Point out that many simple tasks can be done interactively

% Create the CellScan object interactively
cs001 = CellScan();

% Process the CellScan object
cs001.process();

% Produce a plot
cs001.plot();

% Produce a slightly different plot (Look at the options)
cs001.plot('doWholeFrame', false, 'FilledROIs', false, 'AlphaSpec', 1);

% Output the data
cs001.output_data();

% TODO: show the data file(s)

% TODO: show where to retrieve the data from inside the object

%% Create a simple CellScan object without any interaction

% TODO: Point out that basically everything can be scripted

% Specify the filename for the image we want to load
fn001 = '.\tests\res\cellscan_gc_scim.tif';

% Specify what channels are contained in the image
ch001 = struct('Ca_Cyto_Astro', 1);

% Specify and load the calibration for the 20x objective
fnCal20x = '.\tests\res\calibration.mat';
cal20x = Calibration_PixelSize.load(fnCal20x);

% TODO: Discuss the calibration file... what it does/means etc

% Create the SCIM_Tif object with all this information
ri001 = SCIM_Tif(fn001, ch001, cal20x);

% Create a ConfigCellScan object from the two individual configs
configFRD = ConfigFindROIsDummy(); % finding the ROIs
configMRD = ConfigMeasureROIsDummy(); % measuring the ROIs
configCS002 = ConfigCellScan(configFRD, configMRD);

% Specify the name for this object (will use the SCIM_Tif name if empty)
nameCS002 = [];

% Create the CellScan object
cs002 = CellScan(nameCS002, ri001, configCS002);

% Process and plot the CellScan object
cs002.process();
cs002.plot();

%% Use ImageJ Rois

% TODO: Show how to create the ImageJ ROIs

% Create the CellScan object interactively
cs003 = CellScan([], ri001);

% Process and plot the CellScan object
cs003.process();
cs003.plot();

% Look at some alternative plots
cs003.plot('rois');
cs003.plot('traces', 'doWholeFrame', false);

%% Use FLIKA

% Create a ConfigCellScan object from the two individual configs
configFRF = ConfigFindROIsFLIKA.from_preset('ca_cyto_astro'); % finding the ROIs
configMRD = ConfigMeasureROIsDummy(); % measuring the ROIs
configCS004a = ConfigCellScan(configFRF, configMRD);

% Create the CellScan object
cs004a = CellScan([], ri001, configCS004a);

% Process and plot the CellScan object
cs004a.process();
cs004a.plot();

% Look at a different plot
cs004a.plot('images');

%% Use FLIKA with alternate parameters

% Create a ConfigCellScan object from the two individual configs
configFRF_mod = ConfigFindROIsFLIKA.from_preset('ca_cyto_astro', ...
    'threshold_constant', 9); 
configMRD = ConfigMeasureROIsDummy(); % measuring the ROIs
configCS004b = ConfigCellScan(configFRF_mod, configMRD);

% Create the CellScan object
cs004b = CellScan([], ri001, configCS004b);

% Process and plot the CellScan object
cs004b.process();
cs004b.plot();

%% Analyse the FRET ratio

% Load the raw image
ri002_raw = SCIM_Tif([], [], cal20x);
ri002 = copy(ri002_raw);

% Calculate the FRET ratio
ri002.ch_calc_ratio([2, 1])

% Create a ConfigCellScan object from the two individual configs
x_pix = 512;
y_pix = 512;
configFRD_ij = ConfigFindROIsDummy.from_ImageJ([], x_pix, y_pix);
configCS005 = ConfigCellScan(configFRD_ij, configMRD);

% Create the CellScan
cs005 = CellScan([], ri002, configCS005);

% Process and plot the CellScan object
cs005.process();
cs005.plot();

% TODO: What do people do about 0/NaN/inf pixels

%% Motion correct within an image

% Create the SCIM_Tif object with all this information
ri003a = copy(ri002);

% Do the motion correction
ri003a.motion_correct('doPlot', true);

% Create the CellScan object, and do the additional processing
cs006a = CellScan([], ri003a, configCS005);

% Process and plot the CellScan object
cs006a.process
cs006a.plot()

%% Motion correct to a reference image

% Create a copy of the FRET image
ri003b = copy(ri002_raw);

% Import a high resolution reference iamge
ri003_hires = SCIM_Tif();

% Extract a reference image
refImg = mean(ri003_hires.rawdata(:,:,2,:),4);
figure, imagesc(refImg), axis image, axis off, colormap(gray)

% Motion correct using the high resolution image
ri003b.motion_correct('refImg', refImg, 'ch', 2, 'doPlot', true);

% Look at the resulting image
ri003b.plot()

% Calculate the FRET ratio
ri003b.ch_calc_ratio([2, 1], 'FRET_ratio')

% Create the CellScan object, and do the additional processing
cs006b = CellScan([], ri003b, configCS005);

% Process and plot the CellScan object
cs006b.process
cs006b.plot()

%% Assemble a dummy reference image from many single images

% Load the raw image files (each with one frame)
riArray001_mixed = SCIM_Tif();

% Calculate the FRET ratio
riArray001_mixed.ch_calc_ratio([2, 1], 'FRET_ratio')

% Split the array into two 
riArray001_spot1 = riArray001_mixed(1:2:end);
riArray001_spot2 = riArray001_mixed(2:2:end);

% Combine the raw images files into a single object (with many frames)
ri004_spot1 = RawImg.cat_data(riArray001_spot1);
ri004_spot2 = RawImg.cat_data(riArray001_spot2);

% Create the CellScan objects, selecting the appropriate ImageJ ROIs
cs007_spot1 = CellScan([], ri004_spot1);
cs007_spot2 = CellScan([], ri004_spot2);

% Process and plot the CellScan objects
cs007_spot1.process()
cs007_spot1.plot()
cs007_spot2.process()
cs007_spot2.plot()

%% Create an array of CellScans interactively

% Create the array of CellScans (use FLIKA / Simple Measurement)
csArray001 = CellScan();

% Process the CellScan array
csArray001.process(true)
csArray001 = csArray001.process(true);

timeit(@() csArray001.process(true), 1) % around 144.99s per run serial
% TODO: Start the parallel processing pool
timeit(@() csArray001.process(true), 1) % around 87.84s per run parallel

% Plot the CellScan array
csArray001.plot()

% Plot one of the individual raw images
csArray001(3).rawImg.plot()

%% Process in parallel

% Create an array of Raw Images
riArray002 = repmat(ri001, 1, 10);


% Create an array of simple CellScans
csArray002_simple = CellScan([], riArray002, configCS002);

% Time the processing speed
timeit(@() csArray002_simple.process(true), 1) % around 9.94s per run serial
% TODO: Start the parallel processing pool
timeit(@() csArray002_simple.process(true), 1) % around 5.26s per run parallel


% TODO: don't run these, but describe the results, approximately

% Create an array of complex CellScans
csArray002_complex = CellScan([], riArray002, configCS004a);

% Time the processing speed
timeit(@() csArray002_complex.process(true), 1) % around 70.06s per run serial
% TODO: Start the parallel processing pool
timeit(@() csArray002_complex.process(true), 1) % around 43.42s per run parallel

%% Process a multichannel image by splitting the channels

% Load a multi channel image
ri005 = SCIM_Tif([], [], cal20x);

% Split the image into different channels
[ri005_ast, ri005_neu] = ri005.split1(3, [1, 1]);

% Create a CellScan object for each of the 
cs008_ast = CellScan([], ri005_ast, configCS004a);
cs008_neu_imgj = CellScan([], ri005_neu);

% Process and plot the CellScan objects
cs008_ast.process
cs008_ast.plot
cs008_neu_imgj.process
cs008_neu_imgj.plot

% % Create, process, plot a CellScan object for the neurons using FLIKA
% cs008_neu_flika = CellScan([], ri005_neu);
% cs008_neu_flika.process()
% cs008_neu_flika.plot()

%% MultiChImg with XSectScan

% Create a MultiChImg and add to it
mc001 = MultiChImg();
mc001.add()

% Process the MultiChImg
mc001 = mc001.process();

% Process the MultiChImg and plot the results
timeit(@() mc001.process(), 1) % around 30.26s per run serial
% 
timeit(@() mc001.process(), 1) % around 11.94s per run parallel

% Plot the MultiChImg
mc001.plot()

%% LineScanVel

lsv001 = LineScanVel();
lsv001.process()
lsv001.plot()
lsv001.plot('windows')

%% LineScanDiam

lsd001 = LineScanDiam();
lsd001.process()
lsd001.plot()
