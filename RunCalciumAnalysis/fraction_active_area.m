%% Load calibration
% Load the calibration for your microscope
% If you do not have a calibration, check 'doc CalibrationPixelSize' and
% the referenced quick start guide
fNameCal = fullfile('P:\_Group\Software\CHIPS/calibration_20x.mat');
calibration = CalibrationPixelSize.load(fNameCal);

%% Give image path(s)
% This example looks at a single TIFF stack. To demonstrate array usage, we
% use the same image 3 times. Your fnList can either be an empty string
% (MATLAB then prompts you to choose a file manually) a single string or a
% cell array of strings as below.

% Download example images
% utils.download_example_imgs();
fImage = 'cellscan_scim.tif';
fnImage = fullfile(utils.CHIPS_rootdir, 'tests\res\', fImage);
fnList = repmat({fnImage}, 3, 1);  

%% Create channels structure
% We need to tell the program, what is on the imaging channels. Currently,
% FLIKA supports a total of three different channels, as described below.
% You can specify multiple channels but make sure to also specify the index
% of the channel that should be used for FLIKA in the channelToUse
% variable.
%
% SYNTAX: channels = struct('Channel_1_Name', channel_1_Index, ...
%                           'Channel_n_Name, channel_n_Index);
%
% known ChannelNames for FLIKA are :    Ca_Cyto_Astro 
%                                       Ca_Memb_Astro
%                                       Ca_Neuron
%                                       FRET_ratio
%                                       cellular_signal
%
% channelIndex refers to channel number in the tiff file

channels = struct('Ca_Memb_Astro', 1);
channelToUse = 1;

%% Create an array of ScanImage Tiffs
% The following lines will load your image files and create an array of
% SCIM_Tif objects. Inside the objects, all your raw data and metadata is 
% stored.
%
% SYNTAX: filename, channels, calibration, skipImport
%
% filename can either be empty, a string or a cell array of strings
% 
img =  SCIM_Tif(fnList, channels, calibration);

%% OPTIONAL pre-processing steps

%Run motion correction
% You can use this step to apply a 2D FFT-based or line-by-line HMM motion
% correction to your specified channel before continuing.
%
% SYNTAX: imgArray.motion_correct(..., 'attribute', value)
%
% The default method so far is convFFT. You can provide a reference image
% through refImg if you want. This can be a high resolution or low
% resolution image. The program will attempt to rescale it by itself. If
% refImg is not specified the first frame of the image stack will be used
% as reference. doPlot is a logical. Set it to true if you need a plot to
% visualize motion. For all possible options, check the help page for
% utils.motion_correct().
% img = img.motion_correct('ch', channelToUse, 'method', 'convFFT', ...
%     'fillBadData', 'inpaint');
% 
%Exclude frames from image 
% This step can be used to exlude certain frames from the stack, for
% example when you have a flash from optogenetic stimulation or bad first
% frames due to a bad Pockel's cell.
%
% SYNTAX: imgArray.exclude_frames(..., 'attribute', value)
%
% See the help page for RawImg.exclude_frames() for more details.
% img.exclude_frames('badFrames', 1, 'method', 'inpaint')

%% Create config
configFind = ConfigFindROIsFLIKA_2D();
configMeasure = ConfigMeasureROIsDummy();
configDetect = ConfigDetectSigsDummy();

configCS = ConfigCellScan(configFind, configMeasure, configDetect);

%% Create CellScan objects
% Finally, the interesting part! This is how you create CellScans.
% SYNTAX: CSArray = CellScan(name, rawImg, configIn, channelToUseIn)
CSArray = CellScan('', img, configCS, channelToUse);

%% Process the images
% You can use CellScan.process('flag') either with flags 'find', 'measure'
% and 'detect' to run the three steps separately or run CellScan.process()
% without any input to run all steps. 

% Process the CellScan
CSArray = CSArray.process;

%% Find threshold value
% Create mean baseline image
newimg = mean(CSArray(1).rawImg.rawdata,4); 

% Create figure, show baseline image and hold on
hFig = figure;
imagesc(newimg);
hold on

% Ask user for new threshold values and upate figure, until user stops
% execution
while 1
    prompt = 'What is the threshold value (type 0 to finish)?\n';
    x = input(prompt);
    
    if x == 0
        break
    end
    
    tempimg = newimg;
    fluoidx = newimg > x;
    tempimg(~fluoidx) = max(newimg(:));
    imagesc(tempimg)
end
