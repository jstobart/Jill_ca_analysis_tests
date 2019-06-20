close all; clear all;

%% Information about your images

% folder where data should be saved
Settings.ResultsFolder = 'E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results';

Settings.FileNames = {  % folder names where images and roiSet.zip is found
    'E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\T-series1',...
    'E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\T-series2',...
    'E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\T-series3',...  % etc.
    };

Settings.Baseline = 2; % time (s) before the whisker stimulator starts,  2s for short trials (10s long), 5s for long trials

channel = struct('Ca_Neuron',1,'Ca_Membr_Astro',2);  % can change 'blank' to any channel as this doesn't matter

%% Loop through each file and make Cell Scans for LckGCaMP and RCaMP

for iFile= 1:length(Settings.FileNames)
    
    % Get the T-series path from the folder
    testRoot =Settings.FilesNames{iFile};
    expfiles = dir(fullfile(testRoot,'*.xml'));  % look for the xml file
    fnTempList = {expfiles(:).name};
    fnList = fullfile(testRoot, fnTempList);
    
    % Load data with BioFormats
    ImgArray =  BioFormats(fnList, channel);
    
    %% Motion Correction
    refImg=mean(ImgArray.rawdata(:,:,1,1:3),4);
    ImgArray = ImgArray.motion_correct('refImg',refImg, 'ch',1,'maxShift', 10,'minCorr', 0.4);%,'doPlot',true);
    
    
    BL_frames= Settings.Baseline*ImgArray.metadata.frameRate; % number of baseline frames
    
    %% Configs for Lck GCaMP Cell Scan
     %ASTROCYTES
    % 2D automated selection for peaks
    AC_findConf = ConfigFindROIsFLIKA_2D.from_preset('ca_memb_astro', 'baselineFrames',...
        BL_frames,'freqPassBand',0.5,'sigmaXY', 2,...
        'sigmaT', 0.14,'thresholdPuff', 7, 'threshold2D', 0.1,...
        'minRiseTime',0.14, 'maxRiseTime', 1,'minROIArea', 10,...
        'dilateXY', 4, 'dilateT', 0.3,'erodeXY', 1, 'erodeT', 0.1);
    
    % measure ROIs (extract the traces)
        AC_measureConf = ConfigMeasureROIsDummy();
        
    % filter the traces to detect the peaks and get info about them    
        AC_detectConf = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...
        'propagateNaNs', false, 'excludeNaNs', false, 'lpWindowTime', 1.5, 'spFilterOrder', 2,...
        'spPassBandMin',0.05, 'spPassBandMax', 0.5, 'thresholdLP', 3,'thresholdSP', 5);
    
    % Combine the configs into a CellScan config
    AC_configCS= ConfigCellScan(AC_findConf, AC_measureConf, AC_detectConf); %
    
    %% Configs for neuronal RCaMP Cell Scan
    
    % hand selected- peaks from cellular structures
    x_pix= size(ImgArray(1,1).rawdata,2); y_pix= size(ImgArray(1,1).rawdata,1);
    scaleF = 1;
    
    zipPath= fullfile(testRoot,'RoiSet.zip');
    
    findConf = ConfigFindROIsDummy.from_ImageJ(zipPath, x_pix, y_pix, scaleF);
    
    AC_findConf{2} = ConfigFindROIsDummy.from_ImageJ(fullfile(testRoot,'Astrocytes.zip'), x_pix, y_pix,1);
    
    
   
    % NEURONS
    % 2D FLIKA selected for peaks from "dendrites"
    Neur_findConf{1} = ConfigFindROIsFLIKA_2D.from_preset('ca_neuron', 'baselineFrames',...
        BL_frames,'freqPassBand',1,'sigmaXY', 2,...
        'sigmaT', 0.1,'thresholdPuff', 7, 'threshold2D', 0.2,...
        'minRiseTime',0.0845, 'maxRiseTime', 2,'minROIArea', 10,...
        'dilateXY', 3, 'dilateT', 0.5,'erodeXY', 2, 'erodeT', 0.3,...
        'discardBorderROIs',true);
    
    % hand selected for peaks from somata
    x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
    Neur_findConf{2} = ConfigFindROIsDummy.from_ImageJ(fullfile(testRoot,'Neurons.zip'), x_pix, y_pix,1);
    
    % 3D FLIKA selected for time and space estimations
    Neur_findConf{3} = ConfigFindROIsFLIKA_3D.from_preset('ca_neuron', 'baselineFrames',...
        BL_frames,'freqPassBand',1,'sigmaXY', 2,...
        'sigmaT', 0.1,'thresholdPuff', 7, 'threshold2D', 0.2,...
        'minRiseTime',0.0845, 'maxRiseTime', 2,'minROIArea', 10,...
        'dilateXY', 3, 'dilateT', 0.5,'erodeXY', 2, 'erodeT', 0.3,...
        'discardBorderROIs',true);
    
    %% Configuration for measuring ROIs
    % AWAKE astrocyte membrane calcium

    
    % AWAKE neuron calcium
    detectConf{2} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,... 'normMethod','z-score',... 'zIters', 10000,...
        'propagateNaNs', false,'excludeNaNs', false, 'lpWindowTime', 2, 'spFilterOrder', 2,...
        'spPassBandMin',0.1, 'spPassBandMax', 1, 'thresholdLP', 3,'thresholdSP', 5);
    
    
    %% Configs for RCaMP Cell Scan
    % load ROI set for
    x_pix= size(ImgArray(1,1).rawdata,2); y_pix= size(ImgArray(1,1).rawdata,1);
    scaleF = 1;
    zipPath= fullfile(testRoot,'RoiSet.zip');
    findConf = ConfigFindROIsDummy.from_ImageJ(zipPath, x_pix, y_pix, scaleF);
    
    measureConf = ConfigMeasureROIsDummy();
    
    detectConf = ConfigDetectSigsClsfy('baselineFrames', 30,...
        'propagateNaNs', false,'excludeNaNs', false, 'lpWindowTime', 2, 'spFilterOrder', 2,...
        'spPassBandMin',0.1, 'spPassBandMax', 1, 'thresholdLP', 3,'thresholdSP', 5);
    
    
    
    %% Create CellScan object
    neurons = CellScan(fnList, ImgArray, configCS, 2); % peaks from hand clicked ROIs
    
    % Process the images
    neurons =neurons.process();
    
    % Make the debugging plots
    neurons.plot();
    %neurons.opt_config()
    
    
    %% Output data
    OutputFile=fullfile(Settings.ResultsFolder, strcat('FieldofView',num2str(iFile)));
    neurons.output_data(OutputFile)
    
end




