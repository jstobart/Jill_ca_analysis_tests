close all; clear all;

%% Information about your images

% folder where data should be saved for each animal
Settings.ResultsFolder = 'D:\Data\test';

Settings.FileNames = {  % folder names where images and roiSet.zip is found
    'J:\Jill_Stobart\In_vivo_2P_Data\66678_Crazy8\2019_06_14\0Hz\spot1-frequenices-TSeries-06142019-1038-009',...
    'J:\Jill_Stobart\In_vivo_2P_Data\66678_Crazy8\2019_06_14\40Hz\spot1-frequenices-TSeries-06142019-1038-012',...
     % etc.
    };

Settings.Baseline = 2; % time (s) before the whisker stimulator starts,  2s for short trials (10s long), 5s for long trials

channel = struct('Ca_Neuron',1,'Ca_Memb_Astro',2);  % can change 'blank' to any channel as this doesn't matter

%% Loop through each file and make Cell Scans for LckGCaMP and RCaMP

for iFile= 1:length(Settings.FileNames)
    
    % Get the T-series path from the folder
    testRoot =Settings.FileNames{iFile};
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
     % Pericyte RCaMP
    % 2D automated selection for peaks
    findConf1 = ConfigFindROIsFLIKA_2D.from_preset('ca_neuron', 'baselineFrames',...
        BL_frames,'freqPassBand',0.5,'sigmaXY', 2,...
        'sigmaT', 0.14,'thresholdPuff', 7, 'threshold2D', 0.1,...
        'minRiseTime',0.14, 'maxRiseTime', 1,'minROIArea', 10,...
        'dilateXY', 4, 'dilateT', 0.3,'erodeXY', 1, 'erodeT', 0.1);
    
    % hand selected- peaks from cellular structures
    x_pix= size(ImgArray(1,1).rawdata,2); y_pix= size(ImgArray(1,1).rawdata,1);
    scaleF = 1;
    
    zipfiles = dir(fullfile(testRoot,'*.zip'));  % look for the zip ROI file
    fnTempList2 = {zipfiles(:).name};
    zipPath = fullfile(testRoot, fnTempList2);
       
    % load image J ROIs for mask
    findConf2 = ConfigFindROIsDummy.from_ImageJ(zipPath{1,1}, x_pix, y_pix, scaleF);
    
    
    % measure ROIs (extract the traces)
        measureConf = ConfigMeasureROIsDummy('baselineFrames', BL_frames);
        
    % filter the traces to detect the peaks and get info about them    
        detectConf = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...
        'propagateNaNs', false, 'excludeNaNs', false, 'lpWindowTime', 1.5, 'spFilterOrder', 2,...
        'spPassBandMin',0.05, 'spPassBandMax', 0.5, 'thresholdLP', 3,'thresholdSP', 5);
    
    % Combine the configs into a CellScan config for automated ROI
    % selection
    configCS1= ConfigCellScan(findConf1, measureConf, detectConf); %
    
    % Combine the configs into a CellScan config for ImageJ ROIs
    configCS2= ConfigCellScan(findConf2, measureConf, detectConf); %
    
    
    %% Create CellScan objects
    automatedROIs = CellScan(fnList, ImgArray, configCS1, 1); % peaks from automated ROIs
    ImageJROIs = CellScan(fnList, ImgArray, configCS2, 1); % peaks from hand clicked ROIs
    
    % Process the images
    automatedROIs =automatedROIs.process();
    ImageJROIs =ImageJROIs.process();
    
    % Make the debugging plots
    automatedROIs.plot();
    ImageJROIs.plot();
    
    %automatedROIs.opt_config()
    %ImageJROIs.opt_config()
    %ImageJROIs.plot('signals');
    
    %% Output data
    OutputFile1=fullfile(Settings.ResultsFolder, strcat('Astrocytes_FieldofView',num2str(iFile)));
    OutputFile2=fullfile(Settings.ResultsFolder, strcat('Neurons_FieldofView',num2str(iFile)));
    
    automatedROIs.output_data(OutputFile1)
    ImageJROIs.output_data(OutputFile2)
    
end




