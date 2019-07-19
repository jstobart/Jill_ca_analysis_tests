close all; clear all;

%% Information about your images
Settings.ResultsFolder = 'E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results';  % folder where data should be saved

Settings.FileNames = {  % folder names where images and roiSet.zip is found
    'E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\T-series1',...
    'E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\T-series2',...
    'E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\T-series3',...  % etc.
    };

channel = struct('<blank>',1,'Ca_Neuron',2);  % can change 'blank' to any channel as this doesn't matter

% Load calibration file
calibration ='C:\JillsFiles\matlab\CalibrationFiles\calibration_25x_approx.mat';
CalFile = CalibrationPixelSize.load(calibration);

for iFile= 1:length(Settings.FileNames)
    
    % Get image paths
    testRoot =Settings.FilesNames{iFile};    
    expfiles = dir(fullfile(testRoot,'*.xml'));
    fnTempList = {expfiles(:).name};
    fnList = fullfile(testRoot, fnTempList);
    
    % Load data with BioFormats
    ImgArray =  BioFormats(fnList, channel, CalFile);
    
    %ImgArray.plot();
    %% Configs for Cell Scan
    
    % load ROI set for
    x_pix= size(ImgArray(1,1).rawdata,2); y_pix= size(ImgArray(1,1).rawdata,1);
    scaleF = 1;
    zipPath= fullfile(testRoot,'RoiSet.zip');
    findConf = ConfigFindROIsDummy.from_ImageJ(zipPath, x_pix, y_pix, scaleF);
    
    measureConf = ConfigMeasureROIsDummy();
    
    detectConf = ConfigDetectSigsClsfy('baselineFrames', 30,...
        'propagateNaNs', false,'excludeNaNs', false, 'lpWindowTime', 2, 'spFilterOrder', 2,...
        'spPassBandMin',0.1, 'spPassBandMax', 1, 'thresholdLP', 3,'thresholdSP', 5);
    
    % Combine the configs into a CellScan config
    configCS= ConfigCellScan(findConf, measureConf,detectConf); %
    
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




