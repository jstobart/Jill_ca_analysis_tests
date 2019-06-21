close all; clear all;

%% Information about your images

% folder where data should be saved for each animal
Settings.ResultsFolder = 'D:\Data\test';

Settings.FileNames = {  % folder names where images and roiSet.zip is found
    'J:\Jill_Stobart\In_vivo_2P_Data\66678_Crazy8\2019_06_14\spot1',...
    % etc.
    };

Settings.Baseline = 2; % time (s) before the whisker stimulator starts,  2s for short trials (10s long), 5s for long trials

channel = struct('Ca_Neuron',1,'Ca_Memb_Astro',2);  % can change 'blank' to any channel as this doesn't matter

%% Loop through each file and make Cell Scans for LckGCaMP and RCaMP

for iSpot= 1:length(Settings.FileNames)
    
    SpotRoot= Settings.FileNames{iSpot};
    
    % Get a list of all files and folders in this folder.
    % different stimulation conditions
    ConditionFolders = dir(SpotRoot);  % look for folders
    ConditionFolders(ismember( {ConditionFolders.name}, {'.', '..'})) = [];
    % Get a logical vector that tells which is a directory.
    names    = {ConditionFolders.name};
    dirFlags = [ConditionFolders.isdir];
    % Extract only those that are directories.
    subDirsNames = names(dirFlags);
    
    for iCondition= 1:length(subDirsNames)
        FolderName=fullfile(SpotRoot, subDirsNames(iCondition));
        
        % Get a list of all T-series folders in this folder.
        % different trials
        TrialFolders = dir(FolderName{1,1});  % look for folders
        TrialFolders(ismember( {TrialFolders.name}, {'.', '..'})) = [];
        % Get a logical vector that tells which is a directory.
        names2    = {TrialFolders.name};
        dirFlags2 = [TrialFolders.isdir];
        % Extract only those that are directories.
        TrialNames = names2(dirFlags2);
        
        for iTrial= 1:length(TrialNames)
            % Get the T-series path from the folder
            testRoot =fullfile(FolderName, TrialNames{iTrial});
            
            expfiles = dir(fullfile(testRoot{1,1},'*.xml'));  % look for the xml file
            fnList = fullfile(testRoot{1,1}, expfiles.name);
            
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
            
            % Combine the configs into a CellScan config for membrane tagged GCaMP
            AC_configCS= ConfigCellScan(AC_findConf, AC_measureConf, AC_detectConf); %
            
            %% Configs for neuronal RCaMP Cell Scan
            % NEURONS
            % hand selected- peaks from cellular structures
            x_pix= size(ImgArray(1,1).rawdata,2); y_pix= size(ImgArray(1,1).rawdata,1);
            scaleF = 1;
            
            zipfiles = dir(fullfile(testRoot,'*.zip'));  % look for the zip ROI file
            fnTempList2 = {zipfiles(:).name};
            zipPath = fullfile(testRoot, fnTempList2);
            
            % load image J ROIs for mask
            N_findConf = ConfigFindROIsDummy.from_ImageJ(zipPath{1,1}, x_pix, y_pix, scaleF);
            
            % measure ROIs (extract the traces)
            N_measureConf = ConfigMeasureROIsDummy();
            
            % filter the traces to detect the peaks and get info about them
            N_detectConf = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...
                'propagateNaNs', false,'excludeNaNs', false, 'lpWindowTime', 5, 'spFilterOrder', 2,...
                'spPassBandMin',0.1, 'spPassBandMax', 1, 'thresholdLP', 5,'thresholdSP', 4);
            
            
            % Combine the configs into a CellScan config for neuronal RCaMP
            N_configCS= ConfigCellScan(N_findConf, N_measureConf, N_detectConf); %
            
            
            %% Create CellScan objects
            astrocytes = CellScan(fnList, ImgArray, AC_configCS, 2); % peaks from automated ROIs
            neurons = CellScan(fnList, ImgArray, N_configCS, 1); % peaks from hand clicked ROIs
            
            % Process the images
            astrocytes =astrocytes.process();
            neurons =neurons.process();
            
            % Make the debugging plots
            astrocytes.plot();
            neurons.plot();
            
            %astrocytes.opt_config()
            %neurons.opt_config()
            %neurons.plot('signals');
            
            %% Output data
            OutputFile1=fullfile(Settings.ResultsFolder, strcat('Astrocytes_FieldofView',num2str(iFile)));
            OutputFile2=fullfile(Settings.ResultsFolder, strcat('Neurons_FieldofView',num2str(iFile)));
            
            astrocytes.output_data(OutputFile1)
            neurons.output_data(OutputFile2)
            
            % output the different levels of data
            % output traces
            % output ROI masks
            % output onsets
            % output peaks
            % csv and matlab files.....
            
        end
    end
end




