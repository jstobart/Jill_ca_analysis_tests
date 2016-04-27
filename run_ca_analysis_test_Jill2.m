%% Clear workspace
clearvars

%% Questions to work on
% looping through multiple conditions?
% loop through spots

% improve classification parameters?

% output- add in channel, trial to the data table

%% Information about your images
%Settings.MainDir = 'E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f';
%Settings.MainDir = 'E:\Data\Two_Photon_Data\GCaMP_RCaMP\cyto_GCaMP6s';
Settings.MainDir = 'D:\Data\GCaMP_RCaMP\cyto_GCaMP6s';
%Settings.MainDir = 'D:\Data\GCaMP_RCaMP\Lck_GCaMP6f';

Settings.AnimalNames = {
    %'RG16',...
    'RG14',...
    };
Settings.ScoreSheetNames = {
   % 'RG16_Scoresheet_test.xls',...
    'RG14_Scoresheet_test.xls',...
    };
Settings.NameConditions = {'Nostim','Stim','shortstim'};
Settings.AutoOverlay = 1;  %0,1,2 channels for automated ROI selection
Settings.ManualOverlay = 2;

%channel = struct('Ca_Memb_Astro',1,'Ca_Neuron',2);
channel = struct('Ca_Cyto_Astro',1,'Ca_Neuron',2);
        
Classify = 1; %Classify peaks
plotMotion = 1; %Plot motion correction movie
doplots = 1; %Plots for each trial

% final data file name
SaveFiles{1,1} = fullfile(Settings.MainDir, 'Results', 'Ch1_tests.csv');
SaveFiles{1,2} = fullfile(Settings.MainDir, 'Results', 'Ch2_tests.csv');

%% Load calibration file
calibration ='D:\matlab\2p-img-analysis\tests\res\calibration.mat';
CalFile = Calibration_PixelSize.load(calibration);

%% load scoresheet
Settings.ScoreSheetPath = fullfile(Settings.MainDir,'Scoresheets',Settings.ScoreSheetNames);

numAnimals = length(Settings.AnimalNames);
for iAnimal = 1:numAnimals
    CurrentAnimal = Settings.AnimalNames{iAnimal};
    %read Scoresheet of current animal
    CurrentSheet = Settings.ScoreSheetPath{iAnimal};
    Settings = readScoresheet2(CurrentSheet, Settings);
    
    %% Give image paths
    for iPath= 1:length(Settings.LowresPath)
        testRoot =Settings.LowresPath{iPath};
        expfiles = dir(fullfile(testRoot,'lowres*'));
        fnTempList = {expfiles(:).name};
        fnList = fullfile(testRoot, fnTempList);
        
        %% Create an array of ScanImage Tiffs
        ImgArray =  SCIM_Tif(fnList, channel, CalFile);
        
        %% Run motion correction
        expfiles2 = dir(fullfile(testRoot,'highres*'));
        fnTempList2 = {expfiles2(:).name};
        RefImgName = cell2mat(fullfile(testRoot, fnTempList2));
        
        HighRes = SCIM_Tif(RefImgName,channel, CalFile);
        % Extract a reference image
        refImg = mean(HighRes.rawdata(:,:,2,:),4);
        %figure, imagesc(refImg), axis image, axis off, colormap(gray)
        
        ImgArray = ImgArray.motion_correct('refImg',refImg, 'ch',2,'minCorr', 0.4);%,'doPlot',true); 
      
        if plotMotion==1
            ImgArray(:).plot()
        end
        
        %% Create configuration
        % for finding ROIs
        if ~isempty(find(Settings.AutoOverlay==1))
            if isfield(channel,'Ca_Cyto_Astro')
                % automated membrane astrocyte
                configFind_Ch1 = ConfigFindROIsFLIKA.from_preset('ca_cyto_astro',...
                    'baselineFrames', Settings.BL_frames,'sigmaXY', 2.9,...
                    'sigmaT', 0.5, 'threshold_constant', 7,...
                    'min_rise_time',0.1689, 'erosionRadius', 1.4457,...
                    'discardBorderROIs',true);
            elseif isfield(channel,'Ca_Memb_Astro')
                % automated membrane astrocyte
                configFind_Ch1 = ConfigFindROIsFLIKA.from_preset('ca_memb_astro',...
                    'baselineFrames', Settings.BL_frames,'sigmaXY', 1,...
                    'max_rise_time',7,'minPuffArea',2, 'minPuffTime', 0.4,...
                    'dilateXY',2,'threshold_constant', 6,...
                    'dilationRadius', 1.0573,'discardBorderROIs',true);
            end
        elseif ~isempty(find(Settings.ManualOverlay==1))
            % ImageJ ROIs
            x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
            configFind_Ch1 = ConfigFindROIsDummy.from_ImageJ(fullfile(testRoot,'RoiSet.zip'), x_pix, y_pix);
        end
        if ~isempty(find(Settings.AutoOverlay==2))
            % automated neuronal signals
            configFind_Ch2 = ConfigFindROIsFLIKA.from_preset('ca_neuron',...
                'baselineFrames', Settings.BL_frames,'sigmaXY', 2.9,...
                'sigmaT', 0.5, 'max_rise_time',0.5,...
                'dilateXY',2,'threshold_constant', 7,'erosionRadius',2);
        elseif ~isempty(find(Settings.ManualOverlay==2))
            % ImageJ ROIs
            x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
            configFind_Ch2= ConfigFindROIsDummy.from_ImageJ(fullfile(testRoot,'RoiSet.zip'), x_pix, y_pix);
        end
        
        % Measurment configuration (with or without peak classification)
        if Classify ==1
            configMeasure = ConfigMeasureROIsClsfy('baselineFrames', Settings.BL_frames, 'excludeNaNs', false); %measure ROIs with peak classification
        else
            configMeasure = ConfigMeasureROIsDummy('propagateNaNs',true);
        end
        
        %% Combine the two configs
        configCellScan_Ch1 = ConfigCellScan(configFind_Ch1, configMeasure);
        configCellScan_Ch2 = ConfigCellScan(configFind_Ch2, configMeasure);
        
        %% Create CellScan objects
        
        CSArray_Ch1 = CellScan(fnList, ImgArray, configCellScan_Ch1, 1);
        %CSArray_Ch2 = CellScan(fnList, ImgArray, configCellScan_Ch2, 2);
        
        %% Process the images
        
        %CSArray_Ch2.process();
        Ch1_test =CSArray_Ch1.process();
        
        if isfield(channel,'Ca_Cyto_Astro')
            % Use this function to combine the masks
            CSArray_Ch1_combined= combine_masks(Ch1_test);
            %CSArray_Ch2 = utils.combine_masks(CSArray_Ch2);
            
%             configFindDummy = ConfigFindROIsDummy('roiMask', CSArray_Ch1_combined);
%             calcFindDummy = configFindDummy.create_calc();
%             
%             
%             [CSArray_Ch1(:).calcFindROIs] = deal(calcFindDummy);
            
            CSArray_Ch1.process();
        end
        
        %% Make the debugging plots
        if doplots ==1
            CSArray_Ch1.plot;
           % CSArray_Ch2.plot;
        end
        
        %% Output data
        CSArray_Ch1.output_data();%([],SaveFiles{1,1});
        %CSArray_Ch2.output_data([],Files{1,2});
    end
end

