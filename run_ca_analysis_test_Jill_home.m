%% Clear workspace
clearvars

% Questions to work on
% looping through multiple conditions?
% loop through spots
% loop through channels

% check with Kim Re: FLIKA parameters he uses

%% Information about your images

% file names to be saved (csv files)
filename1='test_peaks_26_04_2016.csv';  % peak data
filename2='test_traces_26_04_2016.csv'; % area under the curve data

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

plotMotion = 0; %Plot motion correction movie
doplots = 1; %Plots for each trial

% data cell array
data = [];
traces = [];
All_data=[];

% Load calibration file
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
        
        ImgArray = ImgArray.motion_correct('refImg',refImg, 'ch',2,'minCorr', 0.4,'doPlot',true);
        
        if plotMotion==1
            ImgArray(:).plot()
        end
        
        %% Create configuration and loop through both channels
        for iCh= 1:length(fieldnames(channel))
            Ch_names= fieldnames(channel);
            % for finding ROIs
            if ~isempty(find(Settings.AutoOverlay==1)) && iCh==1
                if strcmp(Ch_names{iCh},'Ca_Cyto_Astro')
                    % automated cytosolic astrocyte
                    configFind{iCh}= ConfigFindROIsFLIKA.from_preset('ca_cyto_astro',...
                        'baselineFrames', Settings.BL_frames,'sigmaXY', 2.9,...
                        'sigmaT', 0.5, 'threshold_constant', 7,...
                        'min_rise_time',0.1689, 'erosionRadius', 1.4457,...
                        'discardBorderROIs',true);
                elseif strcmp(Ch_names{iCh},'Ca_Memb_Astro')
                    % automated membrane astrocyte
                    configFind{iCh} = ConfigFindROIsFLIKA.from_preset('ca_memb_astro',...
                        'baselineFrames', Settings.BL_frames,'sigmaXY', 1,...
                        'max_rise_time',7,'minPuffArea',2, 'minPuffTime', 0.4,...
                        'dilateXY',2,'threshold_constant', 6,...
                        'dilationRadius', 1.0573,'discardBorderROIs',true);
                end
            elseif ~isempty(find(Settings.AutoOverlay==2)) && iCh==2
                % automated neuronal signals
                configFind{iCh} = ConfigFindROIsFLIKA.from_preset('ca_neuron',...
                    'baselineFrames', Settings.BL_frames,'sigmaXY', 2.9,...
                    'sigmaT', 0.5, 'max_rise_time',0.5,...
                    'dilateXY',2,'threshold_constant', 7,'erosionRadius',2);
            end
            
            
            if ~isempty(find(Settings.ManualOverlay==1))&& iCh==1
                % ImageJ ROIs
                x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
                configFind{iCh} = ConfigFindROIsDummy.from_ImageJ(fullfile(testRoot,'RoiSet.zip'), x_pix, y_pix, 1);
            elseif ~isempty(find(Settings.ManualOverlay==2))&& iCh==2
                % ImageJ ROIs
                x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
                configFind{iCh}= ConfigFindROIsDummy.from_ImageJ(fullfile(testRoot,'RoiSet.zip'), x_pix, y_pix, 1);
            end
            
            % Measurment configuration for peak classification
            configMeasure = ConfigMeasureROIsClsfy('baselineFrames', Settings.BL_frames, 'excludeNaNs', false,...
                    'thresholdSD_low', 3, 'thresholdSD_band', 3); %measure ROIs with peak classification
               
            
            % Combine the two configs
            configCellScan{iCh} = ConfigCellScan(configFind{iCh}, configMeasure);
        %end
            % Create CellScan objects
            CSArray{iCh} = CellScan(fnList, ImgArray, configCellScan{iCh}, iCh);
        end
            
            %% Process the images
            
            %CSArray_Ch1.process();
            CSArray.process([], 'find');
            
            if isfield(channel,'Ca_Cyto_Astro')
                % Use this function to combine the masks
                combine_masks(CSArray);
                
                %CSArray_combined.process();
                CSArray.process();
            end
            
            %% Make the debugging plots
            if doplots ==1
                CSArray.plot;
            end
            
            %% Output data
            % loop through all CellScans and concatenate data
            for itrial= 1:length(CSArray)
                for ipeaks = 1:length(CSArray(1,itrial).calcMeasureROIs.data.amplitude)
                    temp{ipeaks,1} = CSArray(1,itrial).calcMeasureROIs.data.amplitude{ipeaks};
                    temp{ipeaks,2} = CSArray(1,itrial).calcMeasureROIs.data.area{ipeaks};
                    temp{ipeaks,3} = CSArray(1,itrial).calcMeasureROIs.data.fullWidth{ipeaks};
                    temp{ipeaks,4} = CSArray(1,itrial).calcMeasureROIs.data.halfWidth{ipeaks};
                    temp{ipeaks,5} = CSArray(1,itrial).calcMeasureROIs.data.numPeaks{ipeaks};
                    temp{ipeaks,6} = CSArray(1,itrial).calcMeasureROIs.data.peakTime{ipeaks};
                    temp{ipeaks,7} = CSArray(1,itrial).calcMeasureROIs.data.peakType{ipeaks};
                    temp{ipeaks,8} = CSArray(1,itrial).calcMeasureROIs.data.prominence{ipeaks};
                    temp{ipeaks,9} = CSArray(1,itrial).calcMeasureROIs.data.ROIname{ipeaks};
                    temp{ipeaks,10} = CSArray(1,itrial).calcMeasureROIs.data.peakAUC{ipeaks};
                    temp{ipeaks,11} = strcat('Trial',num2str(itrial,'%02d'));
                    %temp{ipeaks,12} = Settings.NameConditions{iCond};
                    temp{ipeaks,13} = Settings.SpotIDs;
                    temp{ipeaks,14} = CurrentAnimal;
                    temp{ipeaks,15} = CSArray(1,itrial).channelToUse;
                    % what about method for finding ROIs?
                    
                    temptraces{1,ipeaks}= CSArray(1,itrial).calcMeasureROIs.data.tracesNorm(1:end,ipeaks);
                end
                data= vertcat(data,temp);
                traces = cat(1,traces, temptraces);
            end
        end
    end
%end

names = {'amplitude';'area';'fullWidth';'halfWidth';'numPeaks';'peakTime';'peakType';'prominence';...
    'ROIname';'peakAUC';'Trial';'Condition';'Spot';'Animal';'Channel'};

All_data= vertcat(names, data);

%% save as a CSV
cd(fullfile(Settings.MainDir, 'Results'));
% write date to created file
cell2csv(filename1, All_data);
cell2csv(filename2, traces);

