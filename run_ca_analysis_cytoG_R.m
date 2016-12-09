%% Clear workspace
clearvars

AllData= [];
All_AUC= [];
%% Questions to work on

% spectral unmixing?

% improve classification parameters?

% exclude overlapping ROIs?

%% Information about your images
Settings.MainDir = 'E:\Data\Two_Photon_Data\GCaMP_RCaMP\cyto_GCaMP6s';

Settings.AnimalNames = {
    %'RG16',...
    'RG14',...
    };
Settings.ScoreSheetNames = {
    % 'RG16_Scoresheet_test.xls',...
    'RG14_Scoresheet_test.xls',...
    };
Settings.NameConditions = {'Nostim','Stim'};

channel = struct('Ca_Cyto_Astro',1,'Ca_Neuron',2);

plotMotion = 0; %Plot motion correction movie
doplots = 0; %Plots for each trial

% final data file name
SaveFiles{1,1} = fullfile(Settings.MainDir, 'Results', 'Tests.csv'); % all data
SaveFiles{1,2}= 'CellScan_AC_FLIKA.mat'; % astrocyte FLIKA cell scan
SaveFiles{1,3}= 'CellScan_AC_Hand.mat'; % astrocyte hand click cell scan
SaveFiles{1,4}= 'CellScan_Ne_Hand.mat'; % neuronal hand click cell scan
SaveFiles{1,5}= fullfile(Settings.MainDir, 'Results', 'TraceAUC_10s.csv'); % neuronal hand click cell scan

%% Load calibration file
calibration ='E:\matlab\2p-img-analysis\tests\res\calibration_20x.mat';
CalFile = Calibration_PixelSize.load(calibration);

%% load scoresheet and loop through animal, spot, etc.

Settings.ScoreSheetPath = fullfile(Settings.MainDir,'Scoresheets',Settings.ScoreSheetNames);

numAnimals = length(Settings.AnimalNames);
for iAnimal = 1:numAnimals
    CurrentAnimal = Settings.AnimalNames{iAnimal};
    %read Scoresheet of current animal
    CurrentSheet = Settings.ScoreSheetPath{iAnimal};
    Settings = readScoresheet2(CurrentSheet, Settings);
    
    % Get SpotID
    spots = unique(Settings.SpotIDs);
    
    for iSpot = 1:length(spots) % looks at every 2nd line of SpotIDs-
        
        % Extract spot name
        spotId = spots{iSpot};
        
        % Find the idx of paths matching this spot (all Conditions)
        matchingSpotsIdx = find(~cellfun(@isempty, regexp(Settings.SpotIDs, spotId)));
        SpotPaths = Settings.LowresPath(matchingSpotsIdx);
        CurrentDepth = Settings.Depth(matchingSpotsIdx);
        CurrentBaseline = Settings.BL_frames(matchingSpotsIdx);
        
        if length(SpotPaths) ~= length(Settings.NameConditions)
            error('not all stimulus conditions are present')
        end
        
        for iCond = 1:length(Settings.NameConditions)
            CurrentCondition = Settings.NameConditions{iCond};
            PathIdx = find(~cellfun(@isempty, regexp(SpotPaths, CurrentCondition)));
            BL_frames = CurrentBaseline(PathIdx);
            
            % Get image paths
            testRoot =SpotPaths{PathIdx};
            
            if ~exist(fullfile(testRoot,SaveFiles{1,2}),'file') %if the CellScan does not exist, create a new one
                expfiles = dir(fullfile(testRoot,'lowres*'));
                fnTempList = {expfiles(:).name};
                fnList = fullfile(testRoot, fnTempList);
                
                % Create an array of ScanImage Tiffs
                ImgArray =  SCIM_Tif(fnList, channel, CalFile);
                
                % Spectral Unmixing of GCaMP and RCaMP
                %ImgArray= ImgArray.unmix_ch;
                
                % Run motion correction
                expfiles2 = dir(fullfile(testRoot,'highres*'));
                fnTempList2 = {expfiles2(:).name};
                RefImgName = cell2mat(fullfile(testRoot, fnTempList2));
                
                HighRes = SCIM_Tif(RefImgName,channel, CalFile);
                % Extract a reference image
                refImg = mean(HighRes.rawdata(:,:,2,:),4);
                %figure, imagesc(refImg), axis image, axis off, colormap(gray)
                
                ImgArray = ImgArray.motion_correct('refImg',refImg, 'ch',2,'maxShift', 10,'minCorr', 0.4);%,'doPlot',true);
                % fill in bad data?
                if plotMotion==1
                    ImgArray(1,1).plot()
                end
                
                %% Configs for Finding ROIs
                %Astrocyte Calcium
                % automated selection
                findConf{1} = ConfigFindROIsFLIKA_2D.from_preset('ca_cyto_astro', 'baselineFrames',...
                    BL_frames,'sigmaXY', 2,...
                    'sigmaT', 0.5,'threshold_std', 5,...
                    'min_rise_time',0.1689, 'max_rise_time', 8,...
                    'dilateXY', 4,... %'minPuffArea', 5,...
                    'dilateT', 0.5,'erodeXY', 2,...
                    'inpaintIters', 5, 'discardBorderROIs',false);
                
                % Astrocyte Calcium
                % hand selected- cellular structures
                x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
                findConf{2} = ConfigFindROIsDummy.from_ImageJ(fullfile(testRoot,'Astrocytes.zip'), x_pix, y_pix,1);
                
                % Neuronal Calcium
                % hand selected
                x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
                findConf{3} = ConfigFindROIsDummy.from_ImageJ(fullfile(testRoot,'Neurons.zip'), x_pix, y_pix,1);
                
                %% Configuration for measuring ROIs
                % AWAKE astrocyte cyto calcium
                measureConf{1} = ConfigMeasureROIsClsfy('baselineFrames', BL_frames,...'normMethod', 'z-score',...
                    'propagateNaNs', false, 'excludeNaNs', false, 'lpWindowTime', 5, 'spFilterOrder', 4,...
                    'spPassBand', [0.025, 0.2], 'thresholdSD_low', 3,'thresholdSD_band', 5);
                
                % AWAKE neuron calcium
                measureConf{2} = ConfigMeasureROIsClsfy('baselineFrames', BL_frames,...  'normMethod','z-score','zIters', 10000,...
                    'propagateNaNs', false,'excludeNaNs', false, 'lpWindowTime', 3, 'spFilterOrder', 2,...
                    'spPassBand', [0.1, 2], 'thresholdSD_low', 7,'thresholdSD_band', 7);
                
                % for calculating AUC for each trace
                measureConf{3} = ConfigMeasureROIsDummy('baselineFrames', BL_frames);
                
                %%
                % Combine the configs into a CellScan config
                configCS{1,1} = ConfigCellScan(findConf{1,1}, measureConf{1,1}); % astrocyte FLIKA, peaks
                configCS{1,2} = ConfigCellScan(findConf{1,2}, measureConf{1,1}); % astrocyte hand, peaks
                configCS{1,3} = ConfigCellScan(findConf{1,3}, measureConf{1,2}); % neuronal hand, peaks
                configCS{1,4} = ConfigCellScan(findConf{1,1}, measureConf{1,3}); % astrocyte FLIKA, AUC
                configCS{1,5} = ConfigCellScan(findConf{1,2}, measureConf{1,3}); % astrocyte hand, AUC
                configCS{1,6} = ConfigCellScan(findConf{1,3}, measureConf{1,3}); % neuronal hand, AUC
                
                
                %% Create CellScan objects
                CSArray_Ch1_FLIKA1 = CellScan(fnList, ImgArray, configCS{1,1}, 1);
                CSArray_Ch1_FLIKA2 = CellScan(fnList, ImgArray, configCS{1,4}, 1);
                
                CSArray_Ch1_Hand1 = CellScan(fnList, ImgArray, configCS{1,2}, 1);
                CSArray_Ch1_Hand2 = CellScan(fnList, ImgArray, configCS{1,5}, 1);
                
                CSArray_Ch2_Hand1 = CellScan(fnList, ImgArray, configCS{1,3}, 2);
                CSArray_Ch2_Hand2 = CellScan(fnList, ImgArray, configCS{1,6}, 2);
                
                %% Process the images
                CSArray_Ch1_FLIKA1 =CSArray_Ch1_FLIKA1.process();
                CSArray_Ch1_Hand1 =CSArray_Ch1_Hand1.process();
                CSArray_Ch2_Hand1 =CSArray_Ch2_Hand1.process();
                
                CSArray_Ch1_FLIKA2 =CSArray_Ch1_FLIKA2.process();
                CSArray_Ch1_Hand2 =CSArray_Ch1_Hand2.process();
                CSArray_Ch2_Hand2 =CSArray_Ch2_Hand2.process();
                
                %NO IMAGE COMBINATION
                % if isfield(channel,'Ca_Cyto_Astro')
                % Use this function to combine the masks
                %CSArray_Ch1_combined= combine_masks(Ch1_test);
                %CSArray_Ch1_combined.process();
                %end
                
                %% Make the debugging plots
                if doplots ==1
                    CSArray_Ch1_FLIKA1.plot;
                    CSArray_Ch1_Hand1.plot;
                    CSArray_Ch2_Hand1.plot;
                    
                    %             CSArray_Ch1_FLIKA1.plot('peaks');
                    %             CSArray_Ch1_Hand1.plot('peaks');
                    %            CSArray_Ch2_Hand1.plot('peaks');
                end
            else
                load (fullfile(testRoot,SaveFiles{1,2}));
                load (fullfile(testRoot,SaveFiles{1,3}));
                load (fullfile(testRoot,SaveFiles{1,4}));
            end
            
            %% Output data
            
            % make a giant data table
            
            listFields = {'amplitude', 'area', 'fullWidth', 'halfWidth', ...
                'numPeaks', 'peakTime', 'peakStart', 'peakStartHalf', ...
                'peakType', 'prominence', 'ROIname', 'peakAUC'};
            
            % Astrocyte FLIKA
            for itrial=1:length(CSArray_Ch1_FLIKA1)
                % peak output
                temp=CSArray_Ch1_FLIKA1(1,itrial).calcMeasureROIs.data;
                temp2.trialname ={};
                temp2.animalname = {};
                temp2.channel = {};
                temp2.Spot = {};
                temp2.Cond = {};
                temp2.depth = {};
                
                % extract fields from Class
                for jField = 1:numel(listFields)
                    isFirst = (itrial == 1 );
                    if isFirst
                        data.(listFields{jField}) = {};
                    end
                    data.(listFields{jField}) = [data.(listFields{jField}); ...
                        temp.(listFields{jField})];
                end
                
                % create fields for trial, animal, spot, condition, etc.
                for iData = 1:length(temp.amplitude)
                    temp2.trialname{iData,1}=strcat('trial', num2str(itrial));
                    temp2.channel{iData,1}= 'GCaMP';
                    temp2.Spot{iData,1}= spotId;
                    temp2.animalname{iData,1}= CurrentAnimal;
                    temp2.Cond{iData,1} = CurrentCondition;
                    temp2.depth{iData,1} = CurrentDepth(1,1);
                end
                
                isFirst = (itrial == 1 );
                if isFirst
                    data.Trial = {};
                    data.Animal = {};
                    data.Channel = {};
                    data.Spot = {};
                    data.Condition = {};
                    data.Depth = {};
                end
                data.Trial= [data.Trial; temp2.trialname];
                data.Animal= [data.Animal; temp2.animalname];
                data.Channel= [data.Channel; temp2.channel];
                data.Spot= [data.Spot; temp2.Spot];
                data.Condition= [data.Condition; temp2.Cond];
                data.Depth= [data.Depth; temp2.depth];
                
                % trace AUC output
                traces= CSArray_Ch1_FLIKA2(1,itrial).calcMeasureROIs.data.tracesNorm;
                for iROI = 1:size(traces,2)
                    AUCdata{iROI,1}= trapz(traces(60:60+round(10*11.84),iROI));
                    AUCdata{iROI,2}=CSArray_Ch1_FLIKA2(1,itrial).calcFindROIs.data.roiNames{iROI,1};
                    AUCdata{iROI,3}=strcat('trial', num2str(itrial));
                    AUCdata{iROI,4}= 'GCaMP';
                    AUCdata{iROI,5}= spotId;
                    AUCdata{iROI,6}= CurrentAnimal;
                    AUCdata{iROI,7} = CurrentCondition;
                    AUCdata{iROI,8} = CurrentDepth(1,1);
                end
                All_AUC=vertcat(All_AUC, AUCdata);
                clearvars AUCdata
            end
            clearvars temp temp2 
            
            % Astrocyte Hand clicked ROIs
            for itrial=1:length(CSArray_Ch1_Hand1)
                
                temp=CSArray_Ch1_Hand1(1,itrial).calcMeasureROIs.data;
                temp2.trialname ={};
                temp2.animalname = {};
                temp2.channel = {};
                temp2.Spot = {};
                temp2.Cond = {};
                temp2.depth = {};
                
                % extract fields from Class
                for jField = 1:numel(listFields)
                    data.(listFields{jField}) = [data.(listFields{jField}); ...
                        temp.(listFields{jField})];
                end
                
                % create fields for trial, animal, spot, condition, etc.
                for iData = 1:length(temp.amplitude)
                    temp2.trialname{iData,1}=strcat('trial', num2str(itrial));
                    temp2.channel{iData,1}= 'GCaMP';
                    temp2.Spot{iData,1}= spotId;
                    temp2.animalname{iData,1}= CurrentAnimal;
                    temp2.Cond{iData,1} = CurrentCondition;
                    temp2.depth{iData,1} = CurrentDepth(1,1);
                end
                
                data.Trial= [data.Trial; temp2.trialname];
                data.Animal= [data.Animal; temp2.animalname];
                data.Channel= [data.Channel; temp2.channel];
                data.Spot= [data.Spot; temp2.Spot];
                data.Condition= [data.Condition; temp2.Cond];
                data.Depth= [data.Depth; temp2.depth];
                
                % trace AUC output
                traces= CSArray_Ch1_FLIKA2(1,itrial).calcMeasureROIs.data.tracesNorm;
                for iROI = 1:size(traces,2)
                    AUCdata{iROI,1}= trapz(traces(60:60+round(10*11.84),iROI));
                    AUCdata{iROI,2}=CSArray_Ch1_FLIKA2(1,itrial).calcFindROIs.data.roiNames{iROI,1};
                    AUCdata{iROI,3}=strcat('trial', num2str(itrial));
                    AUCdata{iROI,4}= 'GCaMP';
                    AUCdata{iROI,5}= spotId;
                    AUCdata{iROI,6}= CurrentAnimal;
                    AUCdata{iROI,7} = CurrentCondition;
                    AUCdata{iROI,8} = CurrentDepth(1,1);
                end
                All_AUC=vertcat(All_AUC, AUCdata);
                clearvars AUCdata
            end
            clearvars temp temp2
            
            % Neuronal Hand clicked ROIs
            for itrial=1:length(CSArray_Ch2_Hand1)
                
                temp=CSArray_Ch2_Hand1(1,itrial).calcMeasureROIs.data;
                temp2.trialname ={};
                temp2.animalname = {};
                temp2.channel = {};
                temp2.Spot = {};
                temp2.Cond = {};
                temp2.depth = {};
                
                % extract fields from Class
                for jField = 1:numel(listFields)
                    data.(listFields{jField}) = [data.(listFields{jField}); ...
                        temp.(listFields{jField})];
                end
                
                % create fields for trial, animal, spot, condition, etc.
                for iData = 1:length(temp.amplitude)
                    temp2.trialname{iData,1}=strcat('trial', num2str(itrial));
                    temp2.channel{iData,1}= 'RCaMP';
                    temp2.Spot{iData,1}= spotId;
                    temp2.animalname{iData,1}= CurrentAnimal;
                    temp2.Cond{iData,1} = CurrentCondition;
                    temp2.depth{iData,1} = CurrentDepth(1,1);
                end
                
                data.Trial= [data.Trial; temp2.trialname];
                data.Animal= [data.Animal; temp2.animalname];
                data.Channel= [data.Channel; temp2.channel];
                data.Spot= [data.Spot; temp2.Spot];
                data.Condition= [data.Condition; temp2.Cond];
                data.Depth= [data.Depth; temp2.depth];
                
                % trace AUC output
                traces= CSArray_Ch1_FLIKA2(1,itrial).calcMeasureROIs.data.tracesNorm;
                for iROI = 1:size(traces,2)
                    AUCdata{iROI,1}= trapz(traces(60:60+round(10*11.84),iROI));
                    AUCdata{iROI,2}=CSArray_Ch1_FLIKA2(1,itrial).calcFindROIs.data.roiNames{iROI,1};
                    AUCdata{iROI,3}=strcat('trial', num2str(itrial));
                    AUCdata{iROI,4}= 'GCaMP';
                    AUCdata{iROI,5}= spotId;
                    AUCdata{iROI,6}= CurrentAnimal;
                    AUCdata{iROI,7} = CurrentCondition;
                    AUCdata{iROI,8} = CurrentDepth(1,1);
                end
                All_AUC=vertcat(All_AUC, AUCdata);
                clearvars AUCdata
            end
            
            dataNames=fieldnames(data);
            data2= struct2cell(data);
            data3= [data2{:}];
            
            AllData=vertcat(AllData, data3);
            
            if ~exist(fullfile(testRoot,SaveFiles{1,2}),'file')
                %Save CellScan
                save(fullfile(testRoot,SaveFiles{1,2}), 'CSArray_Ch1_FLIKA1');
                save(fullfile(testRoot,SaveFiles{1,3}), 'CSArray_Ch1_Hand1');
                save(fullfile(testRoot,SaveFiles{1,4}), 'CSArray_Ch2_Hand1');
            end
            
            clearvars data data3 temp temp2 
        end
    end
end

% %% Save all data for R analysis
AllData2= [dataNames';AllData];

AUCNames= {'AUC10s','ROI','Trial','Channel','Spot','Animal','Condition','depth'};
All_AUC2=[AUCNames;All_AUC];
%
% cd(fullfile(Settings.MainDir, 'Results'));
% % write date to created file
% cell2csv(SaveFiles{1,1}, AllData2);
% cell2csv(SaveFiles{1,5}, All_AUC2);
