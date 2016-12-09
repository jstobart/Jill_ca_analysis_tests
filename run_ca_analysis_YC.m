%% Clear workspace
clearvars

AllData= [];

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

channel = struct('Ca_Neuron',2);

plotMotion = 0; %Plot motion correction movie
doplots = 1; %Plots for each trial

% final data file name
SaveFiles{1,1} = fullfile(Settings.MainDir, 'Results', 'Tests.csv'); % all data
SaveFiles{1,2}= 'CellScan_Neurons.mat';

% Load calibration file
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
    
    for iSpot = 1:length(spots)         
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
                
                ImgArray = ImgArray.motion_correct('refImg',refImg, 'ch',2,'minCorr', 0.4);%,'doPlot',true);
                % fill in bad data? 
                if plotMotion==1
                    ImgArray(:).plot()
                end
                
                %% Configs for Finding ROIs 
                % Neuronal Calcium
                % hand selected
                x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
                findConf = ConfigFindROIsDummy.from_ImageJ(fullfile(testRoot,'RoiSet.zip'), x_pix, y_pix,1);
                
                %% Configuration for measuring ROIs               
                % AWAKE neuron calcium
                measureConf = ConfigMeasureROIsClsfy('baselineFrames', BL_frames,... 'normMethod','z-score', 'zIters', 100
                    'propagateNaNs', false,'excludeNaNs', false, 'lpWindowTime', 3, 'spFilterOrder', 2,...
                   'spPassBand', [0.1, 2], 'thresholdSD_low', 5,'thresholdSD_band', 7);
          
                %%
                % Combine the configs into a CellScan config
                configCS = ConfigCellScan(findConf, measureConf);               
                
                %% Create CellScan objects
                CSArray_Ch2 = CellScan(fnList, ImgArray, configCS, 2);
                
                %% Process the images
                CSArray_Ch2 =CSArray_Ch2.process();

                %% Make the debugging plots
                if doplots ==1
                    CSArray_Ch2.plot;
                    % CSArray_Ch2.plot('peaks');
                end
            else
                load (fullfile(testRoot,SaveFiles{1,2}));
            end
            
            %% Output data
            
            % make a giant data table
            
            listFields = {'amplitude', 'area', 'fullWidth', 'halfWidth', ...
                'numPeaks', 'peakTime', 'peakStart', 'peakStartHalf', ...
                'peakType', 'prominence', 'ROIname', 'peakAUC'};
            
            % Astrocyte FLIKA
            for itrial=1:length(CSArray_Ch2)
                
                temp=CSArray_Ch2(1,itrial).calcMeasureROIs.data;
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
                    temp2.channel{iData,1}= 'RCaMP';
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
                
            end
            clearvars temp temp2
            
            dataNames=fieldnames(data);           
            data2= struct2cell(data);
            data3= [data2{:}];
            
            AllData=vertcat(AllData, data3);
            
            if ~exist(fullfile(testRoot,SaveFiles{1,2}),'file')
                %Save CellScan
                save(fullfile(testRoot,SaveFiles{1,2}), 'CSArray_Ch2');
            end
            
            clearvars data data3 temp temp2
        end
    end
end

%% Save all data for R analysis
AllData2= [dataNames';AllData];

cd(fullfile(Settings.MainDir, 'Results'));
% write date to created file
cell2csv(SaveFiles{1,1}, AllData2);