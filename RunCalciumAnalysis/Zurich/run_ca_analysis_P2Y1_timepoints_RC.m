%% 6 cytoGCaMP long stim

clearvars

AllData= [];
All_traces= [];


% NOTE: NO unmixing, NO motion correction

%% Information about your images
Settings.MainDir = 'G:\ZurichData\Astrocyte_Calcium\P2Y1_Mice';

Settings.AnimalNames = {
    'PYRG7',...
    };
Settings.ScoreSheetNames = {
    'PYRG7_Scoresheet_LckGCaMP.xls',...
    };
Settings.NameConditions = {'WPstim','nostim'};

channel = struct('Ca_Memb_Astro',1,'Ca_Neuron',2);
plotMotion =0;

% final data file name
SaveFiles{1,1} = fullfile(Settings.MainDir, 'Results', 'RC_P2Y1_timepoints_peaks_05_2018.csv');%'Control_Peaks_3Conds.csv'); % all data
SaveFiles{1,2} = fullfile(Settings.MainDir, 'Results', 'RC_P2Y1_timepoints_peaks_05_2018.mat');%'Control_Peaks_3Conds.csv'); % all data
SaveFiles{1,3}= fullfile(Settings.MainDir, 'Results','RC_P2Y1_timepoints_traces_05_2018.mat'); %'Control_TraceAUC_20sWindow_3Conds.csv'); % neuronal hand click cell scan

%% Load calibration file
calibration ='C:\JillsFiles\matlab\CalibrationFiles\calibration_20x.mat';
CalFile = CalibrationPixelSize.load(calibration);


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
        
        % Get TimePoints
        CurrentTPs = Settings.TimePoint(matchingSpotsIdx);
        TimePoints = unique(CurrentTPs);
        
        for iTime=1:length(TimePoints)
            % Extract spot name
            TPId = TimePoints{iTime};
            
            % Find the idx of paths matching this spot (all Conditions)
            matchingTimeIdx = find(~cellfun(@isempty, regexp(CurrentTPs, TPId)));
            TimePaths = SpotPaths(matchingTimeIdx);
            
            for iCond = 1:length(Settings.NameConditions)
                if length(TimePaths) < length(Settings.NameConditions)
                    error('not all stimulus conditions are present')
                end
                
                CurrentCondition = Settings.NameConditions{iCond};
                PathIdx = find(~cellfun(@isempty, regexp(TimePaths, CurrentCondition)));
                BL_frames = CurrentBaseline(PathIdx);
                
                % Get image paths
                testRoot =TimePaths{PathIdx};
                
                
                expfiles = dir(fullfile(testRoot,'lowres*'));
                fnTempList = {expfiles(:).name};
                fnList = fullfile(testRoot, fnTempList);
                
                % Create an array of ScanImage Tiffs
                ImgArray =  SCIM_Tif(fnList, channel, CalFile);
                
                %% Configs for Finding ROIs
                % Neurons
                % hand selected- peaks from cellular structures
                x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
                AC_findConf{2} = ConfigFindROIsDummy.from_ImageJ(fullfile(testRoot,'Neurons.zip'), x_pix, y_pix,1);
                
                
                
                %% Configuration for measuring ROIs
                % neuron peaks
                detectConf{1} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...'normMethod', 'z-score','zIters', 100,...
                    'propagateNaNs', false, 'excludeNaNs', false, 'lpWindowTime', 2, 'spFilterOrder', 6,...
                    'spPassBandMin',0.05, 'spPassBandMax', 0.5, 'thresholdLP', 5,'thresholdSP', 7);
                % for calculating AUC for each trace
                measureConf = ConfigMeasureROIsDummy('baselineFrames', BL_frames);
                
                %%
                % Combine the configs into a CellScan config
                configCS{1,2} = ConfigCellScan(AC_findConf{1,2}, measureConf,detectConf{1,1}); % neuron hand, peaks
                
                %% Create CellScan objects
                
                CSArray_Ch2_Hand = CellScan(fnList, ImgArray, configCS{1,2}, 2);
                
                
                
                %% Process the images
                CSArray_Ch2_Hand =CSArray_Ch2_Hand.process();
                
                
                %% Make the debugging plots
                % CSArray_Ch2_FLIKA.plot;
                %                     CSArray_Ch2_Hand.plot;
                %                     CSArray_Ch1_FLIKA.plot;
                %CSArray_Ch2_Hand.opt_config();
                % for working out find peaks parameters
                %CSArray_Ch2_FLIKA(1).opt_config()
                
                % CSArray_Ch1_FLIKA(2).opt_config()
                %CSArray_Ch1_FLIKA(1).opt_config()
                %CSArray_Ch1_FLIKA(4).opt_config()
                %% Output data
                
                % make a giant data table
                listFields = {'amplitude', 'fullWidth', 'halfWidth', ...
                    'numPeaks', 'peakTime', 'peakStart', 'peakStartHalf', ...
                    'peakType', 'prominence', 'roiName', 'peakAUC'};
                
                % Neuron Hand clicked ROIs
                for itrial=1:length(CSArray_Ch2_Hand)
                    
                    temp=CSArray_Ch2_Hand(1,itrial).calcDetectSigs.data;
                    temp2.trialname ={};
                    temp2.animalname = {};
                    temp2.channel = {};
                    temp2.Spot = {};
                    temp2.Cond = {};
                    temp2.depth = {};
                    temp2.overlap = {};
                    temp2.area = {};
                    temp2.pixelsize = {};
                    temp2.TimePoint= {};
                    
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
                    for iPeak = 1:length(temp.amplitude)
                        temp2.trialname{iPeak,1}=strcat('trial', num2str(itrial,'%02d'));
                        temp2.channel{iPeak,1}= 'RCaMP';
                        temp2.Spot{iPeak,1}= spotId;
                        temp2.animalname{iPeak,1}= CurrentAnimal;
                        temp2.Cond{iPeak,1} = CurrentCondition;
                        temp2.depth{iPeak,1} = CurrentDepth(1,1);
                        temp2.area{iPeak,1} = 0;
                        temp2.overlap{iPeak,1} = 0;
                        temp2.pixelsize{iPeak,1} = CSArray_Ch2_Hand(1,itrial).rawImg.metadata.pixelSize;
                        temp2.TimePoint{iPeak,1}= TPId;
                    end
                    
                    
                    isFirst = (itrial == 1 );
                    if isFirst
                        data.Trial = {};
                        data.Animal = {};
                        data.Channel = {};
                        data.Spot = {};
                        data.Condition = {};
                        data.Depth = {};
                        data.area = {};
                        data.overlap = {};
                        data.pixelsize={};
                         data.TimePoint={};
                    end
                    data.Trial= [data.Trial; temp2.trialname];
                    data.Animal= [data.Animal; temp2.animalname];
                    data.Channel= [data.Channel; temp2.channel];
                    data.Spot= [data.Spot; temp2.Spot];
                    data.Condition= [data.Condition; temp2.Cond];
                    data.Depth= [data.Depth; temp2.depth];
                    data.area= [data.area; temp2.area];
                    data.overlap= [data.overlap; temp2.overlap];
                    data.pixelsize= [data.pixelsize; temp2.pixelsize];
                    data.TimePoint= [data.TimePoint; temp2.TimePoint];
                    
                    %traces output
                    traces= CSArray_Ch2_Hand(1,itrial).calcMeasureROIs.data.tracesNorm;
                    %preallocate
                    Trace_data=cell(size(traces,2),10);
                    for iROI = 1:size(traces,2)
                        Trace_data{iROI,1}= CSArray_Ch2_Hand(1,itrial).calcFindROIs.data.roiNames{iROI,1};
                        Trace_data{iROI,2}= strcat('trial', num2str(itrial,'%02d'));
                        Trace_data{iROI,3}= 'RCaMP';
                        Trace_data{iROI,4}= spotId;
                        Trace_data{iROI,5}= CurrentAnimal;
                        Trace_data{iROI,6}= CurrentCondition;
                        Trace_data{iROI,7} = CurrentDepth(1,1);
                        Trace_data{iROI,8} = traces(:,iROI);
                        Trace_data{iROI,9} = CSArray_Ch2_Hand(1,itrial).calcFindROIs.data.centroidX(iROI,1);
                        Trace_data{iROI,10} = CSArray_Ch2_Hand(1,itrial).calcFindROIs.data.centroidY(iROI,1);
                        Trace_data{iROI,11} = CSArray_Ch2_Hand(1,itrial).calcFindROIs.data.roiMask(:,:,iROI);
                        Trace_data{iROI,14} = 0;
                        Trace_data{iROI,12} = CSArray_Ch2_Hand(1,itrial).rawImg.metadata.pixelSize;
                        Trace_data{iROI,13}= TPId;
                    end
                    All_traces=vertcat(All_traces, Trace_data);
                    clearvars Trace_data
                end
                
                
                dataNames=fieldnames(data);
                data2= struct2cell(data);
                data3= [data2{:}];
                
                AllData=vertcat(AllData, data3);
                
                
                clearvars data data3 temp temp2
                
            end
        end
    end
end

% %% Save all data for R analysis
AllData2= [dataNames';AllData];

cd(fullfile(Settings.MainDir, 'Results'));
% write date to created file
cell2csv(SaveFiles{1,1}, AllData2);
save(SaveFiles{1,3}, 'All_traces','-v7.3');
save(SaveFiles{1,2}, 'AllData2','-v7.3');


