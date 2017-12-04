

% NOTE:::::: RUN FROM THE 1D branch!

% finds ROIs and outputs traces from linescans
%
%% 1 lck no stim vs long stim for all mice (independent of genotype)
clearvars

All_traces= [];


%% Information about your images
Settings.MainDir = 'E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f';

Settings.AnimalNames = {
    'ARG2',...
    'IPRG1',...
    'IPRG2',...
    'IPRG3',...
    'IPRG4',...
    };
Settings.ScoreSheetNames = {
    'ARG2_Scoresheet_LineScans.xls',...
    'IPRG1_Scoresheet_LineScans.xls',...
    'IPRG2_Scoresheet_LineScans.xls',...
    'IPRG3_Scoresheet_LineScans.xls',...
    'IPRG4_Scoresheet_LineScans.xls',...
    };
Settings.NameConditions = {'Nostim','Stim'};

channel = struct('Ca_Memb_Astro',1,'Ca_Neuron',2);

% final data file name
SaveFiles{1,1}= fullfile(Settings.MainDir, 'Results','Traces_allMice_Lck_nostim_vs_longstim_12_2017.mat'); %'Control_TraceAUC_20sWindow_3Conds.csv'); % neuronal hand click cell scan

%% Load calibration file
calibration ='E:\matlab\CalibrationFiles\calibration_20x.mat';
CalFile = CalibrationPixelSize.load(calibration);


%% load scoresheet and loop through animal, spot, etc.

Settings.ScoreSheetPath = fullfile(Settings.MainDir,'Scoresheets',Settings.ScoreSheetNames);

numAnimals = length(Settings.AnimalNames);
for iAnimal = 1:numAnimals
    CurrentAnimal = Settings.AnimalNames{iAnimal};
    %read Scoresheet of current animal
    CurrentSheet = Settings.ScoreSheetPath{iAnimal};
    Settings = readScoresheet2(CurrentSheet, Settings);
    
    % load spectral unmixing matrix of RCaMP and GCaMP
    load(fullfile(Settings.MainDir, 'Results','RCaMP_mGCaMP_Matrix_4Ch.mat'));
    
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
        
        for iCond = 1:length(Settings.NameConditions)
            if length(SpotPaths) < length(Settings.NameConditions)
                error('not all stimulus conditions are present')
            end
            
            CurrentCondition = Settings.NameConditions{iCond};
            PathIdx = find(~cellfun(@isempty, regexp(SpotPaths, CurrentCondition)));
            BL_frames = CurrentBaseline(PathIdx);
            
            % Get image paths
            testRoot =SpotPaths{PathIdx};
            
            
            expfiles = dir(fullfile(testRoot,'lowres*'));
            fnTempList = {expfiles(:).name};
            fnList = fullfile(testRoot, fnTempList);
            
            % Create an array of ScanImage Tiffs
            ImgArray =  SCIM_Tif(fnList, channel, CalFile);
            
            
            % Spectral Unmixing of GCaMP and RCaMP
            ImgArray= ImgArray.unmix_chs(false, [], cell2mat(RCaMP_mGCaMP_Matrix));
            
            %ImgArray(1,1).plot();
            
            
            %% Create CellScan objects
            
            % ASTROCYTES
            load E:\matlab\ca-analysis\Jill_ca_analysis_tests\RunCalciumAnalysis\ConfigCellScanLS1D_LckGC.mat
            CSArray_Ch1_FLIKA = CellScan(fnList, ImgArray, confObj, 1);
            
            % NEURONS
            load E:\matlab\ca-analysis\Jill_ca_analysis_tests\RunCalciumAnalysis\ConfigCellScanLS1D_RC.mat
            CSArray_Ch2_FLIKA = CellScan(fnList, ImgArray, confObj, 2);
            
            
            %% Process the images
            CSArray_Ch1_FLIKA =CSArray_Ch1_FLIKA.process();
            CSArray_Ch2_FLIKA =CSArray_Ch2_FLIKA.process();
            
            %             CSArray_Ch1_FLIKA.plot();
            %             CSArray_Ch2_FLIKA.plot();
            
            % CSArray_Ch1_FLIKA.opt_config();
            %% Calculate onset time for first peak and Output data
            
            CellScans=vertcat(CSArray_Ch1_FLIKA, CSArray_Ch2_FLIKA);
            
            
            % loop through cellscans
            for iScan=1:size(CellScans,1)
                for itrial=1:size(CellScans,2)
                    
                    % make a table of trace info
                    if strcmp(CellScans(iScan,itrial).calcFindROIs.data.roiNames{1,1}, 'none')
                        continue
                    else
                        traces= CellScans(iScan,itrial).calcMeasureROIs.data.tracesNorm;
                        %preallocate
                        Trace_data=cell(size(traces,2),10);
                        for iROI = 1:size(traces,2)
                            Trace_data{iROI,1}= CellScans(iScan,itrial).calcFindROIs.data.roiNames{iROI,1};
                            Trace_data{iROI,2}= strcat('trial', num2str(itrial,'%02d'));
                            if iScan==1
                                Trace_data{iROI,3}= 'GCaMP';
                            else
                                Trace_data{iROI,3}= 'RCaMP';
                            end
                            
                            Trace_data{iROI,4}= spotId;
                            Trace_data{iROI,5}= CurrentAnimal;
                            Trace_data{iROI,6}= CurrentCondition;
                            Trace_data{iROI,7} = CurrentDepth(1,1);
                            Trace_data{iROI,8} = CurrentBaseline(1,1);
                            Trace_data{iROI,9} = traces(:,iROI);
                            Trace_data{iROI,10} = CellScans(iScan,itrial).calcFindROIs.data.roiIdxs{iROI,1};
                            Trace_data{iROI,11} = CellScans(iScan,itrial).rawImg.metadata.pixelSize;
                            
                            % line rate- time for one line scan
                            lineTime=CellScans(1,1).rawImg.metadata.lineTime;
                            Trace_data{iROI,12} = lineTime; %lineRate
                            nLines=length(traces(:,iROI));
                            TimeX(1:nLines) = (1:nLines)*lineTime;
                            
                            % Calculate the first peak onset time and AUC after stim
                            
                            baselineCorrectedTime=TimeX-BL_frames;
                            
                            % onset time
                            Onsets=find_first_onset_time(baselineCorrectedTime(10:end), trace(10:end),2.5,5);
                            if isempty(Onsets)
                                Onsets=nan(1,1);
                            end
                            Trace_data{iROI,13}= Onsets;                           
                            
                        end
                        
                    end
                    All_traces=vertcat(All_traces, Trace_data);
                    clearvars Trace_data
                end
            end
            
        end
    end
end



% table for importing into R
    names={'ROI','Trial','Channel','Spot','Animal', 'Condition','depth','baseline',...
        'traces','ROIIdx','PixelSize','lineTime','OnsetTime'};
    All_traces=vertcat(names, All_traces);

% %% Save traces table
cd(fullfile(Settings.MainDir, 'Results'));
% write date to created file
save(SaveFiles{1,1}, 'All_traces','-v7.3');



