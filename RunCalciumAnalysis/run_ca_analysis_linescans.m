

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
Settings.IP3R2KO = {'IPRG1','IPRG4'};
Settings.IP3R2WT = {'IPRG2','IPRG3'};

channel = struct('Ca_Memb_Astro',1,'Ca_Neuron',2);

% final data file name
SaveFiles{1,1}= fullfile(Settings.MainDir, 'Results','FilesforMatlab','LinescanTraces_allMice_Lck_nostim_vs_longstim_12_2017.mat');
SaveFiles{1,2}= fullfile(Settings.MainDir, 'Results','FilesforR','LinescanOnsets_allMice_Lck_nostim_vs_longstim_12_2017.csv');

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
    
    % Get Drug name
    drugs = unique(Settings.Drug);
    
    for iDrug = 1:length(drugs)
        
        % Extract spot name
        drugId = drugs{iDrug};
        
        % Find the idx of paths matching this spot (all Conditions)
        matchingDrugsIdx = find(~cellfun(@isempty, regexp(Settings.Drug, drugId)));
        SpotPaths2 = Settings.LowresPath(matchingDrugsIdx);
        Depths = Settings.Depth(matchingDrugsIdx);
        Baselines = Settings.BL_frames(matchingDrugsIdx);
        
        % Get SpotIDs
        SpotIDs=Settings.SpotIDs(matchingDrugsIdx);
        spots = unique(Settings.SpotIDs(matchingDrugsIdx));
        
        for iSpot = 1:length(spots) % looks at every 2nd line of SpotIDs-
            
            % Extract spot name
            spotId = spots{iSpot};
            
            % Find the idx of paths matching this spot (all Conditions)
            matchingSpotsIdx = find(~cellfun(@isempty, regexp(SpotIDs, spotId)));
            SpotPaths = SpotPaths2(matchingSpotsIdx);
            CurrentDepth = Depths(matchingSpotsIdx);
            CurrentBaseline = Baselines(matchingSpotsIdx);
            
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
                
                ImgArray(1,1).plot();
                
                
                %% Create CellScan objects
                
                % ASTROCYTES
                load E:\matlab\ca-analysis\Jill_ca_analysis_tests\RunCalciumAnalysis\ConfigCellScanLS1D_LckGC.mat
                CSArray_Ch1_FLIKA = CellScan(fnList, ImgArray, confObj, 1);
                CSArray_Ch1_FLIKA =CSArray_Ch1_FLIKA.process();
                
                CSArray_Ch1_FLIKA.opt_config();
                
                % NEURONS
                load E:\matlab\ca-analysis\Jill_ca_analysis_tests\RunCalciumAnalysis\ConfigCellScanLS1D_RC.mat
                CSArray_Ch2_FLIKA = CellScan(fnList, ImgArray, confObj, 2);
                CSArray_Ch2_FLIKA =CSArray_Ch2_FLIKA.process();
                
                
                %             CSArray_Ch1_FLIKA.plot();
                %             CSArray_Ch2_FLIKA.plot();
                

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
                                Trace_data{iROI,7}= drugId;
                                
                                matchingGenotypeIdx = find(~cellfun(@isempty, regexp(Settings.IP3R2KO, CurrentAnimal)));
                                matchingGenotypeIdx2 = find(~cellfun(@isempty, regexp(Settings.IP3R2WT, CurrentAnimal)));
                                if ~isempty(matchingGenotypeIdx)
                                    Trace_data{iROI,8}= 'IP3R2_KO';
                                elseif ~isempty(matchingGenotypeIdx2)
                                    Trace_data{iROI,8}= 'IP3R2_WT';
                                else
                                    Trace_data{iROI,8}= NaN;
                                end
                                Trace_data{iROI,9} = CurrentDepth(1,1);
                                Trace_data{iROI,10} = CurrentBaseline(1,1);
                                Trace_data{iROI,11} = traces(:,iROI);
                                Trace_data{iROI,12} = CellScans(iScan,itrial).calcFindROIs.data.roiIdxs{iROI,1};
                                Trace_data{iROI,13} = CellScans(iScan,itrial).rawImg.metadata.pixelSize;
                                
                                % line rate- time for one line scan
                                lineTime=(CellScans(1,1).rawImg.metadata.lineTime)/1000;  % time in s
                                Trace_data{iROI,14} = lineTime; %lineRate
                                nLines=length(traces(:,iROI));
                                TimeX(1:nLines) = (1:nLines)*lineTime;
                                
                                % Calculate the first peak onset time and AUC after stim
                                BL_time=round(BL_frames/CellScans(1,1).rawImg.metadata.frameRate);  % number of s for baseline
                                nBL_lines=BL_time/lineTime; % number of lines in baseline time
                                stimStart=TimeX(1,round(nBL_lines));  % exact time of stimulus start
                                baselineCorrectedTime=TimeX-stimStart;
                                
                                % onset time  % 2.5SD from baseline and
                                % smoothing trace at 11 points (5 each side
                                % of middle)
                                Onsets=find_first_onset_time(baselineCorrectedTime(10:end), traces(10:end,iROI),2.5,5);
                                if isempty(Onsets)
                                    Onsets=nan(1,1);
                                end
                                Trace_data{iROI,15}= Onsets;
                                
                            end
                            
                        end
                        All_traces=vertcat(All_traces, Trace_data);
                        clearvars Trace_data
                    end
                end
                
            end
        end
    end
end



% table for importing into R
names={'ROI','Trial','Channel','Spot','Animal', 'Condition','drug','Genotype','depth','baseline',...
    'traces','ROIIdx','PixelSize','lineTime','OnsetTime'};
All_traces2=vertcat(names, All_traces);
All_traces2(:,10)=[];
All_traces2(:,10)=[];

% %% Save traces table
cd(fullfile(Settings.MainDir, 'Results'));
% write date to created file
save(SaveFiles{1,1}, 'All_traces','-v7.3');
cell2csv(SaveFiles{1,2},All_traces2);



