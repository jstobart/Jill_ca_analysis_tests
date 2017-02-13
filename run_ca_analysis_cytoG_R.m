%% Clear workspace
clearvars

AllData= [];
All_traces= [];

%% Information about your images
Settings.MainDir = 'E:\Data\Two_Photon_Data\GCaMP_RCaMP\cyto_GCaMP6s';

Settings.AnimalNames = {
    'RG12',...
    'RG14',...
    'RG16',...
    'RG17',...
    'RG18',...
    };
Settings.ScoreSheetNames = {
    'RG12_Scoresheet_LongTrials.xls',...
    'RG14_Scoresheet_LongTrials.xls',...
    'RG16_Scoresheet_LongTrials.xls',...
    'RG17_Scoresheet_LongTrials.xls',...
    'RG18_Scoresheet_LongTrials.xls',...
    };
Settings.NameConditions = {'Nostim','Stim'};
%Settings.NameConditions = {'Stim'};

channel = struct('Ca_Cyto_Astro',1,'Ca_Neuron',2);

plotMotion = 0; %Plot motion correction movie
doplots = 0; %Plots for each trial

% final data file name
SaveFiles{1,1} = fullfile(Settings.MainDir, 'Results', 'LStim_cGC&RC_11_02_2017.csv');%'Control_Peaks_3Conds.csv'); % all data
SaveFiles{1,2} = fullfile(Settings.MainDir, 'Results', 'LStim_cGC&RC_11_02_2017.mat');%'Control_Peaks_3Conds.csv'); % all data
SaveFiles{1,3}= fullfile(Settings.MainDir, 'Results','LStim_cGC&RC_traces_11_02_2017.mat'); %'Control_TraceAUC_20sWindow_3Conds.csv'); % neuronal hand click cell scan

%% Load calibration file
calibration ='E:\matlab\2p-img-analysis\tests\res\calibration_20x.mat';
CalFile = Calibration_PixelSize.load(calibration);

%% make a mixing matrix for spectral unmixing of RCaMP and GCaMP

if ~exist(fullfile(Settings.MainDir, 'Results','RCaMP_cGCaMP_Matrix.mat'),'file')
    % load example high res image
    unmixFile = 'E:\Data\Two_Photon_Data\GCaMP_RCaMP\cyto_GCaMP6s\RG12\Awake\2016_02_04\spot1_long_Nostim\highres_spot1_long096.tif';
    unmixImg = SCIM_Tif(unmixFile, channel, CalFile);
    
    [unmixImg, RCaMP_cGCaMP_Matrix] = unmix_chs(unmixImg); %returns the mixing matrix to be used for all imaging
    
    % save the mixing matrix for loading later
    cd(fullfile(Settings.MainDir, 'Results'));
    % write matrix to created file
    save('RCaMP_cGCaMP_Matrix.mat', 'RCaMP_cGCaMP_Matrix');
else
    load(fullfile(Settings.MainDir, 'Results','RCaMP_cGCaMP_Matrix.mat'));
end


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
            
            expfiles = dir(fullfile(testRoot,'lowres*'));
            fnTempList = {expfiles(:).name};
            fnList = fullfile(testRoot, fnTempList);
            
            % Create an array of ScanImage Tiffs
            ImgArray =  SCIM_Tif(fnList, channel, CalFile);
            
            %exclude frames at the beginning of each trial where pockel
            %cell was dark
            ImgArray.exclude_frames('badframes',(1:2),'method', 'inpaint','inpaintIters',5);
            
            
            % Spectral Unmixing of GCaMP and RCaMP
            ImgArray= ImgArray.unmix_chs(false, [], RCaMP_cGCaMP_Matrix);
            
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
                ImgArray(1,3).plot('video');
            end
            
            %% Configs for Finding ROIs
            %Astrocyte Calcium
            % automated selection
            findConf{1} = ConfigFindROIsFLIKA_2D.from_preset('ca_cyto_astro', 'baselineFrames',...
                BL_frames,'sigmaXY', 2,...
                'sigmaT', 0.5,'threshold_std', 5,...
                'min_rise_time',0.1689, 'max_rise_time', 8,...
                'dilateXY', 4,...
                'dilateT', 0.5,'erodeXY', 2);
            
            
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
            measureConf{1} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...'normMethod', 'z-score','zIters', 100,...
                'propagateNaNs', false, 'excludeNaNs', false, 'lpWindowTime', 3, 'spFilterOrder', 2,...
                'spPassBandMin',0.05, 'spPassBandMax', 0.2, 'thresholdSD_low', 3,'thresholdSD_band', 5);
            
            % AWAKE neuron calcium
            measureConf{2} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,... 'normMethod','z-score',... 'zIters', 10000,...
                'propagateNaNs', false,'excludeNaNs', false, 'lpWindowTime', 2, 'spFilterOrder', 2,...
                'spPassBandMin',0.1, 'spPassBandMax', 1, 'thresholdSD_low', 3,'thresholdSD_band', 5);
            
            % for calculating AUC for each trace
            detectConf = ConfigMeasureROIsDummy('baselineFrames', BL_frames);
            
            %%
            % Combine the configs into a CellScan config
            configCS{1,1} = ConfigCellScan(findConf{1,1}, measureConf{1,1},detectConf); % astrocyte FLIKA, peaks
            configCS{1,2} = ConfigCellScan(findConf{1,2}, measureConf{1,1},detectConf); % astrocyte hand, peaks
            configCS{1,3} = ConfigCellScan(findConf{1,3}, measureConf{1,2},detectConf); % neuronal hand, peaks
            
            %% Create CellScan objects
            CSArray_Ch1_FLIKA = CellScan(fnList, ImgArray, configCS{1,1}, 1);
            
            CSArray_Ch1_Hand = CellScan(fnList, ImgArray, configCS{1,2}, 1);
            
            CSArray_Ch2_Hand = CellScan(fnList, ImgArray, configCS{1,3}, 2);
            
            %% Process the images
            CSArray_Ch1_FLIKA =CSArray_Ch1_FLIKA.process();
            CSArray_Ch1_Hand =CSArray_Ch1_Hand.process();
            CSArray_Ch2_Hand =CSArray_Ch2_Hand.process();
            
            %NO IMAGE COMBINATION
            % if isfield(channel,'Ca_Cyto_Astro')
            % Use this function to combine the masks
            %CSArray_Ch1_combined= combine_masks(Ch1_test);
            %CSArray_Ch1_combined.process();
            %end
            
            %% Make the debugging plots
            if doplots ==1
                CSArray_Ch1_FLIKA.plot;
                CSArray_Ch1_Hand.plot;
                CSArray_Ch2_Hand.plot;
                
                % for working out find peaks parameters
                %CSArray_Ch1_FLIKA(1).plot('parameters','plotROIs',6:15);
                %CSArray_Ch2_Hand(1).plot('parameters','plotROIs',11:15);
                
                % for viewing peaks that were found
                %             CSArray_Ch1_FLIKA.plot('signals');
                %             CSArray_Ch1_Hand.plot('signals');
                %            CSArray_Ch2_Hand.plot('signals');
            end
            
            %% Output data
            
            % make a giant data table
            listFields = {'amplitude', 'fullWidth', 'halfWidth', ...
                'numPeaks', 'peakTime', 'peakStart', 'peakStartHalf', ...
                'peakType', 'prominence', 'ROIname', 'peakAUC'};
            
            % Astrocyte FLIKA
            for itrial=1:length(CSArray_Ch1_FLIKA)
                
                % peak output
                temp=CSArray_Ch1_FLIKA(1,itrial).calcDetectSigs.data;
                temp2.trialname ={};
                temp2.animalname = {};
                temp2.channel = {};
                temp2.Spot = {};
                temp2.Cond = {};
                temp2.depth = {};
                temp2.overlap = {};
                temp2.area = {};
                
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
                    temp2.trialname{iPeak,1}=strcat('trial', num2str(itrial));
                    temp2.channel{iPeak,1}= 'GCaMP';
                    temp2.Spot{iPeak,1}= spotId;
                    temp2.animalname{iPeak,1}= CurrentAnimal;
                    temp2.Cond{iPeak,1} = CurrentCondition;
                    temp2.depth{iPeak,1} = CurrentDepth(1,1);
                    
                    % get the indices  and area for a particular ROI
                    jROIname = temp.ROIname{iPeak};
                    ROIindex= strcmp(CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.roiNames,jROIname);
                    temp2.area{iPeak,1} = CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.area(ROIindex);
                    
                    ROIpuffIdx = CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.puffIdxs(ROIindex);
                    x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
                    [jy,jx,~] = ind2sub([y_pix,x_pix],cell2mat(ROIpuffIdx));
                    jx = round(mean(jx));
                    jy = round(mean(jy));
                    xclose = round(x_pix*0.03); % 3 percent of pixels
                    yclose = round(y_pix*0.03); % 3 percent of pixels
                    
                    % find astrocyte process and soma ROIs that have
                    % similar centroids
                    for kROI= 1:length(CSArray_Ch1_Hand(1,itrial).calcFindROIs.data.roiNames)
                        % get the indices of the handclicked ROIs
                        ROIMask{kROI}= CSArray_Ch1_Hand(1, 1).calcFindROIs.data.roiMask(:,:,kROI);
                        ROI_Idx{kROI} = find(ROIMask{kROI});
                        % mean pixels for soma ROI
                        [ky,kx,~] = ind2sub([y_pix,x_pix],ROI_Idx{kROI});
                        kx = round(mean(kx));
                        ky = round(mean(ky));
                        
                        % check if they are too close (actually same region)
                        if jx<=(kx+xclose) && jx>=(kx-xclose) && jy<=(ky+yclose) && jy>=(ky-yclose) % only %3 of pixels apart
                            spatialcorr(kROI)  = 1;
                        end
                    end
                    
                    if exist('spatialcorr','var')
                        indx = find(spatialcorr>0);
                        temp2.overlap{iPeak,1}= CSArray_Ch1_Hand(1,itrial).calcFindROIs.data.roiNames{indx};
                    else
                        temp2.overlap{iPeak,1} = 0;
                    end
                    clear spatialcorr
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
                end
                data.Trial= [data.Trial; temp2.trialname];
                data.Animal= [data.Animal; temp2.animalname];
                data.Channel= [data.Channel; temp2.channel];
                data.Spot= [data.Spot; temp2.Spot];
                data.Condition= [data.Condition; temp2.Cond];
                data.Depth= [data.Depth; temp2.depth];
                data.area= [data.area; temp2.area];
                data.overlap= [data.overlap; temp2.overlap];
                
                %traces output processes
                if strcmp(CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.roiNames{1,1}, 'none')
                    continue
                else
                    traces= CSArray_Ch1_FLIKA(1,itrial).calcMeasureROIs.data.tracesNorm;
                    %preallocate
                    Trace_data=cell(size(traces,2),10);
                    for iROI = 1:size(traces,2)
                        Trace_data{iROI,1}= CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.roiNames{iROI,1};
                        Trace_data{iROI,2}= strcat('trial', num2str(itrial));
                        Trace_data{iROI,3}= 'GCaMP';
                        Trace_data{iROI,4}= spotId;
                        Trace_data{iROI,5}= CurrentAnimal;
                        Trace_data{iROI,6}= CurrentCondition;
                        Trace_data{iROI,7} = CurrentDepth(1,1);
                        Trace_data{iROI,8} = traces(:,iROI);
                        Trace_data{iROI,9} = CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.centroid{iROI,1};
                        Trace_data{iROI,10} = CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.puffIdxs{iROI,1};
                        
                        % get the indices  and area for a particular ROI
                        ROIpuffIdx = CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.puffIdxs{iROI,1};
                        x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
                        [jy,jx,~] = ind2sub([y_pix,x_pix],ROIpuffIdx);
                        jx = round(mean(jx));
                        jy = round(mean(jy));
                        xclose = round(x_pix*0.03); % 3 percent of pixels
                        yclose = round(y_pix*0.03); % 3 percent of pixels
                        
                        % find astrocyte process and soma ROIs that have
                        % similar centroids
                        for kROI= 1:length(CSArray_Ch1_Hand(1,itrial).calcFindROIs.data.roiNames)
                            % get the indices of the handclicked ROIs
                            ROIMask{kROI}= CSArray_Ch1_Hand(1, 1).calcFindROIs.data.roiMask(:,:,kROI);
                            ROI_Idx{kROI} = find(ROIMask{kROI});
                            % mean pixels for soma ROI
                            [ky,kx,~] = ind2sub([y_pix,x_pix],ROI_Idx{kROI});
                            kx = round(mean(kx));
                            ky = round(mean(ky));
                            
                            % check if they are too close (actually same region)
                            if jx<=(kx+xclose) && jx>=(kx-xclose) && jy<=(ky+yclose) && jy>=(ky-yclose) % only %3 of pixels apart
                                spatialcorr(kROI)  = 1;
                            end
                        end
                        
                        if exist('spatialcorr','var')
                            indx = find(spatialcorr>0);
                            Trace_data{iROI,11}= CSArray_Ch1_Hand(1,itrial).calcFindROIs.data.roiNames{indx};
                        else
                            Trace_data{iROI,11} = 0;
                        end
                        clear spatialcorr
                        
                        
                    end
                    All_traces=vertcat(All_traces, Trace_data);
                    clearvars Trace_data
                end
            end
            clearvars temp temp2
            
            % Astrocyte Hand clicked ROIs
            for itrial=1:length(CSArray_Ch1_Hand)
                
                temp=CSArray_Ch1_Hand(1,itrial).calcDetectSigs.data;
                temp2.trialname ={};
                temp2.animalname = {};
                temp2.channel = {};
                temp2.Spot = {};
                temp2.Cond = {};
                temp2.depth = {};
                temp2.overlap = {};
                temp2.area = {};
                
                % extract fields from Class
                for jField = 1:numel(listFields)
                    data.(listFields{jField}) = [data.(listFields{jField}); ...
                        temp.(listFields{jField})];
                end
                
                % create fields for trial, animal, spot, condition, etc.
                for iPeak = 1:length(temp.amplitude)
                    temp2.trialname{iPeak,1}=strcat('trial', num2str(itrial));
                    temp2.channel{iPeak,1}= 'GCaMP';
                    temp2.Spot{iPeak,1}= spotId;
                    temp2.animalname{iPeak,1}= CurrentAnimal;
                    temp2.Cond{iPeak,1} = CurrentCondition;
                    temp2.depth{iPeak,1} = CurrentDepth(1,1);
                    temp2.area{iPeak,1} = 0;
                    temp2.overlap{iPeak,1} = 0;
                end
                
                data.Trial= [data.Trial; temp2.trialname];
                data.Animal= [data.Animal; temp2.animalname];
                data.Channel= [data.Channel; temp2.channel];
                data.Spot= [data.Spot; temp2.Spot];
                data.Condition= [data.Condition; temp2.Cond];
                data.Depth= [data.Depth; temp2.depth];
                data.area= [data.area; temp2.area];
                data.overlap= [data.overlap; temp2.overlap];
                
                %traces output
                traces= CSArray_Ch1_Hand(1,itrial).calcMeasureROIs.data.tracesNorm;
                %preallocate
                Trace_data=cell(size(traces,2),10);
                for iROI = 1:size(traces,2)
                    Trace_data{iROI,1}= CSArray_Ch1_Hand(1,itrial).calcFindROIs.data.roiNames{iROI,1};
                    Trace_data{iROI,2}= strcat('trial', num2str(itrial));
                    Trace_data{iROI,3}= 'GCaMP';
                    Trace_data{iROI,4}= spotId;
                    Trace_data{iROI,5}= CurrentAnimal;
                    Trace_data{iROI,6}= CurrentCondition;
                    Trace_data{iROI,7} = CurrentDepth(1,1);
                    Trace_data{iROI,8} = traces(:,iROI);
                    Trace_data{iROI,9} = 0;
                    Trace_data{iROI,10} = CSArray_Ch1_Hand(1,itrial).calcFindROIs.data.roiMask(:,:,iROI);
                    Trace_data{iROI,11} = 0;
                end
                All_traces=vertcat(All_traces, Trace_data);
                clearvars Trace_data
            end
            clearvars temp temp2
            
            % Neuronal Hand clicked ROIs
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
                
                % extract fields from Class
                for jField = 1:numel(listFields)
                    data.(listFields{jField}) = [data.(listFields{jField}); ...
                        temp.(listFields{jField})];
                end
                
                % create fields for trial, animal, spot, condition, etc.
                for iPeak = 1:length(temp.amplitude)
                    temp2.trialname{iPeak,1}=strcat('trial', num2str(itrial));
                    temp2.channel{iPeak,1}= 'RCaMP';
                    temp2.Spot{iPeak,1}= spotId;
                    temp2.animalname{iPeak,1}= CurrentAnimal;
                    temp2.Cond{iPeak,1} = CurrentCondition;
                    temp2.depth{iPeak,1} = CurrentDepth(1,1);
                    temp2.area{iPeak,1} = 0;
                    temp2.overlap{iPeak,1} = 0;
                end
                
                data.Trial= [data.Trial; temp2.trialname];
                data.Animal= [data.Animal; temp2.animalname];
                data.Channel= [data.Channel; temp2.channel];
                data.Spot= [data.Spot; temp2.Spot];
                data.Condition= [data.Condition; temp2.Cond];
                data.Depth= [data.Depth; temp2.depth];
                data.area= [data.area; temp2.area];
                data.overlap= [data.overlap; temp2.overlap];
                
                %traces output
                traces= CSArray_Ch2_Hand(1,itrial).calcMeasureROIs.data.tracesNorm;
                %preallocate
                Trace_data=cell(size(traces,2),10);
                for iROI = 1:size(traces,2)
                    Trace_data{iROI,1}= CSArray_Ch2_Hand(1,itrial).calcFindROIs.data.roiNames{iROI,1};
                    Trace_data{iROI,2}= strcat('trial', num2str(itrial));
                    Trace_data{iROI,3}= 'RCaMP';
                    Trace_data{iROI,4}= spotId;
                    Trace_data{iROI,5}= CurrentAnimal;
                    Trace_data{iROI,6}= CurrentCondition;
                    Trace_data{iROI,7} = CurrentDepth(1,1);
                    Trace_data{iROI,8} = traces(:,iROI);
                    Trace_data{iROI,9} = 0;
                    Trace_data{iROI,10} = CSArray_Ch2_Hand(1,itrial).calcFindROIs.data.roiMask(:,:,iROI);
                    Trace_data{iROI,11} = 0;
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

% %% Save all data for R analysis
AllData2= [dataNames';AllData];

%TraceNames= {'ROI','Trial','Channel','Spot','Animal','Condition','depth','traces','centroid','puffIdx'};
%All_traces2=[TraceNames;All_traces];
%
cd(fullfile(Settings.MainDir, 'Results'));
% write date to created file
cell2csv(SaveFiles{1,1}, AllData2);
save(SaveFiles{1,3}, 'All_traces','-v7.3');
save(SaveFiles{1,2}, 'AllData2','-v7.3');


% %% DSP4
% %% Clear workspace
% clearvars
%
% AllData= [];
% All_AUC= [];
%
% %% Information about your images
% Settings.MainDir = 'E:\Data\Two_Photon_Data\GCaMP_RCaMP\cyto_GCaMP6s';
%
% Settings.AnimalNames = {
%     'RG10',...
%     'RG12',...
%     'RG14',...
%     'RG18',...
%     };
% Settings.ScoreSheetNames = {
%     'RG10_Scoresheet_LongTrialsShortStim_DSP4.xls',...
%     'RG12_Scoresheet_LongTrialsShortStim_DSP4.xls',...
%     'RG14_Scoresheet_LongTrialsShortStim_DSP4.xls',...
%     'RG18_Scoresheet_LongTrialsShortStim_DSP4.xls',...
%     };
% Settings.NameConditions = {'Nostim','Stim','shortstim'};
%
% channel = struct('Ca_Cyto_Astro',1,'Ca_Neuron',2);
%
% plotMotion = 0; %Plot motion correction movie
% doplots = 0; %Plots for each trial
%
% % final data file name
% SaveFiles{1,1} = fullfile(Settings.MainDir, 'Results', 'DSP4_cGC&RC_01_2017.csv');%'Control_Peaks_3Conds.csv'); % all data
% SaveFiles{1,2}= 'CellScan_AC_FLIKA.mat'; % astrocyte FLIKA cell scan
% SaveFiles{1,3}= 'CellScan_AC_Hand.mat'; % astrocyte hand click cell scan
% SaveFiles{1,4}= 'CellScan_Ne_Hand.mat'; % neuronal hand click cell scan
% SaveFiles{1,5}= fullfile(Settings.MainDir, 'Results','DSP4_cGC&RC_AUC_01_2017.csv'); %'Control_TraceAUC_20sWindow_3Conds.csv'); % neuronal hand click cell scan
%
% %% Load calibration file
% calibration ='E:\matlab\2p-img-analysis\tests\res\calibration_20x.mat';
% CalFile = Calibration_PixelSize.load(calibration);
%
% %% make a mixing matrix for spectral unmixing of RCaMP and GCaMP
%
% if ~exist(fullfile(Settings.MainDir, 'Results','RCaMP_cGCaMP_Matrix.mat'),'file')
%     % load example high res image
%     unmixFile = 'E:\Data\Two_Photon_Data\GCaMP_RCaMP\cyto_GCaMP6s\RG12\Awake\2016_02_04\spot1_long_Nostim\highres_spot1_long096.tif';
%     unmixImg = SCIM_Tif(unmixFile, channel, CalFile);
%
%     [unmixImg, RCaMP_cGCaMP_Matrix] = unmix_chs(unmixImg); %returns the mixing matrix to be used for all imaging
%
%     % save the mixing matrix for loading later
%     cd(fullfile(Settings.MainDir, 'Results'));
%     % write matrix to created file
%     save('RCaMP_cGCaMP_Matrix.mat', 'RCaMP_cGCaMP_Matrix');
% else
%     load(fullfile(Settings.MainDir, 'Results','RCaMP_cGCaMP_Matrix.mat'));
% end
%
%
% %% load scoresheet and loop through animal, spot, etc.
%
% Settings.ScoreSheetPath = fullfile(Settings.MainDir,'Scoresheets',Settings.ScoreSheetNames);
%
% numAnimals = length(Settings.AnimalNames);
% for iAnimal = 1:numAnimals
%     CurrentAnimal = Settings.AnimalNames{iAnimal};
%     %read Scoresheet of current animal
%     CurrentSheet = Settings.ScoreSheetPath{iAnimal};
%     Settings = readScoresheet2(CurrentSheet, Settings);
%
%     % Get SpotID
%     spots = unique(Settings.SpotIDs);
%
%     for iSpot = 1:length(spots) % looks at every 2nd line of SpotIDs-
%
%         % Extract spot name
%         spotId = spots{iSpot};
%
%         % Find the idx of paths matching this spot (all Conditions)
%         matchingSpotsIdx = find(~cellfun(@isempty, regexp(Settings.SpotIDs, spotId)));
%         SpotPaths = Settings.LowresPath(matchingSpotsIdx);
%         CurrentDepth = Settings.Depth(matchingSpotsIdx);
%         CurrentBaseline = Settings.BL_frames(matchingSpotsIdx);
%
%         if length(SpotPaths) ~= length(Settings.NameConditions)
%             error('not all stimulus conditions are present')
%         end
%
%         for iCond = 1:length(Settings.NameConditions)
%             CurrentCondition = Settings.NameConditions{iCond};
%             PathIdx = find(~cellfun(@isempty, regexp(SpotPaths, CurrentCondition)));
%             BL_frames = CurrentBaseline(PathIdx);
%
%             % Get image paths
%             testRoot =SpotPaths{PathIdx};
%
%             if ~exist(fullfile(testRoot,SaveFiles{1,2}),'file') %if the CellScan does not exist, create a new one
%                 expfiles = dir(fullfile(testRoot,'lowres*'));
%                 fnTempList = {expfiles(:).name};
%                 fnList = fullfile(testRoot, fnTempList);
%
%                 % Create an array of ScanImage Tiffs
%                 ImgArray =  SCIM_Tif(fnList, channel, CalFile);
%
%                 %exclude frames at the beginning of each trial where pockel
%                 %cell was dark
%                 ImgArray.exclude_frames('badframes',(1:2),'method', 'inpaint','inpaintIters',5);
%
%
%                 % Spectral Unmixing of GCaMP and RCaMP
%                 ImgArray= ImgArray.unmix_chs(false, [], RCaMP_cGCaMP_Matrix);
%
%                 % Run motion correction
%                 expfiles2 = dir(fullfile(testRoot,'highres*'));
%                 fnTempList2 = {expfiles2(:).name};
%                 RefImgName = cell2mat(fullfile(testRoot, fnTempList2));
%
%                 HighRes = SCIM_Tif(RefImgName,channel, CalFile);
%                 % Extract a reference image
%                 refImg = mean(HighRes.rawdata(:,:,2,:),4);
%                 %figure, imagesc(refImg), axis image, axis off, colormap(gray)
%
%                 ImgArray = ImgArray.motion_correct('refImg',refImg, 'ch',2,'maxShift', 10,'minCorr', 0.4);%,'doPlot',true);
%                 % fill in bad data?
%                 if plotMotion==1
%                     ImgArray(1,1).plot()
%                 end
%
%                 %% Configs for Finding ROIs
%                 %Astrocyte Calcium
%                 % automated selection
%                 findConf{1} = ConfigFindROIsFLIKA_2D.from_preset('ca_cyto_astro', 'baselineFrames',...
%                     BL_frames,'sigmaXY', 2,...
%                     'sigmaT', 0.5,'threshold_std', 5,...
%                     'min_rise_time',0.1689, 'max_rise_time', 8,...
%                     'dilateXY', 4,...
%                     'dilateT', 0.5,'erodeXY', 2);
%
%                 %   freqPassBand=0.15, minPuffArea = 3.7832,
%                 %   minPuffTime = 1.3514, erodeT = 2.1146, threshold_2D = 0.2
%
%                 % Astrocyte Calcium
%                 % hand selected- cellular structures
%                 x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
%                 findConf{2} = ConfigFindROIsDummy.from_ImageJ(fullfile(testRoot,'Astrocytes.zip'), x_pix, y_pix,1);
%
%                 % Neuronal Calcium
%                 % hand selected
%                 x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
%                 findConf{3} = ConfigFindROIsDummy.from_ImageJ(fullfile(testRoot,'Neurons.zip'), x_pix, y_pix,1);
%
%                 %% Configuration for measuring ROIs
%                 % AWAKE astrocyte cyto calcium
%                 measureConf{1} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...'normMethod', 'z-score','zIters', 100,...
%                     'propagateNaNs', false, 'excludeNaNs', false, 'lpWindowTime', 3, 'spFilterOrder', 2,...
%                     'spPassBand', [0.05, 0.2], 'thresholdSD_low', 3,'thresholdSD_band', 5);
%
%                 % AWAKE neuron calcium
%                 measureConf{2} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,... 'normMethod','z-score',... 'zIters', 10000,...
%                     'propagateNaNs', false,'excludeNaNs', false, 'lpWindowTime', 2, 'spFilterOrder', 2,...
%                     'spPassBand', [0.1, 1], 'thresholdSD_low', 3,'thresholdSD_band', 5);
%
%                 % OLD PARAMS
% %                 % AWAKE astrocyte cyto calcium
% %                 measureConf{1} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...'normMethod', 'z-score',...
% %                     'propagateNaNs', false, 'excludeNaNs', false, 'lpWindowTime', 5, 'spFilterOrder', 4,...
% %                     'spPassBand', [0.025, 0.2], 'thresholdSD_low', 3,'thresholdSD_band', 3);
% %
% %                 % AWAKE neuron calcium
% %                 measureConf{2} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...  'normMethod','z-score','zIters', 10000,...
% %                     'propagateNaNs', false,'excludeNaNs', false, 'lpWindowTime', 3, 'spFilterOrder', 2,...
% %                     'spPassBand', [0.1, 1], 'thresholdSD_low', 7,'thresholdSD_band', 5);
%
%                 % for calculating AUC for each trace
%                 detectConf = ConfigMeasureROIsDummy('baselineFrames', BL_frames);
%
%                 %%
%                 % Combine the configs into a CellScan config
%                 configCS{1,1} = ConfigCellScan(findConf{1,1}, measureConf{1,1},detectConf); % astrocyte FLIKA, peaks
%                 configCS{1,2} = ConfigCellScan(findConf{1,2}, measureConf{1,1},detectConf); % astrocyte hand, peaks
%                 configCS{1,3} = ConfigCellScan(findConf{1,3}, measureConf{1,2},detectConf); % neuronal hand, peaks
%
%                 %% Create CellScan objects
%                 CSArray_Ch1_FLIKA = CellScan(fnList, ImgArray, configCS{1,1}, 1);
%
%                 CSArray_Ch1_Hand = CellScan(fnList, ImgArray, configCS{1,2}, 1);
%
%                 CSArray_Ch2_Hand = CellScan(fnList, ImgArray, configCS{1,3}, 2);
%
%                 %% Process the images
%                 CSArray_Ch1_FLIKA =CSArray_Ch1_FLIKA.process();
%                 CSArray_Ch1_Hand =CSArray_Ch1_Hand.process();
%                 CSArray_Ch2_Hand =CSArray_Ch2_Hand.process();
%
%                 %NO IMAGE COMBINATION
%                 % if isfield(channel,'Ca_Cyto_Astro')
%                 % Use this function to combine the masks
%                 %CSArray_Ch1_combined= combine_masks(Ch1_test);
%                 %CSArray_Ch1_combined.process();
%                 %end
%
%                 %% Make the debugging plots
%                 if doplots ==1
%                     CSArray_Ch1_FLIKA.plot;
%                     CSArray_Ch1_Hand.plot;
%                     CSArray_Ch2_Hand.plot;
%
%                     % for working out find peaks parameters
%                     %CSArray_Ch1_FLIKA(1).plot('parameters','plotROIs',6:15);
%                     %CSArray_Ch2_Hand(1).plot('parameters','plotROIs',11:15);
%
%                     % for viewing peaks that were found
%                     %             CSArray_Ch1_FLIKA.plot('signals');
%                     %             CSArray_Ch1_Hand.plot('signals');
%                     %            CSArray_Ch2_Hand.plot('signals');
%                 end
%             else
%                 load (fullfile(testRoot,SaveFiles{1,2}));
%                 load (fullfile(testRoot,SaveFiles{1,3}));
%                 load (fullfile(testRoot,SaveFiles{1,4}));
%             end
%
%             %% Output data
%
%             % make a giant data table
%             listFields = {'amplitude', 'fullWidth', 'halfWidth', ...
%                 'numPeaks', 'peakTime', 'peakStart', 'peakStartHalf', ...
%                 'peakType', 'prominence', 'ROIname', 'peakAUC'};
%
%             % Astrocyte FLIKA
%             for itrial=1:length(CSArray_Ch1_FLIKA)
%
%                 % peak output
%                 temp=CSArray_Ch1_FLIKA(1,itrial).calcDetectSigs.data;
%                 temp2.trialname ={};
%                 temp2.animalname = {};
%                 temp2.channel = {};
%                 temp2.Spot = {};
%                 temp2.Cond = {};
%                 temp2.depth = {};
%                 temp2.overlap = {};
%                 temp2.area = {};
%
%                 % extract fields from Class
%                 for jField = 1:numel(listFields)
%                     isFirst = (itrial == 1 );
%                     if isFirst
%                         data.(listFields{jField}) = {};
%                     end
%                     data.(listFields{jField}) = [data.(listFields{jField}); ...
%                         temp.(listFields{jField})];
%                 end
%
%                 % create fields for trial, animal, spot, condition, etc.
%                 for iPeak = 1:length(temp.amplitude)
%                     temp2.trialname{iPeak,1}=strcat('trial', num2str(itrial));
%                     temp2.channel{iPeak,1}= 'GCaMP';
%                     temp2.Spot{iPeak,1}= spotId;
%                     temp2.animalname{iPeak,1}= CurrentAnimal;
%                     temp2.Cond{iPeak,1} = CurrentCondition;
%                     temp2.depth{iPeak,1} = CurrentDepth(1,1);
%
%                     % get the indices  and area for a particular ROI
%                     jROIname = temp.ROIname{iPeak};
%                     ROIindex= strcmp(CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.roiNames,jROIname);
%                     temp2.area{iPeak,1} = CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.area(ROIindex);
%
%                     ROIpuffIdx = CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.puffIdxs(ROIindex);
%                     x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
%                     [jy,jx,~] = ind2sub([y_pix,x_pix],cell2mat(ROIpuffIdx));
%                     jx = round(mean(jx));
%                     jy = round(mean(jy));
%                     xclose = round(x_pix*0.03); % 3 percent of pixels
%                     yclose = round(y_pix*0.03); % 3 percent of pixels
%
%                     % find astrocyte process and soma ROIs that have
%                     % similar centroids
%                     for kROI= 1:length(CSArray_Ch1_Hand(1,itrial).calcFindROIs.data.roiNames)
%                         % get the indices of the handclicked ROIs
%                         ROIMask{kROI}= CSArray_Ch1_Hand(1, 1).calcFindROIs.data.roiMask(:,:,kROI);
%                         ROI_Idx{kROI} = find(ROIMask{kROI});
%                         % mean pixels for soma ROI
%                         [ky,kx,~] = ind2sub([y_pix,x_pix],ROI_Idx{kROI});
%                         kx = round(mean(kx));
%                         ky = round(mean(ky));
%
%                         % check if they are too close (actually same region)
%                         if jx<=(kx+xclose) && jx>=(kx-xclose) && jy<=(ky+yclose) && jy>=(ky-yclose) % only %3 of pixels apart
%                             spatialcorr(kROI)  = 1;
%                         end
%                     end
%
%                     if exist('spatialcorr','var')
%                         indx = find(spatialcorr>0);
%                         temp2.overlap{iPeak,1}= CSArray_Ch1_Hand(1,itrial).calcFindROIs.data.roiNames{indx};
%                     else
%                         temp2.overlap{iPeak,1} = 0;
%                     end
%                     clear spatialcorr
%                 end
%
%                 isFirst = (itrial == 1 );
%                 if isFirst
%                     data.Trial = {};
%                     data.Animal = {};
%                     data.Channel = {};
%                     data.Spot = {};
%                     data.Condition = {};
%                     data.Depth = {};
%                     data.area = {};
%                     data.overlap = {};
%                 end
%                 data.Trial= [data.Trial; temp2.trialname];
%                 data.Animal= [data.Animal; temp2.animalname];
%                 data.Channel= [data.Channel; temp2.channel];
%                 data.Spot= [data.Spot; temp2.Spot];
%                 data.Condition= [data.Condition; temp2.Cond];
%                 data.Depth= [data.Depth; temp2.depth];
%                 data.area= [data.area; temp2.area];
%                 data.overlap= [data.overlap; temp2.overlap];
%
%                 % trace AUC output
%                 if strcmp(CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.roiNames{1,1}, 'none')
%                     continue
%                 else
%                     traces= CSArray_Ch1_FLIKA(1,itrial).calcMeasureROIs.data.tracesNorm;
%                     %preallocate
%                     AUCdata=cell(size(traces,2),8);
%                     for iROI = 1:size(traces,2)
%                         AUCdata{iROI,1}= trapz(traces(60:60+round(20*11.84),iROI));
%                         AUCdata{iROI,2}=CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.roiNames{iROI,1};
%                         AUCdata{iROI,3}=strcat('trial', num2str(itrial));
%                         AUCdata{iROI,4}= 'GCaMP';
%                         AUCdata{iROI,5}= spotId;
%                         AUCdata{iROI,6}= CurrentAnimal;
%                         AUCdata{iROI,7} = CurrentCondition;
%                         AUCdata{iROI,8} = CurrentDepth(1,1);
%                     end
%                     All_AUC=vertcat(All_AUC, AUCdata);
%                     clearvars AUCdata
%                 end
%             end
%             clearvars temp temp2
%
%             % Astrocyte Hand clicked ROIs
%             for itrial=1:length(CSArray_Ch1_Hand)
%
%                 temp=CSArray_Ch1_Hand(1,itrial).calcDetectSigs.data;
%                 temp2.trialname ={};
%                 temp2.animalname = {};
%                 temp2.channel = {};
%                 temp2.Spot = {};
%                 temp2.Cond = {};
%                 temp2.depth = {};
%                 temp2.overlap = {};
%                 temp2.area = {};
%
%                 % extract fields from Class
%                 for jField = 1:numel(listFields)
%                     data.(listFields{jField}) = [data.(listFields{jField}); ...
%                         temp.(listFields{jField})];
%                 end
%
%                 % create fields for trial, animal, spot, condition, etc.
%                 for iPeak = 1:length(temp.amplitude)
%                     temp2.trialname{iPeak,1}=strcat('trial', num2str(itrial));
%                     temp2.channel{iPeak,1}= 'GCaMP';
%                     temp2.Spot{iPeak,1}= spotId;
%                     temp2.animalname{iPeak,1}= CurrentAnimal;
%                     temp2.Cond{iPeak,1} = CurrentCondition;
%                     temp2.depth{iPeak,1} = CurrentDepth(1,1);
%                     temp2.area{iPeak,1} = 0;
%                     temp2.overlap{iPeak,1} = 0;
%                 end
%
%                 data.Trial= [data.Trial; temp2.trialname];
%                 data.Animal= [data.Animal; temp2.animalname];
%                 data.Channel= [data.Channel; temp2.channel];
%                 data.Spot= [data.Spot; temp2.Spot];
%                 data.Condition= [data.Condition; temp2.Cond];
%                 data.Depth= [data.Depth; temp2.depth];
%                 data.area= [data.area; temp2.area];
%                 data.overlap= [data.overlap; temp2.overlap];
%
%                 % trace AUC output
%                 traces= CSArray_Ch1_Hand(1,itrial).calcMeasureROIs.data.tracesNorm;
%                 %preallocate
%                 AUCdata=cell(size(traces,2),8);
%                 for iROI = 1:size(traces,2)
%                     AUCdata{iROI,1}= trapz(traces(60:60+round(20*11.84),iROI));
%                     AUCdata{iROI,2}=CSArray_Ch1_Hand(1,itrial).calcFindROIs.data.roiNames{iROI,1};
%                     AUCdata{iROI,3}=strcat('trial', num2str(itrial));
%                     AUCdata{iROI,4}= 'GCaMP';
%                     AUCdata{iROI,5}= spotId;
%                     AUCdata{iROI,6}= CurrentAnimal;
%                     AUCdata{iROI,7} = CurrentCondition;
%                     AUCdata{iROI,8} = CurrentDepth(1,1);
%                 end
%                 All_AUC=vertcat(All_AUC, AUCdata);
%                 clearvars AUCdata
%             end
%             clearvars temp temp2
%
%             % Neuronal Hand clicked ROIs
%             for itrial=1:length(CSArray_Ch2_Hand)
%
%                 temp=CSArray_Ch2_Hand(1,itrial).calcDetectSigs.data;
%                 temp2.trialname ={};
%                 temp2.animalname = {};
%                 temp2.channel = {};
%                 temp2.Spot = {};
%                 temp2.Cond = {};
%                 temp2.depth = {};
%                 temp2.overlap = {};
%                 temp2.area = {};
%
%                 % extract fields from Class
%                 for jField = 1:numel(listFields)
%                     data.(listFields{jField}) = [data.(listFields{jField}); ...
%                         temp.(listFields{jField})];
%                 end
%
%                 % create fields for trial, animal, spot, condition, etc.
%                 for iPeak = 1:length(temp.amplitude)
%                     temp2.trialname{iPeak,1}=strcat('trial', num2str(itrial));
%                     temp2.channel{iPeak,1}= 'RCaMP';
%                     temp2.Spot{iPeak,1}= spotId;
%                     temp2.animalname{iPeak,1}= CurrentAnimal;
%                     temp2.Cond{iPeak,1} = CurrentCondition;
%                     temp2.depth{iPeak,1} = CurrentDepth(1,1);
%                     temp2.area{iPeak,1} = 0;
%                     temp2.overlap{iPeak,1} = 0;
%                 end
%
%                 data.Trial= [data.Trial; temp2.trialname];
%                 data.Animal= [data.Animal; temp2.animalname];
%                 data.Channel= [data.Channel; temp2.channel];
%                 data.Spot= [data.Spot; temp2.Spot];
%                 data.Condition= [data.Condition; temp2.Cond];
%                 data.Depth= [data.Depth; temp2.depth];
%                 data.area= [data.area; temp2.area];
%                 data.overlap= [data.overlap; temp2.overlap];
%
%                 % trace AUC output
%                 traces= CSArray_Ch2_Hand(1,itrial).calcMeasureROIs.data.tracesNorm;
%                 %preallocate
%                 AUCdata=cell(size(traces,2),8);
%                 for iROI = 1:size(traces,2)
%                     AUCdata{iROI,1}= trapz(traces(60:60+round(20*11.84),iROI));
%                     AUCdata{iROI,2}=CSArray_Ch2_Hand(1,itrial).calcFindROIs.data.roiNames{iROI,1};
%                     AUCdata{iROI,3}=strcat('trial', num2str(itrial));
%                     AUCdata{iROI,4}= 'RCaMP';
%                     AUCdata{iROI,5}= spotId;
%                     AUCdata{iROI,6}= CurrentAnimal;
%                     AUCdata{iROI,7} = CurrentCondition;
%                     AUCdata{iROI,8} = CurrentDepth(1,1);
%                 end
%                 All_AUC=vertcat(All_AUC, AUCdata);
%                 clearvars AUCdata
%             end
%
%             dataNames=fieldnames(data);
%             data2= struct2cell(data);
%             data3= [data2{:}];
%
%             AllData=vertcat(AllData, data3);
%
%             if ~exist(fullfile(testRoot,SaveFiles{1,2}),'file')
%                 %Save CellScan
%               %  save(fullfile(testRoot,SaveFiles{1,2}), 'CSArray_Ch1_FLIKA','-v7.3');
%              %   save(fullfile(testRoot,SaveFiles{1,3}), 'CSArray_Ch1_Hand','-v7.3');
%                % save(fullfile(testRoot,SaveFiles{1,4}), 'CSArray_Ch2_Hand','-v7.3');
%             end
%
%             clearvars data data3 temp temp2
%         end
%     end
% end
%
% % %% Save all data for R analysis
% AllData2= [dataNames';AllData];
%
% AUCNames= {'AUC10s','ROI','Trial','Channel','Spot','Animal','Condition','depth'};
% All_AUC2=[AUCNames;All_AUC];
% %
% cd(fullfile(Settings.MainDir, 'Results'));
% % write date to created file
% cell2csv(SaveFiles{1,1}, AllData2);
% cell2csv(SaveFiles{1,5}, All_AUC2);
%
