clearvars

AllData= [];
All_traces= [];


%% Information about your images
Settings.MainDir = 'E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f';

Settings.AnimalNames = {
        'RG14',...
        'RG16',...
        'RG17',...
        'RG18',...
    };
Settings.ScoreSheetNames = {
        'RG14_Scoresheet_LongTrials.xls',...
        'RG16_Scoresheet_LongTrials.xls',...
        'RG17_Scoresheet_LongTrials.xls',...
        'RG18_Scoresheet_LongTrials.xls',...
    };
Settings.NameConditions = {'Nostim','Stim'};
%Settings.NameConditions = {'Stim'};

channel = struct('Ca_Memb_Astro',1,'Ca_Neuron',2);

plotMotion = 0; %Plot motion correction movie
doplots = 0; %Plots for each trial

% final data file name
SaveFiles{1,1} = fullfile(Settings.MainDir, 'Results', 'LStim_LckGC&RC_17_02_2017.csv');%'Control_Peaks_3Conds.csv'); % all data
SaveFiles{1,2} = fullfile(Settings.MainDir, 'Results', 'LStim_LckGC&RC_17_02_2017.mat');%'Control_Peaks_3Conds.csv'); % all data
SaveFiles{1,3}= fullfile(Settings.MainDir, 'Results','LStim_LckGC&RC_traces_17_02_2017.mat'); %'Control_TraceAUC_20sWindow_3Conds.csv'); % neuronal hand click cell scan

%% Load calibration file
calibration ='E:\matlab\2p-img-analysis\tests\res\calibration_20x.mat';
CalFile = Calibration_PixelSize.load(calibration);

%% make a mixing matrix for spectral unmixing of RCaMP and GCaMP

if ~exist(fullfile(Settings.MainDir, 'Results','RCaMP_mGCaMP_Matrix.mat'),'file')
    % load example high res image
    unmixFile = 'E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\RG14\Awake\test\highres_spot1_long021.tif';
    unmixImg = SCIM_Tif(unmixFile, channel, CalFile);
    
    [unmixImg, RCaMP_mGCaMP_Matrix] = unmix_chs(unmixImg); %returns the mixing matrix to be used for all imaging
    
    % save the mixing matrix for loading later
    cd(fullfile(Settings.MainDir, 'Results'));
    % write matrix to created file
    save('RCaMP_mGCaMP_Matrix.mat', 'RCaMP_mGCaMP_Matrix');
else
    load(fullfile(Settings.MainDir, 'Results','RCaMP_mGCaMP_Matrix.mat'));
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
                ImgArray= ImgArray.unmix_chs(false, [], cell2mat(RCaMP_mGCaMP_Matrix));
                
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
                    ImgArray(1,3).plot();
                end
                
                %% Configs for Finding ROIs
                %Astrocyte Calcium
                % automated selection
                findConf{1} = ConfigFindROIsFLIKA_2D.from_preset('ca_memb_astro', 'baselineFrames',...
                    BL_frames,'freqPassBand',1,'sigmaXY', 2,...
                    'sigmaT', 0.1,'threshold_std', 7, 'threshold2D', 0.2,...
                    'min_rise_time',0.0845, 'max_rise_time', 1,'minPuffArea', 10,...
                    'dilateXY', 5, 'dilateT', 0.3,'erodeXY', 1, 'erodeT', 0.1,...
                    'discardBorderROIs',true);
                
                
                % Astrocyte Calcium
                % hand selected- cellular structures
                x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
                findConf{2} = ConfigFindROIsDummy.from_ImageJ(fullfile(testRoot,'Astrocytes.zip'), x_pix, y_pix,1);
                
                % Neuronal Calcium
                % hand selected
                x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
                findConf{3} = ConfigFindROIsDummy.from_ImageJ(fullfile(testRoot,'Neurons.zip'), x_pix, y_pix,1);
                
                % FLIKA selected
                    findConf{4} = ConfigFindROIsFLIKA_3D.from_preset('ca_neuron', 'baselineFrames',...
                    BL_frames,'freqPassBand',1,'sigmaXY', 2,...
                    'sigmaT', 0.1,'threshold_std', 10, 'threshold2D', 0.2,...
                    'min_rise_time',0.0845, 'max_rise_time', 1,'minPuffArea', 10,...
                    'minPuffTime',0.25,'dilateXY', 2, 'dilateT', 0.3,'erodeXY', 1, 'erodeT', 0.1,...
                    'discardBorderROIs',true);
                %% Configuration for measuring ROIs
                % AWAKE astrocyte membrane calcium
                measureConf{1} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...'normMethod', 'z-score','zIters', 100,...
                    'propagateNaNs', false, 'excludeNaNs', false, 'lpWindowTime', 1.5, 'spFilterOrder', 2,...
                    'spPassBandMin',0.05, 'spPassBandMax', 0.5, 'thresholdSD_low', 3,'thresholdSD_band', 5);
                
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
                configCS{1,4} = ConfigCellScan(findConf{1,4}, measureConf{1,2},detectConf); % neuronal hand, peaks
                
                %% Create CellScan objects
                CSArray_Ch1_FLIKA = CellScan(fnList, ImgArray, configCS{1,1}, 1);
                
                CSArray_Ch1_Hand = CellScan(fnList, ImgArray, configCS{1,2}, 1);
                
                CSArray_Ch2_Hand = CellScan(fnList, ImgArray, configCS{1,3}, 2);
                
                CSArray_Ch2_FLIKA = CellScan(fnList, ImgArray, configCS{1,4}, 2);
                
                %% Process the images
                CSArray_Ch1_FLIKA =CSArray_Ch1_FLIKA.process();
                CSArray_Ch1_Hand =CSArray_Ch1_Hand.process();
                CSArray_Ch2_Hand =CSArray_Ch2_Hand.process();
               % CSArray_Ch2_FLIKA =CSArray_Ch2_FLIKA.process();
                
                %NO IMAGE COMBINATION
                % if isfield(channel,'Ca_Cyto_Astro')
                % Use this function to combine the masks
                %CSArray_Ch1_combined= combine_masks(Ch1_test);
                %CSArray_Ch1_combined.process();
                %end
                
                %% Make the debugging plots
                   % CSArray_Ch2_FLIKA.plot;
%                     CSArray_Ch1_Hand.plot;
%                     CSArray_Ch2_Hand.plot;
                    
                    % for working out find peaks parameters
                    %CSArray_Ch2_FLIKA(1).opt_config()


            
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
                temp2.pixelsize = {};
                
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
                    temp2.pixelsize{iPeak,1} = CSArray_Ch1_FLIKA(1,itrial).rawImg.metadata.pixelSize;
                    
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
                    data.pixelsize={};
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
                        Trace_data{iROI,11} = CSArray_Ch1_FLIKA(1,itrial).rawImg.metadata.pixelSize;
                    
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
                            Trace_data{iROI,12}= CSArray_Ch1_Hand(1,itrial).calcFindROIs.data.roiNames{indx};
                        else
                            Trace_data{iROI,12} = 0;
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
                temp2.pixelsize = {};
                
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
                    temp2.pixelsize{iPeak,1} = CSArray_Ch1_FLIKA(1,itrial).rawImg.metadata.pixelSize;
                    
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
                    Trace_data{iROI,12} = 0;
                    Trace_data{iROI,11} = CSArray_Ch1_Hand(1,itrial).rawImg.metadata.pixelSize;
                    
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
                temp2.pixelsize = {};
                
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
                    temp2.pixelsize{iPeak,1} = CSArray_Ch1_FLIKA(1,itrial).rawImg.metadata.pixelSize;
                    
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
                    Trace_data{iROI,12} = 0;
                    Trace_data{iROI,11} = CSArray_Ch2_Hand(1,itrial).rawImg.metadata.pixelSize;
                    
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










%% Lck

% all ROIs  correlation
clearvars

doplot=0;

FileSave{1,1}='E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LongStim_Correlations.mat';
FileSave{1,2}='E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LongStim_Correlations.csv';

% load data traces
load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\S&LStim_LckGC&RC_traces_17_02_2017.mat');

Long_Short=All_traces;

load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LStim_LckGC&RC_traces_17_02_2017.mat');

Long=All_traces;

All_traces=vertcat(Long_Short,Long);

FrameRate=11.84;

searchString='Stim';

for iROI=1:length(All_traces)
    Stim_str{iROI} = strfind(All_traces{iROI,6},searchString);
    str_idx(iROI)=~isempty(Stim_str{iROI});
    
    %find ROITypes
    N_str= strfind(All_traces{iROI, 1},'N');
    S_str= strfind(All_traces{iROI, 1},'S');
    EF_str= strfind(All_traces{iROI, 1},'E');
    P_str= strfind(All_traces{iROI, 1},'r');
    NP_str= strfind(All_traces{iROI, 1},'np');
    
    if ~isempty(N_str)
        All_traces{iROI,13}='Neuron';
    elseif ~isempty(S_str)
        All_traces{iROI,13}='Soma';
    elseif ~isempty(EF_str)
        All_traces{iROI,13}='Endfeet';
    elseif ~isempty(P_str)
        All_traces{iROI,13}='Process';
    elseif ~isempty(NP_str)
        All_traces{iROI,13}='Neuropil';
    end
    
    % make new names
    All_traces{iROI,14}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,1});
    All_traces{iROI,15}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,2}, '_',All_traces{iROI,1});
    All_traces{iROI,16}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,2}, '_',All_traces{iROI,1},'_',All_traces{iROI,6});
    All_traces{iROI,17}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,2}, '_',All_traces{iROI,6});
    
end

data_traces = All_traces(str_idx',:);

for iROI=1:length(data_traces)
    %if a process ROI is overlapping with a soma or process, exclude it
    %from data
    %All_traces{iROI,11}=0;
    nonOverlapIdx(iROI)=~ischar(data_traces{iROI,12});
    
end

% remove overlapping processes
data_traces = data_traces(nonOverlapIdx',:);

%% Correlation of each ROI with all ROIs in the same trial
% CONSIDER WHOLE TRIAL and 25 sec stim window

% only compare ROIs from the same trial together
trial_condition=unique(data_traces(:,17));
All_Corrs=[];
for itrial= 1:length(trial_condition)
    CurrentTrial=trial_condition{itrial};
    
    for iTrace= 1:size(data_traces,1)
        TrialIdx(iTrace)=strcmp(data_traces{iTrace,17}, CurrentTrial);
    end
    
    CurrentData=data_traces(TrialIdx,:);
    
    for iROI=1:size(CurrentData,1)
        CorrData(1,1:6)=CurrentData(iROI,2:7);
        CorrData(1,7)=CurrentData(iROI,17);
        
        CorrData(1,8)= CurrentData(iROI,1);
        CorrData(1,9)= CurrentData(iROI,13);
        CorrData(1,10)= CurrentData(iROI,9);
        CorrData(1,11)= CurrentData(iROI,10);
        
        %whole trace
        TraceX= CurrentData{1,8};
        % stim window
        TraceX_short=TraceX(1:(25*FrameRate),:);
        
        for kROI = (iROI+1):size(CurrentData,1)
            CorrData(1,12)= CurrentData(kROI,3);
            CorrData(1,13)= CurrentData(kROI,1);
            CorrData(1,14)= CurrentData(kROI,13);
            CorrData(1,15)= CurrentData(kROI,9);
            CorrData(1,16)= CurrentData(kROI,10);
            
            %whole trace
            TraceY= CurrentData{kROI,8};
            % stim window
            TraceY_short=TraceY(1:(25*FrameRate),:);
            
            % whole trial correlation
            [R,Ps] = corrcoef(TraceX, TraceY);
            CorrData{1,17}= R(1,2); % correlation coeff
            CorrData{1,18}= Ps(1,2); % p value
            
            % stim window correlation
            [R2,Ps2]= corrcoef(TraceX_short, TraceY_short);
            CorrData{1,19}= R2(1,2); % correlation coeff
            CorrData{1,20}= Ps2(1,2); % p value
            
            % stim window cross correlation
            [Corr_Sequence,lag]  = xcorr(TraceX_short, TraceY_short,'coeff');
            [~,I] = max(abs(Corr_Sequence));
            CorrData{1,21}=  max(abs(Corr_Sequence)); % max of correlation sequence
            CorrData{1,22}= (lag(I))/FrameRate; % lag time
            
            % correlation plots
            %             subplot(3,1,1); plot(TraceX_short); title('s1');
            %             subplot(3,1,2); plot(TraceY_short); title('s2');
            %             subplot(3,1,3); plot(lag,Corr_Sequence);
            %             title('Cross-correlation between s1 and s2')
            
            %% Distance Calculations
            % find the minimium distance between the edges of the two ROIs
            
            %create a binary image with each pair of ROIs
            % for the first ROI
            if isnumeric(CorrData{1,11})
                Image1=zeros(128,128);
                Image1(CorrData{1,11})=1;
                %Image1=im2bw(Image1);
            elseif islogical(CorrData{1,11})
                Image1= double(CorrData{1,11});
            else
                Image1=[];
            end
            
            % for the second ROI
            if isnumeric(CorrData{1,16})
                Image2=zeros(128,128);
                Image2(CorrData{1,16})=1;
                Image2=im2bw(Image2);
            elseif islogical(CorrData{1,16})
                Image2= double(CorrData{1,16});
            else
                Image2=[];
            end
            
            Mask=Image1+Image2;
            Mask=im2bw(Mask);
            
            %Pythaogrean theorem method
            
            % Define object boundaries
            boundaries = bwboundaries(Mask);
            numberOfBoundaries = size(boundaries, 1);
            if numberOfBoundaries==1
                CorrData{1,23} = 0;
            elseif numberOfBoundaries>1
                boundary1 = boundaries{1};
                boundary2 = boundaries{2};
                boundary1x = boundary1(:, 2);
                boundary1y = boundary1(:, 1);
                for k = 1 : length(boundary2)
                    boundary2x = boundary2(k, 2);
                    boundary2y = boundary2(k, 1);
                    % For this blob, compute distances from boundaries to edge.
                    allDistances = sqrt((boundary1x - boundary2x).^2 + (boundary1y - boundary2y).^2);
                    % Find closest point, min distance.
                    [minDistance(k), indexOfMin] = min(allDistances);
                end
                % Find the overall min distance
                CorrData{1,23} = (min(minDistance)*cell2mat(CurrentData(kROI,11)));
            else
                CorrData{1,23} = [];
            end
            CorrData{1,24} = numberOfBoundaries; % number of ROIs in the mask
            
            All_Corrs=vertcat(All_Corrs,CorrData); % concatenate all data into a big matrix
            
            clear boundaries boundary1 boundary2 boundary1x boundary1y boundary2x boundary2y minDistance indexofMin
        end
        
    end
    
end




%% save as a CSV
names = {'Trial','ChannelX','Spot','Animal','Condition','Depth','TrialName','ROI_X',...
    'ROI_X_type','ChannelY','ROI_Y','ROI_Y_type','Long_Corr','LongPvalue',...
    'Short_Corr','ShortPvalue','xCorr','Lag','MinDistance','numROIs'};

CorrelationData=All_Corrs;
CorrelationData(:,10) = [];
CorrelationData(:,10) = [];
CorrelationData(:,13) = [];
CorrelationData(:,13) = [];

AllData=[names;CorrelationData];

% write date to created file
cell2csv(FileSave{1,2}, AllData);
save(FileSave{1,1}, 'All_Corrs','-v7.3');
