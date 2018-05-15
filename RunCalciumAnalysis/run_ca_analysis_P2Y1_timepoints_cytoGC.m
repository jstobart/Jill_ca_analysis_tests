%% 6 cytoGCaMP long stim

clearvars

AllData= [];
All_traces= [];


%% Information about your images
Settings.MainDir = 'G:\ZurichData\Astrocyte_Calcium\P2Y1_Mice';

Settings.AnimalNames = {
    'PYRG7',...
    };
Settings.ScoreSheetNames = {
    'PYRG7_Scoresheet_cytoGCaMP.xls',...
    };
Settings.NameConditions = {'WPstim','nostim'};

channel = struct('Ca_Cyto_Astro',1,'Ca_Neuron',2);
plotMotion =0;

% final data file name
SaveFiles{1,1} = fullfile(Settings.MainDir, 'Results', 'cytoGC_P2Y1_timepoints_peaks_05_2018.csv');%'Control_Peaks_3Conds.csv'); % all data
SaveFiles{1,2} = fullfile(Settings.MainDir, 'Results', 'cytoGC_P2Y1_timepoints_peaks_05_2018.mat');%'Control_Peaks_3Conds.csv'); % all data
SaveFiles{1,3}= fullfile(Settings.MainDir, 'Results','cytoGC_P2Y1_timepoints_traces_05_2018.mat'); %'Control_TraceAUC_20sWindow_3Conds.csv'); % neuronal hand click cell scan

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
        
        % Get TimePoints
        CurrentTPs = Settings.TimePoints(matchingSpotsIdx);
        TimePoints = unique(CurrentTPs);
        
        for iTime=1:length(TimePoints)
            % Extract spot name
            spotId = spots{iSpot};
            
            CurrentDepth = Settings.Depth(matchingSpotsIdx);
            CurrentBaseline = Settings.BL_frames(matchingSpotsIdx);
            
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
                %ASTROCYTES
                % 2D automated selection for peaks
                AC_findConf{1} = ConfigFindROIsFLIKA_2D.from_preset('ca_cyto_astro', 'baselineFrames',...
                    BL_frames,'sigmaXY', 2,...
                    'sigmaT', 0.5,'threshold_std', 5,...
                    'min_rise_time',0.1689, 'max_rise_time', 8,...
                    'dilateXY', 4,...
                    'dilateT', 0.5,'erodeXY', 2);
                
                % hand selected- peaks from cellular structures
                x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
                AC_findConf{2} = ConfigFindROIsDummy.from_ImageJ(fullfile(testRoot,'Astrocytes.zip'), x_pix, y_pix,1);
                
                
                % 3D automated selection for time and space estimations
                AC_findConf{3} = ConfigFindROIsFLIKA_3D.from_preset('ca_cyto_astro', 'baselineFrames',...
                    BL_frames,'sigmaXY', 2,...
                    'sigmaT', 0.5,'threshold_std', 5,...
                    'min_rise_time',0.1689, 'max_rise_time', 8,...
                    'dilateXY', 4,...
                    'dilateT', 0.5,'erodeXY', 2);
                
                % NEURONS
                % 2D FLIKA selected for peaks from "dendrites"
                Neur_findConf{1} = ConfigFindROIsFLIKA_2D.from_preset('ca_neuron', 'baselineFrames',...
                    BL_frames,'freqPassBand',1,'sigmaXY', 2,...
                    'sigmaT', 0.1,'threshold_std', 7, 'threshold_2D', 0.2,...
                    'min_rise_time',0.0845, 'max_rise_time', 2,'minPuffArea', 10,...
                    'dilateXY', 3, 'dilateT', 0.5,'erodeXY', 2, 'erodeT', 0.3,...
                    'discardBorderROIs',true);
                
                % hand selected for peaks from somata
                x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
                Neur_findConf{2} = ConfigFindROIsDummy.from_ImageJ(fullfile(testRoot,'Neurons.zip'), x_pix, y_pix,1);
                
                % 3D FLIKA selected for time and space estimations
                Neur_findConf{3} = ConfigFindROIsFLIKA_3D.from_preset('ca_neuron', 'baselineFrames',...
                    BL_frames,'freqPassBand',1,'sigmaXY', 2,...
                    'sigmaT', 0.1,'threshold_std', 7, 'threshold_2D', 0.2,...
                    'min_rise_time',0.0845, 'max_rise_time', 2,'minPuffArea', 10,...
                    'dilateXY', 3, 'dilateT', 0.5,'erodeXY', 2, 'erodeT', 0.3,...
                    'discardBorderROIs',true);
                
                %% Configuration for measuring ROIs
                % AWAKE astrocyte cyto calcium
                detectConf{1} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...'normMethod', 'z-score','zIters', 100,...
                    'propagateNaNs', false, 'excludeNaNs', false, 'lpWindowTime', 3, 'spFilterOrder', 2,...
                    'spPassBandMin',0.05, 'spPassBandMax', 0.2, 'thresholdSD_low', 3,'thresholdSD_band', 5);
                
                % AWAKE neuron calcium
                detectConf{2} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,... 'normMethod','z-score',... 'zIters', 10000,...
                    'propagateNaNs', false,'excludeNaNs', false, 'lpWindowTime', 2, 'spFilterOrder', 2,...
                    'spPassBandMin',0.1, 'spPassBandMax', 1, 'thresholdSD_low', 3,'thresholdSD_band', 5);
                
                % for 3D FLIKA
                detectConf{3} = ConfigDetectSigsDummy();
                
                % for calculating AUC for each trace
                measureConf = ConfigMeasureROIsDummy('baselineFrames', BL_frames);
                
                %%
                % Combine the configs into a CellScan config
                configCS{1,1} = ConfigCellScan(AC_findConf{1,1}, measureConf,detectConf{1,1}); % astrocyte FLIKA, peaks
                configCS{1,2} = ConfigCellScan(AC_findConf{1,2}, measureConf,detectConf{1,1}); % astrocyte hand, peaks
                configCS{1,3} = ConfigCellScan(AC_findConf{1,3}, measureConf,detectConf{1,3}); % astrocyte 3D FLIKA
                
                %% Create CellScan objects
                CSArray_Ch1_FLIKA = CellScan(fnList, ImgArray, configCS{1,1}, 1);
                
                CSArray_Ch1_Hand = CellScan(fnList, ImgArray, configCS{1,2}, 1);
                
                
                
                %% Process the images
                CSArray_Ch1_FLIKA =CSArray_Ch1_FLIKA.process();
                CSArray_Ch1_Hand =CSArray_Ch1_Hand.process();
                
                
                %% Make the debugging plots
                % CSArray_Ch2_FLIKA.plot;
                %                     CSArray_Ch1_Hand.plot;
                %                     CSArray_Ch1_FLIKA.plot;
                
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
                        temp2.trialname{iPeak,1}=strcat('trial', num2str(itrial,'%02d'));
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
                            Trace_data{iROI,2}= strcat('trial', num2str(itrial,'%02d'));
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
                        temp2.trialname{iPeak,1}=strcat('trial', num2str(itrial,'%02d'));
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
                        Trace_data{iROI,2}= strcat('trial', num2str(itrial,'%02d'));
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


