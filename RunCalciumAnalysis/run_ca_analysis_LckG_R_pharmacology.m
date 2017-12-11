
%% 1 lck no stim vs long stim PHARMACOLOGY
clearvars

AllData= [];
All_traces= [];


%% Information about your images
Settings.MainDir = 'E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f';

Settings.AnimalNames = {
    'ARG2',...
    'IPRG2',...
    'IPRG3',...
    
    };
Settings.ScoreSheetNames = {
    'ARG2_Scoresheet_Pharmacology.xls',...
    'IPRG2_Scoresheet_Pharmacology.xls',...
    'IPRG3_Scoresheet_Pharmacology.xls',...
    };
Settings.NameConditions = {'Nostim','Stim'};

channel = struct('Ca_Memb_Astro',1,'Ca_Neuron',2);
plotMotion =0;

% final data file name
SaveFiles{1,1} = fullfile(Settings.MainDir, 'Results','FilesforR', 'Peaks_pharmacology_Lck_nostim_vs_longstim_12_2017.csv');%'Control_Peaks_3Conds.csv'); % all data
SaveFiles{1,2} = fullfile(Settings.MainDir, 'Results','FilesforMatlab', 'Peaks_pharmacology_Lck_nostim_vs_longstim_12_2017.mat');%'Control_Peaks_3Conds.csv'); % all data
SaveFiles{1,3}= fullfile(Settings.MainDir, 'Results','FilesforMatlab','Traces_pharmacology_Lck_nostim_vs_longstim_12_2017.mat'); %'Control_TraceAUC_20sWindow_3Conds.csv'); % neuronal hand click cell scan
SaveFiles{1,4}= fullfile(Settings.MainDir, 'Results','FilesforR','OnsetTimes_pharmacology_Lck_nostim_vs_longstim_12_2017.csv');

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
    if Settings.PMTnum(1,1)==2
        load(fullfile(Settings.MainDir, 'Results','RCaMP_mGCaMP_Matrix_2Ch.mat'));
    elseif Settings.PMTnum(1,1)==4
        load(fullfile(Settings.MainDir, 'Results','RCaMP_mGCaMP_Matrix_4Ch.mat'));
    end
    
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
                
                %exclude frames at the beginning of each trial where pockel
                %cell was dark
                %ImgArray.exclude_frames('badframes',(1:2),'method', 'inpaint','inpaintIters',5);
                
                
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
                %ASTROCYTES
                % 2D automated selection for peaks
                AC_findConf{1} = ConfigFindROIsFLIKA_2D.from_preset('ca_memb_astro', 'baselineFrames',...
                    BL_frames,'freqPassBand',1,'sigmaXY', 2,...
                    'sigmaT', 0.1,'thresholdPuff', 7, 'threshold2D', 0.2,...
                    'minRiseTime',0.0845, 'maxRiseTime', 1,'minROIArea', 10,...
                    'dilateXY', 5, 'dilateT', 0.3,'erodeXY', 1, 'erodeT', 0.1,...
                    'discardBorderROIs',true);
                
                % hand selected- peaks from cellular structures
                x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
                AC_findConf{2} = ConfigFindROIsDummy.from_ImageJ(fullfile(testRoot,'Astrocytes.zip'), x_pix, y_pix,1);
                
                
                % 3D automated selection for time and space estimations
                AC_findConf{3} = ConfigFindROIsFLIKA_3D.from_preset('ca_memb_astro', 'baselineFrames',...
                    BL_frames,'freqPassBand',1,'sigmaXY', 2,...
                    'sigmaT', 0.1,'thresholdPuff', 7, 'threshold2D', 0.2,...
                    'minRiseTime',0.0845, 'maxRiseTime', 1,'minROIArea', 10,...
                    'dilateXY', 5, 'dilateT', 0.3,'erodeXY', 1, 'erodeT', 0.1,...
                    'discardBorderROIs',true);
                
                % NEURONS
                % 2D FLIKA selected for peaks from "dendrites"
                Neur_findConf{1} = ConfigFindROIsFLIKA_2D.from_preset('ca_neuron', 'baselineFrames',...
                    BL_frames,'freqPassBand',1,'sigmaXY', 2,...
                    'sigmaT', 0.1,'thresholdPuff', 7, 'threshold2D', 0.2,...
                    'minRiseTime',0.0845, 'maxRiseTime', 2,'minROIArea', 10,...
                    'dilateXY', 3, 'dilateT', 0.5,'erodeXY', 2, 'erodeT', 0.3,...
                    'discardBorderROIs',true);
                
                % hand selected for peaks from somata
                x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
                Neur_findConf{2} = ConfigFindROIsDummy.from_ImageJ(fullfile(testRoot,'Neurons.zip'), x_pix, y_pix,1);
                
                % 3D FLIKA selected for time and space estimations
                Neur_findConf{3} = ConfigFindROIsFLIKA_3D.from_preset('ca_neuron', 'baselineFrames',...
                    BL_frames,'freqPassBand',1,'sigmaXY', 2,...
                    'sigmaT', 0.1,'thresholdPuff', 7, 'threshold2D', 0.2,...
                    'minRiseTime',0.0845, 'maxRiseTime', 2,'minROIArea', 10,...
                    'dilateXY', 3, 'dilateT', 0.5,'erodeXY', 2, 'erodeT', 0.3,...
                    'discardBorderROIs',true);
                
                %% Configuration for measuring ROIs
                % AWAKE astrocyte membrane calcium
                detectConf{1} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...'normMethod', 'z-score','zIters', 100,...
                    'propagateNaNs', false, 'excludeNaNs', false, 'lpWindowTime', 1.5, 'spFilterOrder', 2,...
                    'spPassBandMin',0.05, 'spPassBandMax', 0.5, 'thresholdLP', 3,'thresholdSP', 5);
                
                % AWAKE neuron calcium
                detectConf{2} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,... 'normMethod','z-score',... 'zIters', 10000,...
                    'propagateNaNs', false,'excludeNaNs', false, 'lpWindowTime', 2, 'spFilterOrder', 2,...
                    'spPassBandMin',0.1, 'spPassBandMax', 1, 'thresholdLP', 3,'thresholdSP', 5);
                
                % for 3D FLIKA
                detectConf{3} = ConfigDetectSigsDummy();
                
                % for calculating AUC for each trace
                measureConf = ConfigMeasureROIsDummy('baselineFrames', BL_frames);
                
                %%
                % Combine the configs into a CellScan config
                configCS{1,1} = ConfigCellScan(AC_findConf{1,1}, measureConf,detectConf{1,1}); % astrocyte FLIKA, peaks
                configCS{1,2} = ConfigCellScan(AC_findConf{1,2}, measureConf,detectConf{1,1}); % astrocyte hand, peaks
                configCS{1,3} = ConfigCellScan(AC_findConf{1,3}, measureConf,detectConf{1,3}); % astrocyte 3D FLIKA
                configCS{1,4} = ConfigCellScan(Neur_findConf{1,1}, measureConf,detectConf{1,2}); % neuronal FLIKA, peaks
                configCS{1,5} = ConfigCellScan(Neur_findConf{1,2}, measureConf,detectConf{1,2}); % neuronal hand peaks
                configCS{1,6} = ConfigCellScan(Neur_findConf{1,3}, measureConf,detectConf{1,3}); % neuronal 3D FLIKA
                
                %% Create CellScan objects
                CSArray_Ch1_FLIKA = CellScan(fnList, ImgArray, configCS{1,1}, 1);
                
                CSArray_Ch1_Hand = CellScan(fnList, ImgArray, configCS{1,2}, 1);
                
                CSArray_Ch2_FLIKA = CellScan(fnList, ImgArray, configCS{1,4}, 2);
                
                CSArray_Ch2_Hand = CellScan(fnList, ImgArray, configCS{1,5}, 2);
                
                
                
                %% Process the images
                CSArray_Ch1_FLIKA =CSArray_Ch1_FLIKA.process();
                CSArray_Ch1_Hand =CSArray_Ch1_Hand.process();
                CSArray_Ch2_Hand =CSArray_Ch2_Hand.process();
                CSArray_Ch2_FLIKA =CSArray_Ch2_FLIKA.process();
                
                %             CSArray_Ch1_FLIKA.plot();
                %             CSArray_Ch2_FLIKA.plot();
                
                % CSArray_Ch1_FLIKA.opt_config();
                %% Output data
                
                % make a giant data table
                listFields = {'amplitude', 'fullWidth', 'halfWidth', ...
                    'numPeaks', 'peakTime', 'peakStart', 'peakStartHalf', ...
                    'peakType', 'prominence', 'roiName', 'peakAUC'};
                
                CellScans=vertcat(CSArray_Ch1_FLIKA, CSArray_Ch1_Hand,...
                    CSArray_Ch2_FLIKA, CSArray_Ch2_Hand);
                
                
                % loop through cellscans
                for iScan=1:size(CellScans,1)
                    for itrial=1:size(CellScans,2)
                        
                        % get fraction of active MD area
                        neuroMask = CellScans(4,itrial).calcFindROIs.data.roiMask;
                        if ndims(neuroMask) == 3
                            neuroMask = max(neuroMask, [], 3);
                        end
                        [nFluoPix, nActivePix, nTotalPix] = ...
                            getFracActive(CellScans(1,itrial), 'nanMask', ...
                            neuroMask);
                        
                        % peak output
                        temp=CellScans(iScan,itrial).calcDetectSigs.data;
                        temp2.trialname ={};
                        temp2.animalname = {};
                        temp2.drug={};
                        temp2.channel = {};
                        temp2.Spot = {};
                        temp2.Cond = {};
                        temp2.depth = {};
                        temp2.overlap = {};
                        temp2.area = {};
                        temp2.pixelsize = {};
                        temp2.nFluoPix = {};
                        temp2.nActivePix = {};
                        
                        % extract fields from Class
                        for jField = 1:numel(listFields)
                            isFirst = (itrial == 1 && iScan == 1);
                            if isFirst
                                data.(listFields{jField}) = {};
                            end
                            data.(listFields{jField}) = [data.(listFields{jField}); ...
                                temp.(listFields{jField})];
                        end
                        
                        % create fields for trial, animal, spot, condition, etc.
                        for iPeak = 1:length(temp.amplitude)
                            temp2.trialname{iPeak,1}=strcat('trial', num2str(itrial,'%02d'));
                            if iScan==1 || iScan==2
                                temp2.channel{iPeak,1}= 'GCaMP';
                            else
                                temp2.channel{iPeak,1}= 'RCaMP';
                            end
                            temp2.Spot{iPeak,1}= spotId;
                            temp2.animalname{iPeak,1}= CurrentAnimal;
                            temp2.drug{iPeak,1}=drugId;
                            temp2.Cond{iPeak,1} = CurrentCondition;
                            temp2.depth{iPeak,1} = CurrentDepth(1,1);
                            temp2.pixelsize{iPeak,1} = CellScans(iScan,itrial).rawImg.metadata.pixelSize;
                            temp2.nFluoPix{iPeak,1} = nFluoPix;
                            temp2.nActivePix{iPeak,1} = nActivePix;
                            
                            % get the indices  and area for a particular ROI
                            if iScan==1 || iScan==3
                                jROIname = temp.roiName{iPeak};
                                ROIindex= strcmp(CellScans(iScan,itrial).calcFindROIs.data.roiNames,jROIname);
                                temp2.area{iPeak,1} = CellScans(iScan,itrial).calcFindROIs.data.area(ROIindex);
                                
                                jx = CellScans(iScan,itrial).calcFindROIs.data.centroidX(ROIindex,1);
                                jy = CellScans(iScan,itrial).calcFindROIs.data.centroidY(ROIindex,1);
                                xclose = round(x_pix*0.03); % 3 percent of pixels
                                yclose = round(y_pix*0.03); % 3 percent of pixels
                                
                                if isempty(jx)
                                    temp2.overlap{iPeak,1} = 0;
                                else
                                    % find FLIKA ROIs and hand selected ROIS that have
                                    % similar centroids
                                    if iScan==1
                                        compScan=2;
                                    elseif iScan==3;
                                        compScan=4;
                                    end
                                    
                                    for kROI= 1:size(CellScans(compScan,1).calcFindROIs.data.roiNames,1)
                                        % get the centroids of the handclicked ROIs
                                        kx = CellScans(compScan,itrial).calcFindROIs.data.centroidX(kROI,1);
                                        ky = CellScans(compScan,itrial).calcFindROIs.data.centroidY(kROI,1);
                                        
                                        % check if they are too close (actually same region)
                                        if jx<=(kx+xclose) && jx>=(kx-xclose) && jy<=(ky+yclose) && jy>=(ky-yclose) % only %3 of pixels apart
                                            spatialcorr(kROI)  = 1;
                                        end
                                    end
                                    
                                    if exist('spatialcorr','var')
                                        indx = find(spatialcorr>0);
                                        temp2.overlap{iPeak,1}= CellScans(compScan,1).calcFindROIs.data.roiNames{indx};
                                    else
                                        temp2.overlap{iPeak,1} = 0;
                                    end
                                    
                                    clear spatialcorr
                                end
                            else
                                temp2.area{iPeak,1} = 0;
                                temp2.overlap{iPeak,1} = 0;
                            end
                            
                        end
                        
                        
                        
                        %%
                        isFirst = (itrial == 1 && iScan == 1);
                        if isFirst
                            data.Trial = {};
                            data.Animal = {};
                            data.Channel = {};
                            data.Spot = {};
                            data.Drug={};
                            data.Condition = {};
                            data.Depth = {};
                            data.area = {};
                            data.overlap = {};
                            data.pixelsize={};
                            data.nFluoPix = {};
                            data.nActivePix = {};
                        end
                        data.Trial= [data.Trial; temp2.trialname];
                        data.Animal= [data.Animal; temp2.animalname];
                        data.Channel= [data.Channel; temp2.channel];
                        data.Spot= [data.Spot; temp2.Spot];
                        data.Drug= [data.Drug; temp2.drug];
                        data.Condition= [data.Condition; temp2.Cond];
                        data.Depth= [data.Depth; temp2.depth];
                        data.area= [data.area; temp2.area];
                        data.overlap= [data.overlap; temp2.overlap];
                        data.pixelsize= [data.pixelsize; temp2.pixelsize];
                        data.nFluoPix = [data.nFluoPix; temp2.nFluoPix];
                        data.nActivePix = [data.nActivePix; temp2.nActivePix];
                        
                        clearvars temp temp2
                        
                        
                        % make a table of trace info
                        
                        %traces output processes
                        if strcmp(CellScans(iScan,itrial).calcFindROIs.data.roiNames{1,1}, 'none')
                            continue
                        else
                            traces= CellScans(iScan,itrial).calcMeasureROIs.data.tracesNorm;
                            %preallocate
                            Trace_data=cell(size(traces,2),10);
                            for iROI = 1:size(traces,2)
                                Trace_data{iROI,1}= CellScans(iScan,itrial).calcFindROIs.data.roiNames{iROI,1};
                                Trace_data{iROI,2}= strcat('trial', num2str(itrial,'%02d'));
                                if iScan==1 || iScan==2
                                    Trace_data{iROI,3}= 'GCaMP';
                                else
                                    Trace_data{iROI,3}= 'RCaMP';
                                end
                                
                                
                                Trace_data{iROI,4}= spotId;
                                Trace_data{iROI,5}= CurrentAnimal;
                                Trace_data{iROI,6}= CurrentCondition;
                                Trace_data{iROI,7}= drugId;
                                Trace_data{iROI,8} = CurrentDepth(1,1);
                                Trace_data{iROI,9} = CurrentBaseline(1,1);
                                Trace_data{iROI,10} = traces(:,iROI);
                                if iScan==1 || iScan==3
                                    Trace_data{iROI,11} = CellScans(iScan,itrial).calcFindROIs.data.roiIdxs{iROI,1};
                                    Trace_data{iROI,12} = CellScans(iScan,itrial).rawImg.metadata.pixelSize;
                                else
                                    
                                    Trace_data{iROI,11} = CellScans(iScan,itrial).calcFindROIs.data.roiMask(:,:,iROI);
                                    Trace_data{iROI,12} = CellScans(iScan,itrial).rawImg.metadata.pixelSize;
                                end
                                
                                % get the indices  and area for a particular ROI
                                if iScan==1 || iScan==3
                                    jx = CellScans(iScan,itrial).calcFindROIs.data.centroidX(iROI,1);
                                    jy = CellScans(iScan,itrial).calcFindROIs.data.centroidY(iROI,1);
                                    xclose = round(x_pix*0.03); % 3 percent of pixels
                                    yclose = round(y_pix*0.03); % 3 percent of pixels
                                    
                                    if isempty(jx)
                                        Trace_data{iROI,13} = 0;
                                    else
                                        % find FLIKA ROIs and hand selected ROIS that have
                                        % similar centroids
                                        if iScan==1
                                            compScan=2;
                                        elseif iScan==3;
                                            compScan=4;
                                        end
                                        for kROI= 1:length(CellScans(compScan,1).calcFindROIs.data.roiNames)
                                            % get the centroids of the handclicked ROIs
                                            kx = CellScans(compScan,itrial).calcFindROIs.data.centroidX(kROI,1);
                                            ky = CellScans(compScan,itrial).calcFindROIs.data.centroidY(kROI,1);
                                            
                                            % check if they are too close (actually same region)
                                            if jx<=(kx+xclose) && jx>=(kx-xclose) && jy<=(ky+yclose) && jy>=(ky-yclose) % only %3 of pixels apart
                                                spatialcorr(kROI)  = 1;
                                            end
                                        end
                                        
                                        
                                        if exist('spatialcorr','var')
                                            indx = find(spatialcorr>0);
                                            Trace_data{iROI,13}= CellScans(compScan,1).calcFindROIs.data.roiNames{indx};
                                        else
                                            Trace_data{iROI,13} = 0;
                                        end
                                        clear spatialcorr
                                    end
                                else
                                    Trace_data{iROI,13} = 0;
                                end
                             FrameRate= CellScans(1, 1).rawImg.metadata.frameRate;
                            Trace_data{iROI,14} = FrameRate; % frameRate
                            
                            nFrames=length(traces(:,iROI));
                            TimeX(1:nFrames) = (1:nFrames)/FrameRate;
                            
                            % Calculate the first peak onset time and AUC after stim
                            BL_time=round(BL_frames/FrameRate);  % number of s for baseline
                            baselineCorrectedTime=TimeX-BL_time;
                            
                            % onset time  % 2.5SD from baseline and
                            % smoothing trace at 11 points (5 each side
                            % of middle)
                            Onsets=find_first_onset_time(baselineCorrectedTime(10:end), traces(10:end,iROI),2.5,2);
                            if isempty(Onsets)
                                Onsets=nan(1,1);
                            end
                            Trace_data{iROI,15}= Onsets;
                            
                            % trace AUC
                            x2=round(FrameRate*(BL_time+1));
                            x3= round(FrameRate*(BL_time+10));
                            Trace_data{iROI,16}=trapz(traces(BL_frames:x2,iROI));
                            Trace_data{iROI,17}=trapz(traces(BL_frames:x3,iROI));                     
                            
                            end
                        end
                            
                        All_traces=vertcat(All_traces, Trace_data);
                        clearvars Trace_data
                    end
                end
                
                
                dataNames=fieldnames(data);
                data2= struct2cell(data);
                data3= [data2{:}];
                
                AllData=vertcat(AllData, data3);
                
                
                clearvars data data3
            end
        end
    end
end


% %% Save all data for R analysis
AllData2= [dataNames';AllData];

%onsetTimeTable
names={'ROI','Trial','Channel','Spot','Animal', 'Condition','Drug','depth','baseline',...
    'trace','ROIIdx','PixelSize','overlap','FrameRate','OnsetTime','TraceAUC1','TraceAUC10'};
All_traces2=vertcat(names, All_traces);
All_traces2(:,10)=[];
All_traces2(:,10)=[];

cd(fullfile(Settings.MainDir, 'Results'));
% write date to created file
cell2csv(SaveFiles{1,1}, AllData2);
save(SaveFiles{1,3}, 'All_traces','-v7.3');
save(SaveFiles{1,2}, 'AllData2','-v7.3');
cell2csv(SaveFiles{1,4}, All_traces2);



