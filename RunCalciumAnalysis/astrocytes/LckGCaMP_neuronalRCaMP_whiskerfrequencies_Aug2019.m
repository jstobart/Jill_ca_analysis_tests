% new params for extracitng data
% shorter baseline (slightly) to account for changes in fluor and

close all; clear variables;

All_traces=[]; AllData=[]; All_traces2=[]; AllData2=[]; LckData=[]; LckFieldData=[];
%% Information about your images

% folder where data should be saved for each animal
Settings.ResultsFolder = 'G:\Data\GCaMP_RCaMP\WhiskerFrequencies\Results';

% File names for saving
SaveFiles{1,1} = fullfile(Settings.ResultsFolder,'FilesforMatlab', 'Alice_peaks_frequencies.mat');
SaveFiles{1,2} = fullfile(Settings.ResultsFolder,'FilesforR','Alice_peaks_frequencies.csv');
SaveFiles{1,3}= fullfile(Settings.ResultsFolder,'FilesforMatlab','Alice_traces_frequencies.mat');
SaveFiles{1,4}= fullfile(Settings.ResultsFolder,'FilesforR','Alice_onset_time_frequencies.csv');
SaveFiles{1,5}= fullfile(Settings.ResultsFolder,'FilesforMatlab','Alice_Lck_field_frequencies.mat');
SaveFiles{1,6}= fullfile(Settings.ResultsFolder,'FilesforR','Alice_Lck_field_frequencies.csv');

Settings.FileNames = {  % folder names where images and roiSet.zip is found
    'G:\Data\GCaMP_RCaMP\WhiskerFrequencies\Alice\2019_06_14_frequencies\spot1',...
    'G:\Data\GCaMP_RCaMP\WhiskerFrequencies\Alice\2019_06_14_frequencies\spot2',...
    'G:\Data\GCaMP_RCaMP\WhiskerFrequencies\Alice\2019_06_22_frequencies\spot2',...
    'G:\Data\GCaMP_RCaMP\WhiskerFrequencies\Alice\2019_07_02_frequencies\spot1',...
    'G:\Data\GCaMP_RCaMP\WhiskerFrequencies\Alice\2019_07_02_frequencies\spot2',...
    'G:\Data\GCaMP_RCaMP\WhiskerFrequencies\Alice\2019_07_03_frequencies\spot1',...
    'G:\Data\GCaMP_RCaMP\WhiskerFrequencies\Alice\2019_07_03_frequencies\spot2',...
    };

Settings.Baseline = 1.5; % time (s) before the whisker stimulator starts,  2s for short trials (10s long), 5s for long trials
Settings.Animal = 'Alice';

channel = struct('Ca_Neuron',1,'Ca_Memb_Astro',2);  % what's on each channel


%% Loop through each file and make Cell Scans for LckGCaMP and RCaMP

for iSpot= 1:length(Settings.FileNames)
    
    % loop through each spot
    SpotRoot= Settings.FileNames{iSpot};
    SpotId1=SpotRoot(end-27:end-17);
    SpotId2=SpotRoot(end-4:end);
    SpotId=strcat(SpotId1,SpotId2);
    
    % Get a list of all files and folders in this folder.
    % different stimulation conditions
    ConditionFolders = dir(SpotRoot);  % look for folders
    ConditionFolders(ismember( {ConditionFolders.name}, {'.', '..'})) = [];
    % Get a logical vector that tells which is a directory.
    names    = {ConditionFolders.name};
    dirFlags = [ConditionFolders.isdir];
    % Extract only those that are directories.
    subDirsNames = names(dirFlags);
    
    % find folders for each condition (stim vs nostim, different
    % frequencies)
    for iCondition= 1:length(subDirsNames)
        FolderName=fullfile(SpotRoot, subDirsNames(iCondition));
        
        % Get a list of all T-series folders in this folder.
        % different trials
        TrialFolders = dir(FolderName{1,1});  % look for folders
        TrialFolders(ismember( {TrialFolders.name}, {'.', '..'})) = [];
        % Get a logical vector that tells which is a directory.
        names2    = {TrialFolders.name};
        dirFlags2 = [TrialFolders.isdir];
        % Extract only those that are directories.
        TrialNames = names2(dirFlags2);
        
        % load each trial
        for iTrial= 1:length(TrialNames)
            % Get the T-series path from the folder
            testRoot =fullfile(FolderName, TrialNames{iTrial});
            
            expfiles = dir(fullfile(testRoot{1,1},'*.xml'));  % look for the xml file
            fnList = fullfile(testRoot{1,1}, expfiles.name);
            
            % Load data with BioFormats
            ImgArray =  BioFormats(fnList, channel);
            
            %% Motion Correction
            
            refImg=mean(ImgArray.rawdata(:,:,1,1:3),4);
            
            ImgArray = ImgArray.motion_correct('refImg',refImg, 'ch',1,'maxShift', 10,'minCorr', 0.4);%,'doPlot',true);
            
            BL_frames= floor(Settings.Baseline*ImgArray.metadata.frameRate); % number of baseline frames
            
            %% Configs for Lck GCaMP Cell Scan
            %ASTROCYTES
            % 2D automated selection for peaks
            AC_findConf = ConfigFindROIsFLIKA_2D.from_preset('ca_memb_astro', 'baselineFrames',...
                3:BL_frames,'freqPassBand',0.5,'sigmaXY', 2,...
                'sigmaT', 0.14,'thresholdPuff', 7, 'threshold2D', 0.1,...
                'minRiseTime',0.14, 'maxRiseTime', 1,'minROIArea', 10,...
                'dilateXY', 5, 'dilateT', 0.3,'erodeXY', 1, 'erodeT', 0.1);
            
            % measure ROIs (extract the traces)
            AC_measureConf = ConfigMeasureROIsDummy('baselineFrames', 3:BL_frames);
            
            % filter the traces to detect the peaks and get info about them
            AC_detectConf = ConfigDetectSigsClsfy('baselineFrames', 3:BL_frames,...
                'propagateNaNs', false, 'excludeNaNs', false, 'lpWindowTime', 1.5, 'spFilterOrder', 2,...
                'spPassBandMin',0.05, 'spPassBandMax', 0.5, 'thresholdLP', 7,'thresholdSP', 7);
            
            % Combine the configs into a CellScan config for membrane tagged GCaMP
            AC_configCS= ConfigCellScan(AC_findConf, AC_measureConf, AC_detectConf); %
            
            %% Configs for neuronal RCaMP Cell Scan
            % NEURONS
            % hand selected- peaks from cellular structures
            x_pix= size(ImgArray(1,1).rawdata,2); y_pix= size(ImgArray(1,1).rawdata,1);
            scaleF = 1;
            
            zipfiles = dir(fullfile(SpotRoot,'*.zip'));  % look for the zip ROI file
            fnTempList2 = {zipfiles(:).name};
            zipPath = fullfile(SpotRoot, fnTempList2);
            
            % load image J ROIs for mask
            N_findConf{1} = ConfigFindROIsDummy.from_ImageJ(zipPath{1,1}, x_pix, y_pix, scaleF);
            
            % 2D FLIKA selected for peaks from "dendrites"
            N_findConf{2} = ConfigFindROIsFLIKA_2D.from_preset('ca_neuron', 'baselineFrames',...
                BL_frames,'freqPassBand',1,'sigmaXY', 2,...
                'sigmaT', 0.14,'thresholdPuff', 7, 'threshold2D', 0.1,...
                'minRiseTime',0.07, 'maxRiseTime', 1,'minROIArea', 10,...
                'dilateXY', 5, 'dilateT', 0.2,'erodeXY', 1, 'erodeT', 0.1,...
                'discardBorderROIs',false);
            
            % measure ROIs (extract the traces)
            N_measureConf = ConfigMeasureROIsDummy();
            
            % filter the traces to detect the peaks and get info about them
            N_detectConf = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...
                'propagateNaNs', false,'excludeNaNs', false, 'lpWindowTime', 2, 'spFilterOrder', 2,...
                'spPassBandMin',0.1, 'spPassBandMax', 1, 'thresholdLP', 6,'thresholdSP', 5);
            
            
            % Combine the configs into a CellScan config for neuronal RCaMP
            N_configCS_ImageJ= ConfigCellScan(N_findConf{1}, N_measureConf, N_detectConf); %
            N_configCS_FLIKA= ConfigCellScan(N_findConf{2}, N_measureConf, N_detectConf); %
            
            
            %% Create CellScan objects
            astrocytes = CellScan(fnList, ImgArray, AC_configCS, 2); % peaks from automated ROIs
            neurons1 = CellScan(fnList, ImgArray, N_configCS_ImageJ, 1); % peaks from hand clicked ROIs
            neurons2 = CellScan(fnList, ImgArray, N_configCS_FLIKA, 1); % peaks from automated ROIs
            
            % Process the images
            astrocytes =astrocytes.process();
            neurons1 =neurons1.process();
            neurons2 =neurons2.process();
            
            % Make the debugging plots
            %             astrocytes.plot();
            %             neurons1.plot();
            %             neurons2.plot();
            
            %astrocytes.opt_config()
            %neurons2.opt_config()
            %neurons2.plot('signals');
            
            %% Extract/Calculate data we want and create table to output data
            %             if bad
            %                 continue
            %             else
            % make a giant data table for the signal peaks
            listFields = {'amplitude', 'fullWidth', 'halfWidth', ...
                'numPeaks', 'peakTime', 'peakStart', 'peakStartHalf', ...
                'peakType', 'prominence', 'roiName', 'peakAUC'};
            
            CellScans=vertcat(astrocytes, neurons1, neurons2);
            
            
            % loop through cellscans and pull out the data
            for iScan=1:size(CellScans,1)
                
                % peak output
                temp=CellScans(iScan).calcDetectSigs.data;
                temp2.trialname ={};
                temp2.animalname = {};
                temp2.channel = {};
                temp2.Spot = {};
                temp2.Cond = {};
                temp2.area = {};
                temp2.pixelsize = {};
                temp2.overlap ={};
                
                % extract fields from Class
                for jField = 1:numel(listFields)
                    isFirst = (iTrial == 1 && iScan == 1);
                    if isFirst
                        data.(listFields{jField}) = {};
                    end
                    data.(listFields{jField}) = [data.(listFields{jField}); ...
                        temp.(listFields{jField})];
                end
                
                % create fields for trial, animal, spot, condition, etc.
                for iPeak = 1:length(temp.amplitude)
                    temp2.trialname{iPeak,1}=strcat('trial', num2str(iTrial,'%02d'));
                    if iScan==1
                        temp2.channel{iPeak,1}= 'GCaMP';
                    else
                        temp2.channel{iPeak,1}= 'RCaMP';
                    end
                    temp2.Spot{iPeak,1}= SpotId;
                    temp2.animalname{iPeak,1}= Settings.Animal;
                    temp2.Cond{iPeak,1} = subDirsNames{1,iCondition};
                    temp2.pixelsize{iPeak,1} = CellScans(iScan).rawImg.metadata.pixelSize;
                    
                    % get the indices  and area for a particular ROI
                    if iScan==3
                        jROIname = temp.roiName{iPeak};
                        ROIindex= strcmp(CellScans(iScan,1).calcFindROIs.data.roiNames,jROIname);
                        temp2.area{iPeak,1} = CellScans(iScan,1).calcFindROIs.data.area(ROIindex);
                        
                        jx = CellScans(iScan,1).calcFindROIs.data.centroidX(ROIindex,1);
                        jy = CellScans(iScan,1).calcFindROIs.data.centroidY(ROIindex,1);
                        xclose = round(x_pix*0.03); % 3 percent of pixels
                        yclose = round(y_pix*0.03); % 3 percent of pixels
                        
                        if isempty(jx)
                            temp2.overlap{iPeak,1} = 0;
                        else
                            % find NEURONAL FLIKA ROIs and hand selected ROIS that have
                            % similar centroids
                            
                            for kROI= 1:size(CellScans(2,1).calcFindROIs.data.roiNames,1)
                                % get the centroids of the handclicked ROIs
                                kx = CellScans(2,1).calcFindROIs.data.centroidX(kROI,1);
                                ky = CellScans(2,1).calcFindROIs.data.centroidY(kROI,1);
                                
                                % check if they are too close (actually same region)
                                if jx<=(kx+xclose) && jx>=(kx-xclose) && jy<=(ky+yclose) && jy>=(ky-yclose) % only %3 of pixels apart
                                    spatialcorr(kROI)  = 1;
                                end
                            end
                            
                            if exist('spatialcorr','var')
                                indx = find(spatialcorr>0);
                                temp2.overlap{iPeak,1}= CellScans(2,1).calcFindROIs.data.roiNames{indx};
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
                isFirst = (iTrial == 1 && iScan == 1);
                if isFirst
                    data.Trial = {};
                    data.Animal = {};
                    data.Channel = {};
                    data.Spot = {};
                    data.Condition = {};
                    data.area = {};
                    data.pixelsize={};
                    data.overlap ={};
                end
                data.Trial= [data.Trial; temp2.trialname];
                data.Animal= [data.Animal; temp2.animalname];
                data.Channel= [data.Channel; temp2.channel];
                data.Spot= [data.Spot; temp2.Spot];
                data.Condition= [data.Condition; temp2.Cond];
                data.area= [data.area; temp2.area];
                data.pixelsize= [data.pixelsize; temp2.pixelsize];
                data.overlap= [data.overlap; temp2.overlap];
                
                clearvars temp temp2
                
                
                %% make a table of ROI info and traces
                
                %traces output processes
                if strcmp(CellScans(iScan).calcFindROIs.data.roiNames{1,1}, 'none')
                    continue
                else
                    traces= CellScans(iScan).calcMeasureROIs.data.tracesNorm;
                    
                    %preallocate
                    Trace_data=cell(size(traces,2),10);
                    mask = zeros(x_pix, y_pix);
                    
                    for iROI = 1:size(traces,2)
                        Trace_data{iROI,1}= CellScans(iScan).calcFindROIs.data.roiNames{iROI,1};
                        Trace_data{iROI,2}= strcat('trial', num2str(iTrial,'%02d'));
                        if iScan==1
                            Trace_data{iROI,3}= 'GCaMP';
                        else
                            Trace_data{iROI,3}= 'RCaMP';
                        end
                        
                        Trace_data{iROI,4}= SpotId;
                        Trace_data{iROI,5}= Settings.Animal;
                        Trace_data{iROI,6}= subDirsNames{1,iCondition};
                        Trace_data{iROI,7} = Settings.Baseline;
                        Trace_data{iROI,8} = traces(:,iROI);
                        if iScan==1 || iScan==3
                            Trace_data{iROI,9} = CellScans(iScan).calcFindROIs.data.roiIdxs{iROI,1};
                            Trace_data{iROI,10} = CellScans(iScan).rawImg.metadata.pixelSize;
                        else
                            
                            Trace_data{iROI,9} = CellScans(iScan).calcFindROIs.data.roiMask(:,:,iROI);
                            Trace_data{iROI,10} = CellScans(iScan).rawImg.metadata.pixelSize;
                        end
                        
                        FrameRate= CellScans(1).rawImg.metadata.frameRate;
                        Trace_data{iROI,11} = FrameRate; % frameRate
                        
                        nFrames=length(traces(:,iROI));
                        TimeX(1:nFrames) = (1:nFrames)/FrameRate;
                        
                        % Calculate the first peak onset time and AUC after stim
                        BL_time=Settings.Baseline;  % number of s for baseline
                        baselineCorrectedTime=TimeX-BL_time;
                        
                        % onset time  % 2.5SD from baseline and
                        % smoothing trace at 11 points (5 each side
                        % of middle)
                        Onsets=find_first_onset_time(baselineCorrectedTime(10:end), traces(10:end,iROI),2.5,2);
                        if isempty(Onsets)
                            Onsets=nan(1,1);
                        end
                        Trace_data{iROI,12}= Onsets;
                        
                        % trace AUC in the 8 s follow stimulation
                        x2=round(FrameRate*(BL_time+8));
                        Trace_data{iROI,13}=trapz(traces(BL_frames:x2,iROI));
                        
                        % create a ROI mask of all the astrocyte ROIs
                        if iScan==1
                            roiIdx = CellScans(1).calcFindROIs.data.roiIdxs{iROI,1};
                            mask(roiIdx) = true;
                        end
                        
                        % get the indices  and area for a particular ROI
                        if iScan==3
                            jx = CellScans(iScan,1).calcFindROIs.data.centroidX(iROI,1);
                            jy = CellScans(iScan,1).calcFindROIs.data.centroidY(iROI,1);
                            xclose = round(x_pix*0.03); % 3 percent of pixels
                            yclose = round(y_pix*0.03); % 3 percent of pixels
                            
                            if isempty(jx)
                                Trace_data{iROI,14} = 0;
                            else
                                % find FLIKA ROIs and hand selected ROIS that have
                                % similar centroids
                                for kROI= 1:length(CellScans(2,1).calcFindROIs.data.roiNames)
                                    % get the centroids of the handclicked ROIs
                                    kx = CellScans(2,1).calcFindROIs.data.centroidX(kROI,1);
                                    ky = CellScans(2,1).calcFindROIs.data.centroidY(kROI,1);
                                    
                                    % check if they are too close (actually same region)
                                    if jx<=(kx+xclose) && jx>=(kx-xclose) && jy<=(ky+yclose) && jy>=(ky-yclose) % only %3 of pixels apart
                                        spatialcorr(kROI)  = 1;
                                    end
                                end
                                
                                
                                if exist('spatialcorr','var')
                                    indx = find(spatialcorr>0);
                                    Trace_data{iROI,14}= CellScans(2,1).calcFindROIs.data.roiNames{indx};
                                else
                                    Trace_data{iROI,14} = 0;
                                end
                                clear spatialcorr
                            end
                        else
                            Trace_data{iROI,14} = 0;
                        end
                    end
                    
                end
                
                All_traces=vertcat(All_traces, Trace_data);
                clearvars Trace_data
                
                %% Lck field of view specific data
                if iScan==1
                    Lck.trialname{iTrial,1} =strcat('trial', num2str(iTrial,'%02d'));
                    Lck.Spot{iTrial,1}= SpotId;
                    Lck.animalname{iTrial,1}= Settings.Animal;
                    Lck.Cond{iTrial,1} = subDirsNames{1,iCondition};
                    Lck.pixelsize{iTrial,1} = CellScans(iScan).rawImg.metadata.pixelSize;
                    
                    % get fraction of active MD area
                    neuroMask = CellScans(2).calcFindROIs.data.roiMask;
                    if ndims(neuroMask) == 3
                        neuroMask = max(neuroMask, [], 3);
                    end
                    [nFluoPix, nActivePix, nTotalPix] = ...
                        getFracActive(CellScans(1), 'nanMask', ...
                        neuroMask);
                    Lck.nFluoPix{iTrial,1} = nFluoPix;
                    Lck.nActivePix{iTrial,1} = nActivePix;
                    Lck.nTotalPix{iTrial,1} = nTotalPix;
                    
                    % ROI masks output
                    Lck.Trial_ROIMask{iTrial,1} = mask;
                    
                    % all astrocyte ROI masks together
                    trialmasks(:,:,iTrial) = mask;
                    
                end
            end
            
            %             end
        end
        
        % append peak data
        dataNames=fieldnames(data);
        data2= struct2cell(data);
        data3= [data2{:}];
        
        AllData=vertcat(AllData, data3);
        
        clearvars data data2 data3
        
        % Mask of Active Pixels (proportional to number of trials)
        sumImg =sum(trialmasks,3);
        fracImg = sum(trialmasks,3)./numel(TrialNames);
        
        % calculate the score for responses (normalized to the "threshold"
        thresh = 1/numel(TrialNames); % or 1/numel(trials)
        activePxIdx = fracImg > 0;
        activePx = sum(activePxIdx(:));
        fracImg(fracImg <= thresh) = NaN;
        score = nansum(nansum(fracImg)) / activePx;
        
        % save the mask of active pixels
        FigFileName1 = fullfile(FolderName,'Sum_Pixel_Mask.tif');
        FigFileName2 = fullfile(FolderName,'FracActive_Pixel_Mask.tif');
        
        f = figure('visible', 'off');
        imagesc(sumImg);
        axis square
        caxis([0 1])
        colormap('jet')
        colorbar
        saveas(gcf,FigFileName1{1,1})
        close(f)
        
        
        f = figure('visible', 'off');
        imagesc(fracImg);
        axis square
        caxis([0 1])
        colormap('jet')
        colorbar
        saveas(gcf,FigFileName2{1,1})
        close(f)
        
        % store Lck field of view information
        
        for iLck= 1:length(Lck.trialname)
            Lck.Total_ROIMask{iLck,1} = fracImg;
            Lck.Response_Score{iLck,1}=score;
        end
        
        LckNames=fieldnames(Lck);
        Lckdata2= struct2cell(Lck);
        Lckdata3= [Lckdata2{:}];
        
        LckData=vertcat(LckData, Lckdata3);
        
        clearvars Lck Lckdata2 Lckdata3
    end
    
end


% %% Save all data for R analysis
% peak data
AllData2= [dataNames';AllData];

%Field of View data
LckFieldData= [LckNames';LckData];
LckFieldData2=LckFieldData;
LckFieldData(:,9:10)=[];

%onsetTimeTable (ROI table)
names={'ROI','Trial','Channel','Spot','Animal', 'Condition','baseline',...
    'trace','ROIIdx','PixelSize','FrameRate','OnsetTime','TraceAUC1','TraceAUC10'};
All_traces2=vertcat(names, All_traces);
All_traces2(:,8:9)=[];



cd(fullfile(Settings.ResultsFolder));

% write date to created file
save(SaveFiles{1,1}, 'AllData2','-v7.3');
cell2csv(SaveFiles{1,2}, AllData2);
save(SaveFiles{1,3}, 'All_traces','-v7.3');
cell2csv(SaveFiles{1,4}, All_traces2);
save(SaveFiles{1,5}, 'LckFieldData2','-v7.3');
cell2csv(SaveFiles{1,6}, LckFieldData);

close all; clear variables;

%%
All_traces=[]; AllData=[]; All_traces2=[]; AllData2=[]; LckData=[]; LckFieldData=[];
%% Information about your images

% folder where data should be saved for each animal
% folder where data should be saved for each animal
Settings.ResultsFolder = 'G:\Data\GCaMP_RCaMP\WhiskerFrequencies\Results';

% File names for saving
SaveFiles{1,1} = fullfile(Settings.ResultsFolder,'FilesforMatlab', 'Crazy8_peaks_frequencies.mat');
SaveFiles{1,2} = fullfile(Settings.ResultsFolder,'FilesforR','Crazy8_peaks_frequencies.csv');
SaveFiles{1,3}= fullfile(Settings.ResultsFolder,'FilesforMatlab','Crazy8_traces_frequencies.mat');
SaveFiles{1,4}= fullfile(Settings.ResultsFolder,'FilesforR','Crazy8_onset_time_frequencies.csv');
SaveFiles{1,5}= fullfile(Settings.ResultsFolder,'FilesforMatlab','Crazy8_Lck_field_frequencies.mat');
SaveFiles{1,6}= fullfile(Settings.ResultsFolder,'FilesforR','Crazy8_Lck_field_frequencies.csv');

Settings.FileNames = {  % folder names where images and roiSet.zip is found
    'G:\Data\GCaMP_RCaMP\WhiskerFrequencies\Crazy8\2019_06_14_frequencies\spot1',...
    'G:\Data\GCaMP_RCaMP\WhiskerFrequencies\Crazy8\2019_06_14_frequencies\spot2',...
    'G:\Data\GCaMP_RCaMP\WhiskerFrequencies\Crazy8\2019_06_22_frequencies\spot1',...
    'G:\Data\GCaMP_RCaMP\WhiskerFrequencies\Crazy8\2019_06_22_frequencies\spot2',...
    'G:\Data\GCaMP_RCaMP\WhiskerFrequencies\Crazy8\2019_07_02_frequencies\spot1',...
    'G:\Data\GCaMP_RCaMP\WhiskerFrequencies\Crazy8\2019_07_02_frequencies\spot2',...
    'G:\Data\GCaMP_RCaMP\WhiskerFrequencies\Crazy8\2019_07_03_frequencies\spot1',...
    'G:\Data\GCaMP_RCaMP\WhiskerFrequencies\Crazy8\2019_07_03_frequencies\spot2',...
    };
Settings.Baseline = 2; % time (s) before the whisker stimulator starts,  2s for short trials (10s long), 5s for long trials
Settings.Animal = 'Crazy8';

channel = struct('Ca_Neuron',1,'Ca_Memb_Astro',2);  % what's on each channel


%% Loop through each file and make Cell Scans for LckGCaMP and RCaMP

for iSpot= 1:length(Settings.FileNames)
    
    % loop through each spot
    SpotRoot= Settings.FileNames{iSpot};
    SpotId=SpotRoot(end-15:end);
    SpotId(11)='_';
    
    % Get a list of all files and folders in this folder.
    % different stimulation conditions
    ConditionFolders = dir(SpotRoot);  % look for folders
    ConditionFolders(ismember( {ConditionFolders.name}, {'.', '..'})) = [];
    % Get a logical vector that tells which is a directory.
    names    = {ConditionFolders.name};
    dirFlags = [ConditionFolders.isdir];
    % Extract only those that are directories.
    subDirsNames = names(dirFlags);
    
    % find folders for each condition (stim vs nostim, different
    % frequencies)
    for iCondition= 1:length(subDirsNames)
        FolderName=fullfile(SpotRoot, subDirsNames(iCondition));
        
        % Get a list of all T-series folders in this folder.
        % different trials
        TrialFolders = dir(FolderName{1,1});  % look for folders
        TrialFolders(ismember( {TrialFolders.name}, {'.', '..'})) = [];
        % Get a logical vector that tells which is a directory.
        names2    = {TrialFolders.name};
        dirFlags2 = [TrialFolders.isdir];
        % Extract only those that are directories.
        TrialNames = names2(dirFlags2);
        
        % load each trial
        for iTrial= 1:length(TrialNames)
            % Get the T-series path from the folder
            testRoot =fullfile(FolderName, TrialNames{iTrial});
            
            expfiles = dir(fullfile(testRoot{1,1},'*.xml'));  % look for the xml file
            fnList = fullfile(testRoot{1,1}, expfiles.name);
            
            % Load data with BioFormats
            ImgArray =  BioFormats(fnList, channel);
            
            %% Motion Correction
            
            refImg=mean(ImgArray.rawdata(:,:,1,1:3),4);
            
            ImgArray = ImgArray.motion_correct('refImg',refImg, 'ch',1,'maxShift', 10,'minCorr', 0.4);%,'doPlot',true);
            
            BL_frames= floor(Settings.Baseline*ImgArray.metadata.frameRate); % number of baseline frames
            
            %% Configs for Lck GCaMP Cell Scan
            %ASTROCYTES
            % 2D automated selection for peaks
            AC_findConf = ConfigFindROIsFLIKA_2D.from_preset('ca_memb_astro', 'baselineFrames',...
                BL_frames,'freqPassBand',0.5,'sigmaXY', 2,...
                'sigmaT', 0.14,'thresholdPuff', 7, 'threshold2D', 0.1,...
                'minRiseTime',0.14, 'maxRiseTime', 1,'minROIArea', 10,...
                'dilateXY', 5, 'dilateT', 0.3,'erodeXY', 1, 'erodeT', 0.1);
            
            % measure ROIs (extract the traces)
            AC_measureConf = ConfigMeasureROIsDummy();
            
            % filter the traces to detect the peaks and get info about them
            AC_detectConf = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...
                'propagateNaNs', false, 'excludeNaNs', false, 'lpWindowTime', 1.5, 'spFilterOrder', 2,...
                'spPassBandMin',0.05, 'spPassBandMax', 0.5, 'thresholdLP', 7,'thresholdSP', 7);
            
            % Combine the configs into a CellScan config for membrane tagged GCaMP
            AC_configCS= ConfigCellScan(AC_findConf, AC_measureConf, AC_detectConf); %
            
            %% Configs for neuronal RCaMP Cell Scan
            % NEURONS
            % hand selected- peaks from cellular structures
            x_pix= size(ImgArray(1,1).rawdata,2); y_pix= size(ImgArray(1,1).rawdata,1);
            scaleF = 1;
            
            zipfiles = dir(fullfile(SpotRoot,'*.zip'));  % look for the zip ROI file
            fnTempList2 = {zipfiles(:).name};
            zipPath = fullfile(SpotRoot, fnTempList2);
            
            % load image J ROIs for mask
            N_findConf{1} = ConfigFindROIsDummy.from_ImageJ(zipPath{1,1}, x_pix, y_pix, scaleF);
            
            % 2D FLIKA selected for peaks from "dendrites"
            N_findConf{2} = ConfigFindROIsFLIKA_2D.from_preset('ca_neuron', 'baselineFrames',...
                BL_frames,'freqPassBand',1,'sigmaXY', 2,...
                'sigmaT', 0.14,'thresholdPuff', 7, 'threshold2D', 0.1,...
                'minRiseTime',0.07, 'maxRiseTime', 1,'minROIArea', 10,...
                'dilateXY', 5, 'dilateT', 0.2,'erodeXY', 1, 'erodeT', 0.1,...
                'discardBorderROIs',false);
            
            % measure ROIs (extract the traces)
            N_measureConf = ConfigMeasureROIsDummy();
            
            % filter the traces to detect the peaks and get info about them
            N_detectConf = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...
                'propagateNaNs', false,'excludeNaNs', false, 'lpWindowTime', 5, 'spFilterOrder', 2,...
                'spPassBandMin',0.1, 'spPassBandMax', 1, 'thresholdLP', 6, 'thresholdSP', 5);
            
            
            % Combine the configs into a CellScan config for neuronal RCaMP
            N_configCS_ImageJ= ConfigCellScan(N_findConf{1}, N_measureConf, N_detectConf); %
            N_configCS_FLIKA= ConfigCellScan(N_findConf{2}, N_measureConf, N_detectConf); %
            
            
            %% Create CellScan objects
            astrocytes = CellScan(fnList, ImgArray, AC_configCS, 2); % peaks from automated ROIs
            neurons1 = CellScan(fnList, ImgArray, N_configCS_ImageJ, 1); % peaks from hand clicked ROIs
            neurons2 = CellScan(fnList, ImgArray, N_configCS_FLIKA, 1); % peaks from automated ROIs
            
            % Process the images
            astrocytes =astrocytes.process();
            neurons1 =neurons1.process();
            neurons2 =neurons2.process();
            
            % Make the debugging plots
            %             astrocytes.plot();
            %             neurons1.plot();
            %             neurons2.plot();
            
            %astrocytes.opt_config()
            %neurons2.opt_config()
            %neurons2.plot('signals');
            
            %% Extract/Calculate data we want and create table to output data
            
            % make a giant data table for the signal peaks
            listFields = {'amplitude', 'fullWidth', 'halfWidth', ...
                'numPeaks', 'peakTime', 'peakStart', 'peakStartHalf', ...
                'peakType', 'prominence', 'roiName', 'peakAUC'};
            
            CellScans=vertcat(astrocytes, neurons1, neurons2);
            
            
            % loop through cellscans and pull out the data
            for iScan=1:size(CellScans,1)
                
                % peak output
                temp=CellScans(iScan).calcDetectSigs.data;
                temp2.trialname ={};
                temp2.animalname = {};
                temp2.channel = {};
                temp2.Spot = {};
                temp2.Cond = {};
                temp2.area = {};
                temp2.pixelsize = {};
                temp2.overlap ={};
                
                % extract fields from Class
                for jField = 1:numel(listFields)
                    isFirst = (iTrial == 1 && iScan == 1);
                    if isFirst
                        data.(listFields{jField}) = {};
                    end
                    data.(listFields{jField}) = [data.(listFields{jField}); ...
                        temp.(listFields{jField})];
                end
                
                % create fields for trial, animal, spot, condition, etc.
                for iPeak = 1:length(temp.amplitude)
                    temp2.trialname{iPeak,1}=strcat('trial', num2str(iTrial,'%02d'));
                    if iScan==1
                        temp2.channel{iPeak,1}= 'GCaMP';
                    else
                        temp2.channel{iPeak,1}= 'RCaMP';
                    end
                    temp2.Spot{iPeak,1}= SpotId;
                    temp2.animalname{iPeak,1}= Settings.Animal;
                    temp2.Cond{iPeak,1} = subDirsNames{1,iCondition};
                    temp2.pixelsize{iPeak,1} = CellScans(iScan).rawImg.metadata.pixelSize;
                    
                    % get the indices  and area for a particular ROI
                    if iScan==3
                        jROIname = temp.roiName{iPeak};
                        ROIindex= strcmp(CellScans(iScan,1).calcFindROIs.data.roiNames,jROIname);
                        temp2.area{iPeak,1} = CellScans(iScan,1).calcFindROIs.data.area(ROIindex);
                        
                        jx = CellScans(iScan,1).calcFindROIs.data.centroidX(ROIindex,1);
                        jy = CellScans(iScan,1).calcFindROIs.data.centroidY(ROIindex,1);
                        xclose = round(x_pix*0.03); % 3 percent of pixels
                        yclose = round(y_pix*0.03); % 3 percent of pixels
                        
                        if isempty(jx)
                            temp2.overlap{iPeak,1} = 0;
                        else
                            % find NEURONAL FLIKA ROIs and hand selected ROIS that have
                            % similar centroids
                            
                            for kROI= 1:size(CellScans(2,1).calcFindROIs.data.roiNames,1)
                                % get the centroids of the handclicked ROIs
                                kx = CellScans(2,1).calcFindROIs.data.centroidX(kROI,1);
                                ky = CellScans(2,1).calcFindROIs.data.centroidY(kROI,1);
                                
                                % check if they are too close (actually same region)
                                if jx<=(kx+xclose) && jx>=(kx-xclose) && jy<=(ky+yclose) && jy>=(ky-yclose) % only %3 of pixels apart
                                    spatialcorr(kROI)  = 1;
                                end
                            end
                            
                            if exist('spatialcorr','var')
                                indx = find(spatialcorr>0);
                                temp2.overlap{iPeak,1}= CellScans(2,1).calcFindROIs.data.roiNames{indx};
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
                isFirst = (iTrial == 1 && iScan == 1);
                if isFirst
                    data.Trial = {};
                    data.Animal = {};
                    data.Channel = {};
                    data.Spot = {};
                    data.Condition = {};
                    data.area = {};
                    data.pixelsize={};
                    data.overlap ={};
                end
                data.Trial= [data.Trial; temp2.trialname];
                data.Animal= [data.Animal; temp2.animalname];
                data.Channel= [data.Channel; temp2.channel];
                data.Spot= [data.Spot; temp2.Spot];
                data.Condition= [data.Condition; temp2.Cond];
                data.area= [data.area; temp2.area];
                data.pixelsize= [data.pixelsize; temp2.pixelsize];
                data.overlap= [data.overlap; temp2.overlap];
                
                clearvars temp temp2
                
                
                %% make a table of ROI info and traces
                
                %traces output processes
                if strcmp(CellScans(iScan).calcFindROIs.data.roiNames{1,1}, 'none')
                    continue
                else
                    traces= CellScans(iScan).calcMeasureROIs.data.tracesNorm;
                    
                    %preallocate
                    Trace_data=cell(size(traces,2),10);
                    mask = zeros(x_pix, y_pix);
                    
                    for iROI = 1:size(traces,2)
                        Trace_data{iROI,1}= CellScans(iScan).calcFindROIs.data.roiNames{iROI,1};
                        Trace_data{iROI,2}= strcat('trial', num2str(iTrial,'%02d'));
                        if iScan==1
                            Trace_data{iROI,3}= 'GCaMP';
                        else
                            Trace_data{iROI,3}= 'RCaMP';
                        end
                        
                        Trace_data{iROI,4}= SpotId;
                        Trace_data{iROI,5}= Settings.Animal;
                        Trace_data{iROI,6}= subDirsNames{1,iCondition};
                        Trace_data{iROI,7} = Settings.Baseline;
                        Trace_data{iROI,8} = traces(:,iROI);
                        if iScan==1 || iScan==3
                            Trace_data{iROI,9} = CellScans(iScan).calcFindROIs.data.roiIdxs{iROI,1};
                            Trace_data{iROI,10} = CellScans(iScan).rawImg.metadata.pixelSize;
                        else
                            
                            Trace_data{iROI,9} = CellScans(iScan).calcFindROIs.data.roiMask(:,:,iROI);
                            Trace_data{iROI,10} = CellScans(iScan).rawImg.metadata.pixelSize;
                        end
                        
                        FrameRate= CellScans(1).rawImg.metadata.frameRate;
                        Trace_data{iROI,11} = FrameRate; % frameRate
                        
                        nFrames=length(traces(:,iROI));
                        TimeX(1:nFrames) = (1:nFrames)/FrameRate;
                        
                        % Calculate the first peak onset time and AUC after stim
                        BL_time=Settings.Baseline;  % number of s for baseline
                        baselineCorrectedTime=TimeX-BL_time;
                        
                        % onset time  % 2.5SD from baseline and
                        % smoothing trace at 11 points (5 each side
                        % of middle)
                        Onsets=find_first_onset_time(baselineCorrectedTime(10:end), traces(10:end,iROI),2.5,2);
                        if isempty(Onsets)
                            Onsets=nan(1,1);
                        end
                        Trace_data{iROI,12}= Onsets;
                        
                        % trace AUC in the 8 s follow stimulation
                        x2=round(FrameRate*(BL_time+8));
                        Trace_data{iROI,13}=trapz(traces(BL_frames:x2,iROI));
                        
                        % create a ROI mask of all the astrocyte ROIs
                        if iScan==1
                            roiIdx = CellScans(1).calcFindROIs.data.roiIdxs{iROI,1};
                            mask(roiIdx) = true;
                        end
                        
                        % get the indices  and area for a particular ROI
                        if iScan==3
                            jx = CellScans(iScan,1).calcFindROIs.data.centroidX(iROI,1);
                            jy = CellScans(iScan,1).calcFindROIs.data.centroidY(iROI,1);
                            xclose = round(x_pix*0.03); % 3 percent of pixels
                            yclose = round(y_pix*0.03); % 3 percent of pixels
                            
                            if isempty(jx)
                                Trace_data{iROI,14} = 0;
                            else
                                % find FLIKA ROIs and hand selected ROIS that have
                                % similar centroids
                                for kROI= 1:length(CellScans(2,1).calcFindROIs.data.roiNames)
                                    % get the centroids of the handclicked ROIs
                                    kx = CellScans(2,1).calcFindROIs.data.centroidX(kROI,1);
                                    ky = CellScans(2,1).calcFindROIs.data.centroidY(kROI,1);
                                    
                                    % check if they are too close (actually same region)
                                    if jx<=(kx+xclose) && jx>=(kx-xclose) && jy<=(ky+yclose) && jy>=(ky-yclose) % only %3 of pixels apart
                                        spatialcorr(kROI)  = 1;
                                    end
                                end
                                
                                
                                if exist('spatialcorr','var')
                                    indx = find(spatialcorr>0);
                                    Trace_data{iROI,14}= CellScans(2,1).calcFindROIs.data.roiNames{indx};
                                else
                                    Trace_data{iROI,14} = 0;
                                end
                                clear spatialcorr
                            end
                        else
                            Trace_data{iROI,14} = 0;
                        end
                    end
                    
                end
                
                All_traces=vertcat(All_traces, Trace_data);
                clearvars Trace_data
                
                %% Lck field of view specific data
                if iScan==1
                    Lck.trialname{iTrial,1} =strcat('trial', num2str(iTrial,'%02d'));
                    Lck.Spot{iTrial,1}= SpotId;
                    Lck.animalname{iTrial,1}= Settings.Animal;
                    Lck.Cond{iTrial,1} = subDirsNames{1,iCondition};
                    Lck.pixelsize{iTrial,1} = CellScans(iScan).rawImg.metadata.pixelSize;
                    
                    % get fraction of active MD area
                    neuroMask = CellScans(2).calcFindROIs.data.roiMask;
                    if ndims(neuroMask) == 3
                        neuroMask = max(neuroMask, [], 3);
                    end
                    [nFluoPix, nActivePix, nTotalPix] = ...
                        getFracActive(CellScans(1), 'nanMask', ...
                        neuroMask);
                    Lck.nFluoPix{iTrial,1} = nFluoPix;
                    Lck.nActivePix{iTrial,1} = nActivePix;
                    Lck.nTotalPix{iTrial,1} = nTotalPix;
                    
                    % ROI masks output
                    Lck.Trial_ROIMask{iTrial,1} = mask;
                    
                    % all astrocyte ROI masks together
                    trialmasks(:,:,iTrial) = mask;
                    
                end
            end
            
            
        end
        
        % append peak data
        dataNames=fieldnames(data);
        data2= struct2cell(data);
        data3= [data2{:}];
        
        AllData=vertcat(AllData, data3);
        
        clearvars data data2 data3
        
        % Mask of Active Pixels (proportional to number of trials)
        sumImg =sum(trialmasks,3);
        fracImg = sum(trialmasks,3)./numel(TrialNames);
        
        % calculate the score for responses (normalized to the "threshold"
        thresh = 1/numel(TrialNames); % or 1/numel(trials)
        activePxIdx = fracImg > 0;
        activePx = sum(activePxIdx(:));
        fracImg(fracImg <= thresh) = NaN;
        score = nansum(nansum(fracImg)) / activePx;
        
        % save the mask of active pixels
        FigFileName1 = fullfile(FolderName,'Sum_Pixel_Mask.tif');
        FigFileName2 = fullfile(FolderName,'FracActive_Pixel_Mask.tif');
        
        f = figure('visible', 'off');
        imagesc(sumImg);
        axis square
        caxis([0 1])
        colormap('jet')
        colorbar
        saveas(gcf,FigFileName1{1,1})
        close(f)
        
        
        f = figure('visible', 'off');
        imagesc(fracImg);
        axis square
        caxis([0 1])
        colormap('jet')
        colorbar
        saveas(gcf,FigFileName2{1,1})
        close(f)
        
        % store Lck field of view information
        
        for iLck= 1:length(Lck.trialname)
            Lck.Total_ROIMask{iLck,1} = fracImg;
            Lck.Response_Score{iLck,1}=score;
        end
        
        LckNames=fieldnames(Lck);
        Lckdata2= struct2cell(Lck);
        Lckdata3= [Lckdata2{:}];
        
        LckData=vertcat(LckData, Lckdata3);
        
        clearvars Lck Lckdata2 Lckdata3
    end
end


% %% Save all data for R analysis
% peak data
AllData2= [dataNames';AllData];

%Field of View data
LckFieldData= [LckNames';LckData];
LckFieldData2=LckFieldData;
LckFieldData(:,9:10)=[];

%onsetTimeTable (ROI table)
names={'ROI','Trial','Channel','Spot','Animal', 'Condition','baseline',...
    'trace','ROIIdx','PixelSize','FrameRate','OnsetTime','TraceAUC1','TraceAUC10'};
All_traces2=vertcat(names, All_traces);
All_traces2(:,8:9)=[];



cd(fullfile(Settings.ResultsFolder));

% write date to created file
save(SaveFiles{1,1}, 'AllData2','-v7.3');
cell2csv(SaveFiles{1,2}, AllData2);
save(SaveFiles{1,3}, 'All_traces','-v7.3');
cell2csv(SaveFiles{1,4}, All_traces2);
save(SaveFiles{1,5}, 'LckFieldData2','-v7.3');
cell2csv(SaveFiles{1,6}, LckFieldData);


close all; clear variables;

