close all; clear all;

All_traces=[]; AllData=[]; All_traces2=[]; AllData2=[];
%% Information about your images

% folder where data should be saved for each animal
Settings.ResultsFolder = 'D:\Data\test';

% File names for saving
SaveFiles{1,1} = fullfile(Settings.ResultsFolder,'FilesforR', 'Crazy8_peaks_longtrials_07_2019.csv');%'Control_Peaks_3Conds.csv'); % all data
SaveFiles{1,2} = fullfile(Settings.ResultsFolder,'FilesforMatlab', 'Crazy8_peaks_longtrials_07_2019.mat');%'Control_Peaks_3Conds.csv'); % all data
SaveFiles{1,3}= fullfile(Settings.ResultsFolder,'FilesforMatlab','Crazy8_traces_longtrials_07_2019.mat'); %'Control_TraceAUC_20sWindow_3Conds.csv'); % neuronal hand click cell scan
SaveFiles{1,4}= fullfile(Settings.ResultsFolder,'FilesforR','Crazy_onset_time_longtrials_07_2019.csv');

Settings.FileNames = {  % folder names where images and roiSet.zip is found
    'J:\Jill_Stobart\In_vivo_2P_Data\66678_Crazy8\2019_06_14\spot1',...
    'J:\Jill_Stobart\In_vivo_2P_Data\66678_Crazy8\2019_06_14\spot2',...
    % etc.
    };

Settings.Baseline = 2; % time (s) before the whisker stimulator starts,  2s for short trials (10s long), 5s for long trials
Settings.Animal = 'Alice';

channel = struct('Ca_Neuron',1,'Ca_Memb_Astro',2);  % can change 'blank' to any channel as this doesn't matter


%% Loop through each file and make Cell Scans for LckGCaMP and RCaMP

for iSpot= 1:length(Settings.FileNames)
    
    SpotRoot= Settings.FileNames{iSpot};
    SpotId=SpotRoot(end-15:end);
    
    % Get a list of all files and folders in this folder.
    % different stimulation conditions
    ConditionFolders = dir(SpotRoot);  % look for folders
    ConditionFolders(ismember( {ConditionFolders.name}, {'.', '..'})) = [];
    % Get a logical vector that tells which is a directory.
    names    = {ConditionFolders.name};
    dirFlags = [ConditionFolders.isdir];
    % Extract only those that are directories.
    subDirsNames = names(dirFlags);
    
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
        
        for iTrial= 1:length(TrialNames)
            % Get the T-series path from the folder
            testRoot =fullfile(FolderName, TrialNames{iTrial});
            
            expfiles = dir(fullfile(testRoot{1,1},'*.xml'));  % look for the xml file
            fnList = fullfile(testRoot{1,1}, expfiles.name);
            
            % Load data with BioFormats
            ImgArray =  BioFormats(fnList, channel);
            
            %% Motion Correction
            if iTrial==1 && iCondition==1
                refImg=mean(ImgArray.rawdata(:,:,1,1:3),4);
            end
            
            ImgArray = ImgArray.motion_correct('refImg',refImg, 'ch',1,'maxShift', 10,'minCorr', 0.4);%,'doPlot',true);
            
            
            BL_frames= Settings.Baseline*ImgArray.metadata.frameRate; % number of baseline frames
            
            %% Configs for Lck GCaMP Cell Scan
            %ASTROCYTES
            % 2D automated selection for peaks
            AC_findConf = ConfigFindROIsFLIKA_2D.from_preset('ca_memb_astro', 'baselineFrames',...
                BL_frames,'freqPassBand',0.5,'sigmaXY', 2,...
                'sigmaT', 0.14,'thresholdPuff', 7, 'threshold2D', 0.1,...
                'minRiseTime',0.14, 'maxRiseTime', 1,'minROIArea', 10,...
                'dilateXY', 4, 'dilateT', 0.3,'erodeXY', 1, 'erodeT', 0.1);
            
            % measure ROIs (extract the traces)
            AC_measureConf = ConfigMeasureROIsDummy();
            
            % filter the traces to detect the peaks and get info about them
            AC_detectConf = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...
                'propagateNaNs', false, 'excludeNaNs', false, 'lpWindowTime', 1.5, 'spFilterOrder', 2,...
                'spPassBandMin',0.05, 'spPassBandMax', 0.5, 'thresholdLP', 3,'thresholdSP', 5);
            
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
                'dilateXY', 4, 'dilateT', 0.2,'erodeXY', 1, 'erodeT', 0.1,...
                'discardBorderROIs',false);
            
            % measure ROIs (extract the traces)
            N_measureConf = ConfigMeasureROIsDummy();
            
            % filter the traces to detect the peaks and get info about them
            N_detectConf = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...
                'propagateNaNs', false,'excludeNaNs', false, 'lpWindowTime', 5, 'spFilterOrder', 2,...
                'spPassBandMin',0.1, 'spPassBandMax', 1, 'thresholdLP', 5,'thresholdSP', 4);
            
            
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
            astrocytes.plot();
            neurons1.plot();
            neurons2.plot();
            
            %astrocytes.opt_config()
            %neurons2.opt_config()
            %neurons2.plot('signals');
            
            %% Extract/Calculate data we want
            %% Output data
            
            % make a giant data table
            listFields = {'amplitude', 'fullWidth', 'halfWidth', ...
                'numPeaks', 'peakTime', 'peakStart', 'peakStartHalf', ...
                'peakType', 'prominence', 'roiName', 'peakAUC'};
            
            CellScans=vertcat(astrocytes, neurons1, neurons2);
            
            
            % loop through cellscans
            for iScan=1:size(CellScans,1)
                
                % get fraction of active MD area
                neuroMask = CellScans(2).calcFindROIs.data.roiMask;
                if ndims(neuroMask) == 3
                    neuroMask = max(neuroMask, [], 3);
                end
                [nFluoPix, nActivePix, nTotalPix] = ...
                    getFracActive(CellScans(1), 'nanMask', ...
                    neuroMask);
                
                % peak output
                temp=CellScans(iScan).calcDetectSigs.data;
                temp2.trialname ={};
                temp2.animalname = {};
                temp2.channel = {};
                temp2.Spot = {};
                temp2.Cond = {};
                temp2.area = {};
                temp2.pixelsize = {};
                temp2.nFluoPix = {};
                temp2.nActivePix = {};
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
                    temp2.Cond{iPeak,1} = subDirsNames(iCondition);
                    temp2.pixelsize{iPeak,1} = CellScans(iScan).rawImg.metadata.pixelSize;
                    temp2.nFluoPix{iPeak,1} = nFluoPix;
                    temp2.nActivePix{iPeak,1} = nActivePix;
                    
                    % get the indices  and area for a particular ROI
                            if iScan==3 || iScan==1
                                jROIname = temp.roiName{iPeak};
                                ROIindex= strcmp(CellScans(iScan,1).calcFindROIs.data.roiNames,jROIname);
                                temp2.area{iPeak,1} = CellScans(iScan,1).calcFindROIs.data.area(ROIindex);
                                
                                if iScan==3
                                jx = CellScans(iScan,1).calcFindROIs.data.centroidX(ROIindex,1);
                                jy = CellScans(iScan,1).calcFindROIs.data.centroidY(ROIindex,1);
                                xclose = round(x_pix*0.03); % 3 percent of pixels
                                yclose = round(y_pix*0.03); % 3 percent of pixels
                                
                                if isempty(jx)
                                    temp2.overlap{iPeak,1} = 0;
                                else
                                    % find FLIKA ROIs and hand selected ROIS that have
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
                                    temp2.overlap{iPeak,1} = 0;
                                end
                            else
                                temp2.area{iPeak,1} = 0;
                                
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
                    data.nFluoPix = {};
                    data.nActivePix = {};
                    data.overlap ={};
                end
                data.Trial= [data.Trial; temp2.trialname];
                data.Animal= [data.Animal; temp2.animalname];
                data.Channel= [data.Channel; temp2.channel];
                data.Spot= [data.Spot; temp2.Spot];
                data.Condition= [data.Condition; temp2.Cond];
                data.area= [data.area; temp2.area];
                data.pixelsize= [data.pixelsize; temp2.pixelsize];
                data.nFluoPix = [data.nFluoPix; temp2.nFluoPix];
                data.nActivePix = [data.nActivePix; temp2.nActivePix];
                data.overlap= [data.overlap; temp2.overlap];
                
                clearvars temp temp2
                
                
                %% make a table of trace info
                
                %traces output processes
                if strcmp(CellScans(iScan).calcFindROIs.data.roiNames{1,1}, 'none')
                    continue
                else
                    traces= CellScans(iScan).calcMeasureROIs.data.tracesNorm;
                    %preallocate
                    Trace_data=cell(size(traces,2),10);
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
                        Trace_data{iROI,6}= subDirsNames(iCondition);
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
                    
                      % look for overlap of neuronal ROIs
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
                  All_traces=vertcat(All_traces, Trace_data);
                clearvars Trace_data                      
                end
                                
            end
            
            
            dataNames=fieldnames(data);
            data2= struct2cell(data);
            data3= [data2{:}];
            
            AllData=vertcat(AllData, data3);
            
            
            clearvars data data3 data2
        end
    end
end


% %% Save all data for R analysis
AllData2= [dataNames';AllData];

%onsetTimeTable
names={'ROI','Trial','Channel','Spot','Animal', 'Condition','baseline',...
    'trace','ROIIdx','PixelSize','FrameRate','OnsetTime','TraceAUC8','overlap'};
All_traces2=vertcat(names, All_traces);
All_traces2(:,8)=[];
All_traces2(:,8)=[];

cd(fullfile(Settings.ResultsFolder));
% write date to created file
cell2csv(SaveFiles{1,1}, AllData2);
save(SaveFiles{1,3}, 'All_traces','-v7.3');
save(SaveFiles{1,2}, 'AllData2','-v7.3');
cell2csv(SaveFiles{1,4}, All_traces2);


