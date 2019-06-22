close all; clear all;

%% Information about your images

% folder where data should be saved for each animal
Settings.ResultsFolder = 'D:\Data\test';

Settings.FileNames = {  % folder names where images and roiSet.zip is found
    'J:\Jill_Stobart\In_vivo_2P_Data\66678_Crazy8\2019_06_14\spot1',...
    % etc.
    };

Settings.Baseline = 2; % time (s) before the whisker stimulator starts,  2s for short trials (10s long), 5s for long trials
Settings.Animal = 'Alice';

channel = struct('Ca_Neuron',1,'Ca_Memb_Astro',2);  % can change 'blank' to any channel as this doesn't matter


% final data file name
SaveFiles{1,1} = fullfile(Settings.MainDir, 'Results','FilesforR', 'Peaks_pharmacology_Lck_nostim_vs_longstim_12_2017.csv');%'Control_Peaks_3Conds.csv'); % all data
SaveFiles{1,2} = fullfile(Settings.MainDir, 'Results','FilesforMatlab', 'Peaks_pharmacology_Lck_nostim_vs_longstim_12_2017.mat');%'Control_Peaks_3Conds.csv'); % all data
SaveFiles{1,3}= fullfile(Settings.MainDir, 'Results','FilesforMatlab','Traces_pharmacology_Lck_nostim_vs_longstim_12_2017.mat'); %'Control_TraceAUC_20sWindow_3Conds.csv'); % neuronal hand click cell scan
SaveFiles{1,4}= fullfile(Settings.MainDir, 'Results','FilesforR','OnsetTimes_pharmacology_Lck_nostim_vs_longstim_12_2017.csv');

%% Loop through each file and make Cell Scans for LckGCaMP and RCaMP

for iSpot= 1:length(Settings.FileNames)
    
    SpotRoot= Settings.FileNames{iSpot};
    SpotId=SpotRoot(last-15:end);
    
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
            
            zipfiles = dir(fullfile(testRoot,'*.zip'));  % look for the zip ROI file
            fnTempList2 = {zipfiles(:).name};
            zipPath = fullfile(testRoot, fnTempList2);
            
            % load image J ROIs for mask
            N_findConf = ConfigFindROIsDummy.from_ImageJ(zipPath{1,1}, x_pix, y_pix, scaleF);
            
            % measure ROIs (extract the traces)
            N_measureConf = ConfigMeasureROIsDummy();
            
            % filter the traces to detect the peaks and get info about them
            N_detectConf = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...
                'propagateNaNs', false,'excludeNaNs', false, 'lpWindowTime', 5, 'spFilterOrder', 2,...
                'spPassBandMin',0.1, 'spPassBandMax', 1, 'thresholdLP', 5,'thresholdSP', 4);
            
            
            % Combine the configs into a CellScan config for neuronal RCaMP
            N_configCS= ConfigCellScan(N_findConf, N_measureConf, N_detectConf); %
            
            
            %% Create CellScan objects
            astrocytes = CellScan(fnList, ImgArray, AC_configCS, 2); % peaks from automated ROIs
            neurons = CellScan(fnList, ImgArray, N_configCS, 1); % peaks from hand clicked ROIs
            
            % Process the images
            astrocytes =astrocytes.process();
            neurons =neurons.process();
            
            % Make the debugging plots
            astrocytes.plot();
            neurons.plot();
            
            %astrocytes.opt_config()
            %neurons.opt_config()
            %neurons.plot('signals');
            
            %% Output data
%             OutputFile1=fullfile(Settings.ResultsFolder, strcat('Astrocytes_FieldofView',num2str(iFile)));
%             OutputFile2=fullfile(Settings.ResultsFolder, strcat('Neurons_FieldofView',num2str(iFile)));
%             
%             astrocytes.output_data(OutputFile1)
%             neurons.output_data(OutputFile2)
            
            % output the different levels of data
            % output traces
            % output ROI masks
            % output onsets
            % output peaks
            % csv and matlab files.....
            
            % get fraction of active astrocyte MD area in each field of view
            
            % first create a mask of where the neurons are
            neuroMask = neurons.calcFindROIs.data.roiMask;
            if ndims(neuroMask) == 3
                neuroMask = max(neuroMask, [], 3);
            end
            
            % which astrocyte pixels are active?
            [nFluoPix, nActivePix, nTotalPix] = ...
                getFracActive(astrocytes, 'nanMask', ...
                neuroMask);
            
            % output the peak data
            
            % peak output
            peakstempAC=astrocytes.calcDetectSigs.data;
            peakstempN= neurons.calcDetectSigs.data;
            
            peaks.trialname ={};
            peaks.channel={};
            peaks.animalname = {};
            peaks.Spot = {};
            peaks.Cond = {};
            peaks.area = {};
            peaks.pixelsize = {};
            peaks.nFluoPix = {};
            peaks.nActivePix = {};
            
            % extract fields from Class
            for jField = 1:numel(listFields)
                isFirst = (iTrial == 1);
                if isFirst
                    peakdata.(listFields{jField}) = {};
                end
                peakdata.(listFields{jField}) = [peakdata.(listFields{jField}); ...
                    tempAC.(listFields{jField})];
            end
            
            % create fields for trial, animal, spot, condition, etc.
            for iPeak = 1:length(tempAC.amplitude)
                peaks.trialname{iPeak,1}=strcat('trial', num2str(iTrial,'%02d'));
                peaks.channel='GCaMP';
                peaks.Spot{iPeak,1}= SpotId;
                peaks.animalname{iPeak,1}= Settings.Animal;
                peaks.Cond{iPeak,1} = subDirsNames(iCondition);
                peaks.pixelsize{iPeak,1} = astrocytes.rawImg.metadata.pixelSize;
                peaks.nFluoPix{iPeak,1} = nFluoPix;
                peaks.nActivePix{iPeak,1} = nActivePix;
                
                % get the area for a particular ROI
                
                    jROIname = tempAC.roiName{iPeak};
                    ROIindex= strcmp(astrocytes.calcFindROIs.data.roiNames,jROIname);
                    peaks.area{iPeak,1} = astrocytes.calcFindROIs.data.area(ROIindex);
                    
                
            end       
            
            %%
            isFirst = (iTrial == 1);
            if isFirst
                peakdata.Trial = {};
                peakdata.Animal = {};
                peakdata.Channel = {};
                peakdata.Spot = {};
                peakdata.Condition = {};
                peakdata.area = {};
                peakdata.pixelsize={};
                peakdata.nFluoPix = {};
                peakdata.nActivePix = {};
            end
            peakdata.Trial= [peakdata.Trial; peaks.trialname];
            peakdata.Animal= [peakdata.Animal; peaks.animalname];
            peakdata.Channel= [peakdata.Channel; peaks.channel];
            peakdata.Spot= [peakdata.Spot; peaks.Spot];
            peakdata.Condition= [peakdata.Condition; peaks.Cond];
            peakdata.area= [peakdata.area; peaks.area];
            peakdata.pixelsize= [peakdata.pixelsize; peaks.pixelsize];
            peakdata.nFluoPix = [peakdata.nFluoPix; peaks.nFluoPix];
            peakdata.nActivePix = [peakdata.nActivePix; peaks.nActivePix];
            
            clearvars tempAC peaks
            
            
            % make a table of trace info
            
            %traces output processes
            if astrocytes.calcFindROIs.data.roiNames{1,1}, 'none')
                continue
            else
                traces= astrocytes.calcMeasureROIs.data.tracesNorm;
                %preallocate
                Trace_data=cell(size(traces,2),10);
                for iROI = 1:size(traces,2)
                    Trace_data{iROI,1}= CellScans(iScan,itrial).calcFindROIs.data.roiNames{iROI,1};
                    Trace_data{iROI,2}= strcat('trial', num2str(iTrial,'%02d'));
                    Trace_data{iROI,3}= 'GCaMP';                  
                    Trace_data{iROI,4}= SpotId;
                    Trace_data{iROI,5}= Settings.tAnimal;
                    Trace_data{iROI,6}= subDirsNames(iCondition);
                    Trace_data{iROI,9} = Settings.Baseline;
                    Trace_data{iROI,10} = traces(:,iROI);
                    Trace_data{iROI,11} = astrocytes.calcFindROIs.data.roiIdxs{iROI,1};
                    Trace_data{iROI,12} = astrocytes.rawImg.metadata.pixelSize;
                                          
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
            
                    Trace_data{iROI,11} = neurons.calcFindROIs.data.roiMask(:,:,iROI);
                    
            All_traces=vertcat(All_traces, Trace_data);
            clearvars Trace_data
        end
    end
    
    
    dataNames=fieldnames(peakdata);
    data2= struct2cell(peakdata);
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




end
end
end




