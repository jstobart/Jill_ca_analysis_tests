close all; clear all;

All_traces=[]; AllData=[]; All_traces2=[]; AllData2=[];
%% Information about your images

% folder where data should be saved for each animal
Settings.ResultsFolder = 'E:\Jill\Data\Winnipeg\Pericytes\Brain Slice 2P Data\Results';

% File names for saving
SaveFiles{1,1} = fullfile(Settings.ResultsFolder,'GCaMPSlice_peaks_02_02_2020_baselineNimodipine.csv');
SaveFiles{1,2} = fullfile(Settings.ResultsFolder,'GCaMPSlice_peaks_02_02_2020_baselineNimodipine.mat');
SaveFiles{1,3}= fullfile(Settings.ResultsFolder,'GCaMPSlice_traces_02_02_2020_baselineNimodipine.mat');

Drug='Nimodipine';
Settings.FileNames = {  % folder names where images and roiSet.zip is found
    'E:\Jill\Data\Winnipeg\Pericytes\Brain Slice 2P Data\PDGFRB-GCaMP\2019_07_11_nimodipine_bleach\baseline',...
    'E:\Jill\Data\Winnipeg\Pericytes\Brain Slice 2P Data\PDGFRB-GCaMP\2019_07_11_nimodipine_bleach\nimodipine',...
    'E:\Jill\Data\Winnipeg\Pericytes\Brain Slice 2P Data\PDGFRB-GCaMP\2019_07_12_nimodipine\baseline',...
    'E:\Jill\Data\Winnipeg\Pericytes\Brain Slice 2P Data\PDGFRB-GCaMP\2019_07_12_nimodipine\nimodipine',...
    'E:\Jill\Data\Winnipeg\Pericytes\Brain Slice 2P Data\PDGFRB-GCaMP\2019_07_15_nimodipine\baseline',...
    'E:\Jill\Data\Winnipeg\Pericytes\Brain Slice 2P Data\PDGFRB-GCaMP\2019_07_15_nimodipine\nimodipine',...
    'E:\Jill\Data\Winnipeg\Pericytes\Brain Slice 2P Data\PDGFRB-GCaMP\2019_07_16_nimodipine\baseline',...
    'E:\Jill\Data\Winnipeg\Pericytes\Brain Slice 2P Data\PDGFRB-GCaMP\2019_07_16_nimodipine\nimodipine',...
    'E:\Jill\Data\Winnipeg\Pericytes\Brain Slice 2P Data\PDGFRB-GCaMP\2019_07_17_nimodipine\baseline',...
    'E:\Jill\Data\Winnipeg\Pericytes\Brain Slice 2P Data\PDGFRB-GCaMP\2019_07_17_nimodipine\nimodipine',...
    };

Settings.Baseline = 5; % time (s) before the whisker stimulator starts,  2s for short trials (10s long), 5s for long trials

channel = struct('Ca_Cyto_Astro',1);  % can change 'blank' to any channel as this doesn't matter


%% Loop through each file and make Cell Scans for LckRCaMP and RCaMP

for iSpot= 1:length(Settings.FileNames)
    
    SpotRoot= Settings.FileNames{iSpot};
    
    %is this a baseline trial?
    if ~isempty(strfind(SpotRoot,'baseline'))
        Condition='baseline';
    else
        Condition=Drug;
    end
    
    startId=strfind(SpotRoot,'2019');
    SpotId=SpotRoot(66:75);  % MAKE SURE FOLDERS ARE ALL NAMED THE SAME WAY
    
    % Get a list of all files and folders in this folder.
    TrialFolders = dir(SpotRoot);  % look for folders
    TrialFolders(ismember( {TrialFolders.name}, {'.', '..'})) = [];
    % Get a logical vector that tells which is a directory.
    names2    = {TrialFolders.name};
    dirFlags2 = [TrialFolders.isdir];
    % Extract only those that are directories.
    TrialNames = names2(dirFlags2);
    
    for iTrial= 1:length(TrialNames)
        % Get the T-series path from the folder
        testRoot =fullfile(SpotRoot, TrialNames{iTrial});
        
        expfiles = dir(fullfile(testRoot,'*.xml'));  % look for the xml file
        fnList = fullfile(testRoot, expfiles.name);
        
        % Load data with BioFormats
        ImgArray =  BioFormats(fnList, channel);
        
        %% Motion Correction
        if iTrial==1
            refImg=mean(ImgArray.rawdata(:,:,1,1:3),4);
        end
        
       % ImgArray = ImgArray.motion_correct('refImg',refImg, 'ch',1,'maxShift', 10,'minCorr', 0.4, 'inpaintIters', 10);%,'doPlot',true);
        
        
        %frames1=round(35*ImgArray.metadata.frameRate);
        %frames2=round(50*ImgArray.metadata.frameRate);
        %frames3=round(100*ImgArray.metadata.frameRate);
        %frames4=round(125*ImgArray.metadata.frameRate);
        
        %         ImgArray = ImgArray.exclude_frames([frames1:frames2]);
        %
        %[Img1, ~] = split1(ImgArray, 4, [frames1 size(ImgArray.rawdata,4)-frames1]);
        %[~, Img2] = split1(ImgArray, 4, [frames2 size(ImgArray.rawdata,4)-frames2]);
        
        
        
        BL_frames= round(Settings.Baseline*ImgArray.metadata.frameRate); % number of baseline frames
        
        
        %% Configs for RCaMP Cell Scan
        % RCaMP
        % hand selected- peaks from cellular structures
        x_pix= size(ImgArray(1,1).rawdata,2); y_pix= size(ImgArray(1,1).rawdata,1);
        scaleF = 1;
        
        zipfiles = dir(fullfile(SpotRoot,'*.zip'));  % look for the zip ROI file
        fnTempList2 = {zipfiles(:).name};
        zipPath = fullfile(SpotRoot, fnTempList2);
        
        % load image J ROIs for mask
        findConf{1} = ConfigFindROIsDummy.from_ImageJ(zipPath{1,1}, x_pix, y_pix, scaleF);
        
        % 2D FLIKA selected for peaks
        %          findConf{2} = ConfigFindROIsFLIKA_2D.from_preset('ca_cyto_astro', 'baselineFrames',...
        %             BL_frames, 'freqPassBand',0.05,'sigmaXY', 1,...
        %             'sigmaT', 1,'thresholdPuff', 3, 'threshold2D', 0,...
        %             'minRiseTime',0.15, 'maxRiseTime', 1,'minROIArea', 10,...
        %             'dilateXY', 10, 'dilateT', 2,'erodeXY', 5, 'erodeT', 1,...
        %             'discardBorderROIs',false);
        
        % measure ROIs (extract the traces)
        measureConf = ConfigMeasureROIsDummy('baselineFrames', BL_frames);
        
        % filter the traces to detect the peaks and get info about them
        detectConf = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...
            'propagateNaNs', true,'excludeNaNs', false, 'lpWindowTime', 5,...
            'thresholdLP', 5,'thresholdSP', 5,'spPassBandMin', 0.025, 'spPassBandMax', 0.4);
        
        % GCaMP parameters BandMin = 0.025, BandMax = 0.6
        
        % Combine the configs into a CellScan config for neuronal RCaMP
        configCS_ImageJ= ConfigCellScan(findConf{1}, measureConf, detectConf); %
        %configCS_FLIKA= ConfigCellScan(findConf{2}, measureConf, detectConf); %
        
        
        %% Create CellScan objects
        RCaMP1 = CellScan(fnList, ImgArray, configCS_ImageJ, 1); % peaks from hand clicked ROIs
        %RCaMP2 = CellScan(fnList, Img4, configCS_ImageJ, 1);
        
        %RCaMP2 = CellScan(fnList, ImgArray, configCS_FLIKA, 1); % peaks from automated ROIs
        
        % Process the images
        RCaMP1 =RCaMP1.process();
        %RCaMP2 =RCaMP2.process();
        
        % Make the debugging plots
        RCaMP1.plot();
        %RCaMP2.plot();
        
        %RCaMP1.opt_config()
        %RCaMP2.opt_config()
        %RCaMP2.plot('signals');
        
        %% Extract/Calculate data we want
        %% Output data
        
        % make a giant data table
        listFields = {'amplitude', 'fullWidth', 'halfWidth', ...
            'numPeaks', 'peakTime', 'peakStart', 'peakStartHalf', ...
            'peakType', 'prominence', 'roiName', 'peakAUC'};
        
        %CellScans=vertcat(RCaMP1, RCaMP2);
        CellScans=RCaMP1;
        
        % loop through cellscans
        for iScan=1:size(CellScans,1)
            
            
            % peak output
            temp=CellScans(iScan).calcDetectSigs.data;
            temp2.trialname ={};
            temp2.Condition = {};
            temp2.channel = {};
            temp2.Spot = {};
            temp2.area = {};
            temp2.pixelsize = {};
            
            % extract fields from Class
            for jField = 1:numel(listFields)
                isFirst = iScan ==1;
                if isFirst
                    data.(listFields{jField}) = {};
                end
                data.(listFields{jField}) = [data.(listFields{jField}); ...
                    temp.(listFields{jField})];
            end
            
            % create fields for trial, animal, spot, condition, etc.
            for iPeak = 1:length(temp.amplitude)
                temp2.trialname{iPeak,1}=strcat('trial', num2str(iTrial,'%02d'));
                temp2.channel{iPeak,1}= 'GCaMP';
                temp2.Spot{iPeak,1}= SpotId;
                temp2.Condition{iPeak,1}= Condition;
                temp2.pixelsize{iPeak,1} = CellScans(iScan).rawImg.metadata.pixelSize;
                
                % get the indices  and area for a particular ROI
                if iScan==2
                    jROIname = temp.roiName{iPeak};
                    ROIindex= strcmp(CellScans(iScan,1).calcFindROIs.data.roiNames,jROIname);
                    temp2.area{iPeak,1} = CellScans(iScan,1).calcFindROIs.data.area(ROIindex);
                else
                    temp2.area{iPeak,1} = 0;
                    
                end
                
            end
            
            
            
            %%
            isFirst = iScan == 1;
            if isFirst
                data.Trial = {};
                data.Condition = {};
                data.Channel = {};
                data.Spot = {};
                data.area = {};
                data.pixelsize={};
            end
            data.Trial= [data.Trial; temp2.trialname];
            data.Condition= [data.Condition; temp2.Condition];
            data.Channel= [data.Channel; temp2.channel];
            data.Spot= [data.Spot; temp2.Spot];
            data.area= [data.area; temp2.area];
            data.pixelsize= [data.pixelsize; temp2.pixelsize];
            
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
                    Trace_data{iROI,3}= 'GCaMP';
                    
                    Trace_data{iROI,4}= SpotId;
                    Trace_data{iROI,5}= Condition;
                    Trace_data{iROI,6} = Settings.Baseline;
                    Trace_data{iROI,7} = traces(:,iROI);
                    if iScan==2
                        Trace_data{iROI,8} = CellScans(iScan).calcFindROIs.data.roiIdxs{iROI,1};
                        Trace_data{iROI,9} = CellScans(iScan).rawImg.metadata.pixelSize;
                    else
                        Trace_data{iROI,8} = CellScans(iScan).calcFindROIs.data.roiMask(:,:,iROI);
                        Trace_data{iROI,9} = CellScans(iScan).rawImg.metadata.pixelSize;
                    end
                    
                    FrameRate= CellScans(1).rawImg.metadata.frameRate;
                    Trace_data{iROI,10} = FrameRate; % frameRate
                    
                    % trace AUC
                    Trace_data{iROI,11}=trapz(traces(:,iROI));
                    
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


% %% Save all data for R analysis
AllData2= [dataNames';AllData];

%Traces Table
names={'ROI','Trial','Channel','Spot','Condition', 'baseline',...
    'trace','ROIIdx','PixelSize','FrameRate','TraceAUC'};

cd(fullfile(Settings.ResultsFolder));
% write date to created file
cell2csv(SaveFiles{1,1}, AllData2);
save(SaveFiles{1,3}, 'All_traces','-v7.3');
save(SaveFiles{1,2}, 'AllData2','-v7.3');



