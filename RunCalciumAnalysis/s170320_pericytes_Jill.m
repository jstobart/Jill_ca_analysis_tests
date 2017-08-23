close all; clear all;

AllData= [];
All_traces= [];
AllPropData=[];
%% Example script for CellScan application
% Add path to 2p-img-analysis
%addpath('D:\Code\Matlab\2p-img-analysis');
% Add path to experiment directory
%addpath(genpath('D:\Experiments\Tendon Biomechanics'));

%cd(utils.CHIPS_rootdir)


%% Information about your images
Settings.MainDir = 'E:\Data\Pericyte_project\Two-photon-data\Calcium';

Settings.AnimalNames = {
    'Mouse1',...
    };
Settings.ScoreSheetNames = {
    'Mouse1_CalciumFilesScoresheet.xlsx',...
    };

channel = struct('Ca_Memb_Astro', 1, 'blood_plasma', 2);

% final data file name
SaveFiles{1,1} = fullfile(Settings.MainDir, 'Results', 'test_peaks_18_08_2017.csv');% all peak data
SaveFiles{1,2} = fullfile(Settings.MainDir, 'Results', 'test_prop_18_08_2017.csv');% propagation data (duration, distance, etc.)
SaveFiles{1,3} = fullfile(Settings.MainDir, 'Results', 'test_prop_ROIdata_18_08_2017.mat');%'propagation ROI info
SaveFiles{1,4}= fullfile(Settings.MainDir, 'Results','test_traces_18_08_2017.mat'); % normalized traces

%% Load calibration file
calibration ='E:\matlab\2p-img-analysis\tests\res\calibration_20x.mat';
CalFile = Calibration_PixelSize.load(calibration);

%fNameCal = fullfile(utils.CHIPS_rootdir, 'tests/res', 'calibration_20x.mat');
%CalFile = Calibration_PixelSize.load(fNameCal);

%% make a mixing matrix for spectral unmixing of RCaMP and GCaMP

if ~exist(fullfile(Settings.MainDir, 'Results','eGFP_TexasRed_Matrix.mat'),'file')
    % load example high res image
    unmixFile = 'E:\Data\Pericyte_project\Two-photon-data\Calcium\Mouse1\cell2\cell2_940nm_zoom15_highres009.tif';
    unmixImg = SCIM_Tif(unmixFile, channel, CalFile);
    
    [unmixImg, eGFP_TexasRed_Matrix] = unmix_chs(unmixImg); %returns the mixing matrix to be used for all imaging
    
    % save the mixing matrix for loading later
    cd(fullfile(Settings.MainDir, 'Results'));
    % write matrix to created file
    save('eGFP_TexasRed_Matrix.mat', 'eGFP_TexasRed_Matrix');
else
    load(fullfile(Settings.MainDir, 'Results','eGFP_TexasRed_Matrix.mat'));
end


%% load scoresheet and loop through animal, spot, etc.

Settings.ScoreSheetPath = fullfile(Settings.MainDir,Settings.ScoreSheetNames);

numAnimals = length(Settings.AnimalNames);
for iAnimal = 1:numAnimals
    CurrentAnimal = Settings.AnimalNames{iAnimal};
    %read Scoresheet of current animal
    CurrentSheet = Settings.ScoreSheetPath{iAnimal};
    Settings = readScoresheet(CurrentSheet, Settings);
    
    % Get SpotID
    spots = Settings.SpotIDs;
    
    for iSpot = 1:length(spots)
        
        %% load data
        % Extract spot name
        spotId = spots{iSpot};
        
        % Find the idx of paths matching this spot
        CurrentDepth = Settings.Depth(iSpot); %depth
        CurrentCell = Settings.CellType(iSpot); %pericyte type
        
        % Get image paths
        testRoot =Settings.LowresPath{iSpot};
        
        expfiles = dir(fullfile(testRoot,'*.tif'));
        fnTempList = {expfiles(:).name};
        fnList = fullfile(testRoot, fnTempList);
        
        % Create an array of ScanImage Tiffs
        ImgArray =  SCIM_Tif(fnList, channel, CalFile);
        
        %% Run preprocessing steps
        % Spectral Unmixing
        ImgArray= ImgArray.unmix_chs(false, [], cell2mat(eGFP_TexasRed_Matrix));
        
        %Motion correction
        channelToUseMC = 2; % which channel to use
        refImg = squeeze(mean(ImgArray(1,1).rawdata(:,:,channelToUseMC, 5:10),4));
        ImgArray=ImgArray.motion_correct( 'refImg', refImg,'ch', channelToUseMC,'minCorr', 0.4);
        
        
        %% Configs for Finding ROIs
        % Pericyte calcium
        
        % automated selection for peak detection
%         findConf{1} = ConfigFindROIsFLIKA_2p5D.from_preset('ca_memb_astro','freqPassBand',1,'sigmaXY', 2,...
%             'sigmaT', 0.1,'threshold_std', 7, 'threshold2D', 0.2,...
%             'min_rise_time',0.16, 'max_rise_time', 1.5,'minPuffArea', 4,...
%             'minPuffTime', 0.25,'dilateXY', 2, 'dilateT', 0.1,'erodeXY', 1, 'erodeT', 0.1,...
%             'discardBorderROIs',false);
        
        % automated selection for signal propagation
        findConf{2} = ConfigFindROIsFLIKA_3D.from_preset('ca_memb_astro', 'freqPassBand',1,'sigmaXY', 2,...
            'sigmaT', 0.1,'threshold_std', 7, 'threshold2D', 0.2,...
            'min_rise_time',0.16, 'max_rise_time', 1.5,'minPuffArea', 4,...
            'minPuffTime', 0.25,'dilateXY', 2, 'dilateT', 0.1,'erodeXY', 1, 'erodeT', 0.1,...
            'discardBorderROIs',false);
        
        % hand selected
        x_pix= 128; y_pix= 127;
        scaleF = 1;
        % zipPath = 'D:/....';
        zipPath= fullfile(testRoot,'RoiSet.zip');
        findConf{3} = ConfigFindROIsDummy.from_ImageJ(zipPath, x_pix, y_pix, scaleF);
        
        
        %% Configs for measuring ROIs
        % AWAKE astrocyte membrane calcium
        detectConf = ConfigDetectSigsClsfy('propagateNaNs', false, 'excludeNaNs', false,...
            'lpWindowTime', 5, 'spFilterOrder', 2,'spPassBandMin',0.025, 'spPassBandMax', 1,...
            'thresholdSD_low', 7,'thresholdSD_band', 3);
        
        % for calculating AUC for each trace
        measureConf = ConfigMeasureROIsDummy();
        %         measureConf = ConfigMeasureROIsZScore(...
        %             'zIters', 10, ...
        %             'zSDs', 2, ...
        %             'propagateNaNs', false);
        
        % Combine the configs into a CellScan config
       % configCS{1,1} = ConfigCellScan(findConf{1,1}, measureConf,detectConf); % pericyte FLIKA, peaks
        configCS{1,2} = ConfigCellScan(findConf{1,2}, measureConf,detectConf); % pericyte FLIKA, signal propagation
        configCS{1,3} = ConfigCellScan(findConf{1,3}, measureConf,detectConf); % pericyte hand selected
        
        %% Create CellScan objects
        Hand_peaks = CellScan(fnList, ImgArray, configCS{1,3}, 1); % peaks from hand clicked ROIs
        
        FLIKA_prop = CellScan(fnList, ImgArray, configCS{1,2}, 1); % ROIs from 3D FLIKA
        
        
        %% Process the images
        Hand_peaks =Hand_peaks.process();
        FLIKA_prop =FLIKA_prop.process();
        
        
        %% Make the debugging plots
        %Hand_peaks.plot;
        %FLIKA_prop.plot;
        
        % for working out find peaks parameters
        %Hand_peaks(1).opt_config()
        
        %FLIKA_prop(1).opt_config()
        
        
        
        %% Output data
        
        % make a giant data table
        listFields1 = {'amplitude', 'fullWidth', 'halfWidth', ...
            'numPeaks', 'peakTime', 'peakStart', 'peakStartHalf', ...
            'peakType', 'prominence', 'ROIname', 'peakAUC'};
        
        % pericyte peaks
        for itrial=1:length(Hand_peaks)
            
            % peak output
            temp=Hand_peaks(1,itrial).calcDetectSigs.data;
            temp2.trialname ={};
            temp2.animalname = {};
            temp2.Spot = {};
            temp2.depth = {};
            temp2.celltype = {};
            temp2.pixelsize = {};
            
            % extract fields from Class
            for jField = 1:numel(listFields1)
                isFirst = (itrial == 1 );
                if isFirst
                    data.(listFields1{jField}) = {};
                end
                data.(listFields1{jField}) = [data.(listFields1{jField}); ...
                    temp.(listFields1{jField})];
            end
            
            % create fields for trial, animal, spot, condition, etc.
            for iPeak = 1:length(temp.amplitude)
                temp2.trialname{iPeak,1}=strcat('trial', num2str(itrial));
                temp2.Spot{iPeak,1}= spotId;
                temp2.celltype{iPeak,1}= CurrentCell(1,1);
                temp2.animalname{iPeak,1}= CurrentAnimal;
                temp2.depth{iPeak,1} = CurrentDepth(1,1);
                temp2.pixelsize{iPeak,1} = Hand_peaks(1,itrial).rawImg.metadata.pixelSize;
            end
            
            isFirst = (itrial == 1 );
            if isFirst
                data.Trial = {};
                data.Animal = {};
                data.Spot = {};
                data.CellType={};
                data.Depth = {};
                data.pixelsize={};
            end
            data.Trial= [data.Trial; temp2.trialname];
            data.Animal= [data.Animal; temp2.animalname];
            data.Spot= [data.Spot; temp2.Spot];
            data.CellType= [data.CellType; temp2.celltype];
            data.Depth= [data.Depth; temp2.depth];
            data.pixelsize= [data.pixelsize; temp2.pixelsize];
            
            %% extract traces from ROIs for correlations etc.
            traces= Hand_peaks(1,itrial).calcMeasureROIs.data.tracesNorm;
            %preallocate
            Trace_data=cell(size(traces,2),7);
            for iROI = 1:size(traces,2)
                Trace_data{iROI,1}= Hand_peaks(1,itrial).calcFindROIs.data.roiNames{iROI,1};
                Trace_data{iROI,2}= strcat('trial', num2str(itrial));
                Trace_data{iROI,3}= spotId;
                Trace_data{iROI,4}= CurrentAnimal;
                Trace_data{iROI,5}= CurrentCell;
                Trace_data{iROI,6} = CurrentDepth(1,1);
                Trace_data{iROI,7} = traces(:,iROI);
                Trace_data{iROI,8} = Hand_peaks(1,itrial).calcFindROIs.data.roiMask(:,:,iROI);
                Trace_data{iROI,9} = Hand_peaks(1,itrial).rawImg.metadata.pixelSize;
            end
            All_traces=vertcat(All_traces, Trace_data);
            clearvars Trace_data
            clearvars temp temp2
        end
        
        
        %% Pericyte signal propagation
        
        for itrial=1:length(FLIKA_prop)
            temp3=FLIKA_prop(1,itrial).calcFindROIs.data;
            
            % create fields for trial, animal, spot, condition, etc.
            for iROI = 1:length(temp3.distance)
                temp4.roiNames{iROI,1}=temp3.roiNames{iROI,1};
                temp4.trialname{iROI,1}=strcat('trial', num2str(itrial));
                temp4.Spot{iROI,1}= spotId;
                temp4.animalname{iROI,1}= CurrentAnimal;
                temp4.celltype{iROI,1}=CurrentCell;
                temp4.depth{iROI,1} = CurrentDepth(1,1);
                temp4.propRate{iROI,1} = temp3.distance(iROI,1)/temp3.duration(iROI,1); % rate of propagation= distance of propagation/duration (um/s)
                temp4.pixelsize{iROI,1} = Hand_peaks(1,itrial).rawImg.metadata.pixelSize;
                temp4.centroid{iROI,1}=temp3.centroid{iROI,1};
                temp4.puffIdxs{iROI,1}=temp3.puffIdxs{iROI,1};
                temp4.distance{iROI,1}=temp3.distance(iROI,1);
                temp4.duration{iROI,1}=temp3.duration(iROI,1);
                temp4.onset{iROI,1}=temp3.onset(iROI,1);
                temp4.volume{iROI,1}=temp3.volume(iROI,1);
                temp4.area{iROI,1}=temp3.area(iROI,1);
            end
            isFirst = (itrial == 1 );
            if isFirst
                propdata.Trial = {};
                propdata.Animal = {};
                propdata.Spot = {};
                propdata.CellType={};
                propdata.Depth = {};
                propdata.PropRate = {};
                propdata.pixelsize={};
                propdata.distance={};
                propdata.duration={};
                propdata.onset={};
                propdata.volume={};
                propdata.area={};
            end
            propdata.Trial= [propdata.Trial; temp4.trialname];
            propdata.Animal= [propdata.Animal; temp4.animalname];
            propdata.Spot= [propdata.Spot; temp4.Spot];
            propdata.PropRate= [propdata.PropRate; temp4.propRate];
            propdata.CellType= [propdata.CellType; temp4.celltype];
            propdata.Depth= [propdata.Depth; temp4.depth];
            propdata.pixelsize= [propdata.pixelsize; temp4.pixelsize];
            propdata.distance= [propdata.distance; temp4.distance];
            propdata.duration= [propdata.duration; temp4.duration];
            propdata.onset= [propdata.onset; temp4.onset];
            propdata.volume= [propdata.volume; temp4.volume];
            propdata.area= [propdata.area; temp4.area];
            clearvars temp3 temp4
        end
        
        % concatenate peak data
        dataNames1=fieldnames(data);
        data2= struct2cell(data);
        data3= [data2{:}];
        AllData=vertcat(AllData, data3);
        
        % concatenate peak data
        dataNames2=fieldnames(propdata);
        propdata2= struct2cell(propdata);
        propdata3= [propdata2{:}];
        AllPropData=vertcat(AllPropData, propdata3);
        
        
        clearvars data data3 propdata propdata3
        
    end
end


% %% Save all data for R analysis
AllData2= [dataNames1';AllData];
AllPropData2= [dataNames2';AllPropData];

%TraceNames= {'ROI','Trial','Channel','Spot','Animal','Condition','depth','traces','centroid','puffIdx'};
%All_traces2=[TraceNames;All_traces];
%
cd(fullfile(Settings.MainDir, 'Results'));
% write date to created file
cell2csv(SaveFiles{1,1}, AllData2);
cell2csv(SaveFiles{1,2}, AllPropData2);
save(SaveFiles{1,3}, 'AllPropData2','-v7.3');
save(SaveFiles{1,4}, 'All_traces','-v7.3');


