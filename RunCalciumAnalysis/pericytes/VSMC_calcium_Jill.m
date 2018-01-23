close all; clear all;

AllData= [];
All_traces= [];

%% Information about your images
Settings.MainDir = 'E:\Data\Two_Photon_Data\Acta2-RCaMP';

Settings.ScoreSheetNames = {
    'Control_CalciumFilesScoresheet.xlsx',...
    };

channel = struct('Ca_Memb_Astro', 1);

Settings.NameConditions = {'Nostim'};

% final data file name
SaveFiles{1,1} = fullfile(Settings.MainDir, 'Results','FilesforR', 'Peaks_VSMC_01_2018.csv'); % peak data for R
SaveFiles{1,2} = fullfile(Settings.MainDir, 'Results','FilesforMatlab', 'Peaks_VSMC_01_2018.mat'); % peak dat for matlab
SaveFiles{1,3}= fullfile(Settings.MainDir, 'Results','FilesforMatlab','Traces_VSMC_01_2018.mat');
SaveFiles{1,4}= fullfile(Settings.MainDir, 'Results','FilesforR','OnsetTimes_VSMC_01_2018.csv');


doPlot=1;
BL_frames=30;
%% load scoresheet and loop through animal, spot, etc.

Settings.ScoreSheetPath = fullfile(Settings.MainDir,Settings.ScoreSheetNames);

numDrugs = length(Settings.ScoreSheetNames);
for iDrug = 1:numDrugs
    %read Scoresheet of current animal
    CurrentSheet = Settings.ScoreSheetPath{iDrug};
    Settings = readScoresheet(CurrentSheet, Settings);
    
    %%
    % Get SpotID
    spots = Settings.SpotIDs;
    
    for iSpot = 1:length(spots)
        
        %% load data
        % Extract spot name
        spotId = spots{iSpot};
        
        % Find the idx of paths matching this spot
        CurrentDepth = Settings.Depth(iSpot); %depth
        CurrentCell = Settings.CellType(iSpot); %pericyte type
        CurrentDrug = Settings.Drug(iSpot); %drug treatment
        CurrentAnimal = Settings.AnimalNames{iSpot};
        CurrentCondition= 'Nostim';
        
        
        %% Load calibration file
        if strcmp(Settings.Objective(iSpot),'20x')
            %calibration ='E:\matlab\2p-img-analysis\tests\res\calibration_20x.mat';
            calibration ='E:\matlab\CalibrationFiles\calibration_20x.mat';
        else
            calibration ='E:\matlab\CalibrationFiles\calibration_25x_approx.mat';
        end
        CalFile = CalibrationPixelSize.load(calibration);
        
        
        % Get image paths
        testRoot =Settings.LowresPath{iSpot};
        
        
        
        expfiles = dir(fullfile(testRoot,'lowres*'));
        fnTempList = {expfiles(:).name};
        fnList = fullfile(testRoot, fnTempList);
        
        % Create an array of ScanImage Tiffs
        ImgArray =  SCIM_Tif(fnList, channel, CalFile);
        
        
        channelToUseMC = Settings.MotionCorrChannel(iSpot); % which channel to use
        refImg = squeeze(mean(ImgArray(1,1).rawdata(:,:,channelToUseMC, 5:10),4));
        ImgArray=ImgArray.motion_correct('refImg', refImg,'ch', channelToUseMC,'minCorr', 0.4);
        
        
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
        AC_findConf{2} = ConfigFindROIsDummy.from_ImageJ(fullfile(testRoot,'RoiSet.zip'), x_pix, y_pix,1);
        
        
        % 3D automated selection for time and space estimations
        AC_findConf{3} = ConfigFindROIsFLIKA_3D.from_preset('ca_memb_astro', 'baselineFrames',...
            BL_frames,'freqPassBand',1,'sigmaXY', 2,...
            'sigmaT', 0.1,'thresholdPuff', 7, 'threshold2D', 0.2,...
            'minRiseTime',0.0845, 'maxRiseTime', 1,'minROIArea', 10,...
            'dilateXY', 5, 'dilateT', 0.3,'erodeXY', 1, 'erodeT', 0.1,...
            'discardBorderROIs',true);
        
        %         % NEURONS
        %         % 2D FLIKA selected for peaks from "dendrites"
        %         Neur_findConf{1} = ConfigFindROIsFLIKA_2D.from_preset('ca_neuron', 'baselineFrames',...
        %             BL_frames,'freqPassBand',1,'sigmaXY', 2,...
        %             'sigmaT', 0.1,'thresholdPuff', 7, 'threshold2D', 0.2,...
        %             'minRiseTime',0.0845, 'maxRiseTime', 2,'minROIArea', 10,...
        %             'dilateXY', 3, 'dilateT', 0.5,'erodeXY', 2, 'erodeT', 0.3,...
        %             'discardBorderROIs',true);
        %
        %         % hand selected for peaks from somata
        %         x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
        %         Neur_findConf{2} = ConfigFindROIsDummy.from_ImageJ(fullfile(testRoot,'Neurons.zip'), x_pix, y_pix,1);
        %
        %         % 3D FLIKA selected for time and space estimations
        %         Neur_findConf{3} = ConfigFindROIsFLIKA_3D.from_preset('ca_neuron', 'baselineFrames',...
        %             BL_frames,'freqPassBand',1,'sigmaXY', 2,...
        %             'sigmaT', 0.1,'thresholdPuff', 7, 'threshold2D', 0.2,...
        %             'minRiseTime',0.0845, 'maxRiseTime', 2,'minROIArea', 10,...
        %             'dilateXY', 3, 'dilateT', 0.5,'erodeXY', 2, 'erodeT', 0.3,...
        %             'discardBorderROIs',true);
        
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
        %         configCS{1,4} = ConfigCellScan(Neur_findConf{1,1}, measureConf,detectConf{1,2}); % neuronal FLIKA, peaks
        %         configCS{1,5} = ConfigCellScan(Neur_findConf{1,2}, measureConf,detectConf{1,2}); % neuronal hand peaks
        %         configCS{1,6} = ConfigCellScan(Neur_findConf{1,3}, measureConf,detectConf{1,3}); % neuronal 3D FLIKA
        %
        %% Create CellScan objects
        CSArray_Ch1_FLIKA = CellScan(fnList, ImgArray, configCS{1,1}, 1);
        
        CSArray_Ch1_Hand = CellScan(fnList, ImgArray, configCS{1,2}, 1);
        
        %         CSArray_Ch2_FLIKA = CellScan(fnList, ImgArray, configCS{1,4}, 2);
        
        %         CSArray_Ch2_Hand = CellScan(fnList, ImgArray, configCS{1,5}, 2);
        
        
        
        %% Process the images
        CSArray_Ch1_FLIKA =CSArray_Ch1_FLIKA.process();
        CSArray_Ch1_Hand =CSArray_Ch1_Hand.process();
        %         CSArray_Ch2_Hand =CSArray_Ch2_Hand.process();
        %         CSArray_Ch2_FLIKA =CSArray_Ch2_FLIKA.process();
        
        %             CSArray_Ch1_FLIKA.plot();
        %             CSArray_Ch2_FLIKA.plot();
        
        % CSArray_Ch1_FLIKA.opt_config();
        %% Output data
        
        % make a giant data table
        listFields = {'amplitude', 'fullWidth', 'halfWidth', ...
            'numPeaks', 'peakTime', 'peakStart', 'peakStartHalf', ...
            'peakType', 'prominence', 'roiName', 'peakAUC'};
        
        CellScans=vertcat(CSArray_Ch1_FLIKA, CSArray_Ch1_Hand)
        
        % loop through cellscans
        for iScan=1:size(CellScans,1)
            for itrial=1:size(CellScans,2)
                
                
                % peak output
                temp=CellScans(iScan,itrial).calcDetectSigs.data;
                temp2.trialname ={};
                temp2.animalname = {};
                temp2.channel = {};
                temp2.Spot = {};
                temp2.Cond = {};
                temp2.depth = {};
                temp2.area = {};
                temp2.pixelsize = {};
                
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
                    temp2.channel{iPeak,1}= 'RCaMP';
                    temp2.Spot{iPeak,1}= spotId;
                    temp2.animalname{iPeak,1}= CurrentAnimal;
                    temp2.Cond{iPeak,1} = CurrentCondition;
                    
                    temp2.depth{iPeak,1} = CurrentDepth(1,1);
                    temp2.pixelsize{iPeak,1} = CellScans(iScan,itrial).rawImg.metadata.pixelSize;
                    temp2.area{iPeak,1} = 0;
                end
                
                
                
                
                
                %%
                isFirst = (itrial == 1 && iScan == 1);
                if isFirst
                    data.Trial = {};
                    data.Animal = {};
                    data.Channel = {};
                    data.Spot = {};
                    data.Condition = {};
                    data.Depth = {};
                    data.area = {};
                    data.pixelsize={};
                end
                data.Trial= [data.Trial; temp2.trialname];
                data.Animal= [data.Animal; temp2.animalname];
                data.Channel= [data.Channel; temp2.channel];
                data.Spot= [data.Spot; temp2.Spot];
                data.Condition= [data.Condition; temp2.Cond];
                data.Depth= [data.Depth; temp2.depth];
                data.area= [data.area; temp2.area];
                data.pixelsize= [data.pixelsize; temp2.pixelsize];
                
                
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
                        Trace_data{iROI,3}= 'RCaMP';
                        
                        Trace_data{iROI,4}= spotId;
                        Trace_data{iROI,5}= CurrentAnimal;
                        Trace_data{iROI,6}= CurrentCondition;
                        Trace_data{iROI,7} = CurrentDepth(1,1);
                        %                        Trace_data{iROI,9} = CurrentBaseline(1,1);
                        Trace_data{iROI,8} = traces(:,iROI);
                        if iScan==1
                            Trace_data{iROI,9} = CellScans(iScan,itrial).calcFindROIs.data.roiIdxs{iROI,1};
                            Trace_data{iROI,10} = CellScans(iScan,itrial).rawImg.metadata.pixelSize;
                        else
                            Trace_data{iROI,9} = CellScans(iScan,itrial).calcFindROIs.data.roiMask(:,:,iROI);
                            Trace_data{iROI,10} = CellScans(iScan,itrial).rawImg.metadata.pixelSize;
                        end
                        
                        FrameRate= CellScans(1, 1).rawImg.metadata.frameRate;
                        Trace_data{iROI,11} = FrameRate; % frameRate
                        
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
                        Trace_data{iROI,12}= Onsets;
                        
                        % trace AUC
                        x2=round(FrameRate*(BL_time+1));
                        x3= round(FrameRate*(BL_time+10));
                        Trace_data{iROI,13}=trapz(traces(BL_frames:x2,iROI));
                        Trace_data{iROI,14}=trapz(traces(BL_frames:x3,iROI));
                        
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





% %% Save all data for R analysis
AllData2= [dataNames';AllData];

%onsetTimeTable
names={'ROI','Trial','Channel','Spot','Animal', 'Condition','depth',...
    'trace','ROIIdx','PixelSize','FrameRate','OnsetTime','TraceAUC1','TraceAUC10'};
All_traces2=vertcat(names, All_traces);
%All_traces2(:,10)=[];
%All_traces2(:,10)=[];

cd(fullfile(Settings.MainDir, 'Results'));
% write date to created file
cell2csv(SaveFiles{1,1}, AllData2);
save(SaveFiles{1,3}, 'All_traces','-v7.3');
save(SaveFiles{1,2}, 'AllData2','-v7.3');
cell2csv(SaveFiles{1,4}, All_traces2);


