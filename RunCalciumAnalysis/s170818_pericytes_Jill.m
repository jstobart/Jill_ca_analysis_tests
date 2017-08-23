close all; clear all;

AllData= [];
All_traces= [];
AllPropData=[];


tblRaw = table();
tblSummary = table();
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
        
        expfiles = dir(fullfile(testRoot,'*.tif'));
        fnTempList = {expfiles(:).name};
        fnList = fullfile(testRoot, fnTempList);
        
        % Create an array of ScanImage Tiffs
        ImgArray =  SCIM_Tif(fnList, channel, CalFile);
        
        %% Run preprocessing steps
        % Spectral Unmixing
        ImgArray= ImgArray.unmix_chs(false, [], cell2mat(eGFP_TexasRed_Matrix));
        
        %Motion correction
        channelToUseMC = 1; % which channel to use
        refImg = squeeze(mean(ImgArray(1,1).rawdata(:,:,channelToUseMC, 5:10),4));
        ImgArray=ImgArray.motion_correct( 'refImg', refImg,'ch', channelToUseMC,'minCorr', 0.4);
        
        %ImgArray.plot();
        
        % only use part of the data
        [test2, ~] = split1(ImgArray, 4, [200 size(ImgArray.rawdata, 4) - 200]);
        ImgArray=test2;
        
        
        %% Configs for Finding ROIs
        % Pericyte calcium
        
        % automated selection for peak detection
        findConf{1} = ConfigFindROIsFLIKA_2D.from_preset('ca_memb_astro','baselineFrames',1:15,...
            'freqPassBand',2,'sigmaXY', 1,...
            'sigmaT', 0.1,'thresholdPuff', 7, 'threshold2D', 0.1,...
            'minRiseTime',0.16, 'maxRiseTime', 4,'minROIArea', 1.6,... 'maxROIArea'
            'minROITime', 1,'dilateXY', 0, 'dilateT', 0,'erodeXY', 0, 'erodeT', 0,...
            'discardBorderROIs',true);
        
        findConf{2} = ConfigFindROIsFLIKA_2p5D.from_preset('ca_memb_astro','baselineFrames',1:15,...
            'freqPassBand',2,'sigmaXY', 1,...
            'sigmaT', 0.1,'thresholdPuff', 7, 'threshold2D', 0.1,...
            'minRiseTime',0.16, 'maxRiseTime', 4,'minROIArea', 1.6,...'maxROIArea'
            'minROITime', 0.425,'dilateXY', 0, 'dilateT', 0,'erodeXY', 0, 'erodeT', 0,...
            'discardBorderROIs',true);
        
        % automated selection for signal propagation
        findConf{3} = ConfigFindROIsFLIKA_3D.from_preset('ca_memb_astro','baselineFrames',1:15,...
            'freqPassBand',2,'sigmaXY', 1,...
            'sigmaT', 0.1,'thresholdPuff', 7,... 'threshold2D', 0.2,...
            'minRiseTime',0.16, 'maxRiseTime', 4,'minROIArea', 2.5,...'maxROIArea'
            'minROITime', 0.425,'dilateXY', 0, 'dilateT', 0,'erodeXY', 0, 'erodeT', 0,...
            'discardBorderROIs',true);
        
        % hand selected ROIs
        x_pix= size(ImgArray.rawdata,2); y_pix= size(ImgArray.rawdata,1);
        scaleF = 1;
        % zipPath = 'D:/....';
        zipPath= fullfile(testRoot,'RoiSet.zip');
        findConf{4} = ConfigFindROIsDummy.from_ImageJ(zipPath, x_pix, y_pix, scaleF);
        
        
        %% Configs for measuring ROIs
        % AWAKE astrocyte membrane calcium
        detectConf = ConfigDetectSigsClsfy('propagateNaNs', false, 'excludeNaNs', false,...
            'lpWindowTime', 5, 'spFilterOrder', 2,'spPassBandMin',0.025, 'spPassBandMax', 1,...
            'thresholdLP', 7,'thresholdSP', 3);
        
        % for calculating AUC for each trace
        measureConf = ConfigMeasureROIsDummy();
        
        % Combine the configs into a CellScan config
        configCS{1,1} = ConfigCellScan(findConf{1,1}, measureConf,detectConf); % pericyte FLIKA, 2D
        configCS{1,2} = ConfigCellScan(findConf{1,2}, measureConf,detectConf); % pericyte FLIKA 2.5D
        configCS{1,3} = ConfigCellScan(findConf{1,3}, measureConf,detectConf); % pericyte FLIKA 3D
        configCS{1,4} = ConfigCellScan(findConf{1,4}, measureConf,detectConf); % pericyte hand selected
        
        %% Create CellScan objects
        somata = CellScan(fnList, ImgArray, configCS{1,4}, 1); % peaks from hand clicked ROIs
        %FLIKA_2D = CellScan(fnList, ImgArray, configCS{1,1}, 1); % ROIs from 2D FLIKA
        %FLIKA_2p5D = CellScan(fnList, ImgArray, configCS{1,2}, 1); % ROIs from 2.5D FLIKA
        FLIKA_3D = CellScan(fnList, ImgArray, configCS{1,3}, 1); % ROIs from 3D FLIKA
        
        %% Process the images
        somata =somata.process();
        %          boundary =boundary.process();
        %     FLIKA_2D =FLIKA_2D.process();
        %FLIKA_2p5D =FLIKA_2p5D.process();
        FLIKA_3D =FLIKA_3D.process();
        
        
        %% Make the debugging plots
        somata(1,1).plot();
        
        %         FLIKA_2D(1,1).plot();
        %         FLIKA_2D(1,1).plot('video');
        %
        %         FLIKA_2p5D(1,1).plot();
        %         FLIKA_2p5D(1,1).plot('video');
        
        FLIKA_3D(1,1).plot();
        %         FLIKA_3D(1,1).plot('video');
        
        % for working out find peaks parameters
        %         %somata(1).opt_config()
        %          FLIKA_2D(1).opt_config()
        %          FLIKA_2p5D(1).opt_config()
        %         FLIKA_3D(1).opt_config()
        
        
        
        %testCS = CellScan('', ImgArray,
        %ConfigCellScan(ConfigFindROIsCellSort()))      % independent component
        %analysis
        %% Exclude 'noise' ROIs that are outside the pericyte then sort somata and processes
        
        % boundary ROI mask (this is the area to be analyzed- ROIs must be
        % inside this area)
        BoundaryMaskRaw = any(somata(1).calcFindROIs.data.roiMask, 3);
        BoundaryMask = imresize(BoundaryMaskRaw, size(ImgArray(1).rawdata(:,:,1,1)));
        
        % somata ROI mask
        nSomata=length(somata(1).calcFindROIs.data.roiNames);
        SomaMask = false([size(BoundaryMask(:,:)), 1]);
        for iSoma=1:nSomata
            maskROI = somata(1).calcFindROIs.data.roiMask(:,:,iSoma);
            % is it the border ROI based on name?
            if isempty(strfind(somata(1).calcFindROIs.data.roiNames{iSoma,1}, 'B')) % look for border ROI
                SomaMask=SomaMask + maskROI;
            end
        end
        
        nStacks = numel(ImgArray);
        isInside = cell(1, nStacks);
        for iStacks = 1:nStacks
            % Extract the complete mask
            roiMask = FLIKA_3D(iStacks).calcFindROIs.data.roiMask;
            cc = bwconncomp(roiMask);
            ll = labelmatrix(cc);
            
            % Create the artificial 3D masks for only the in bound ROIs (i.e. no
            % garbage noise ROIs
            nROIs = cc.NumObjects;
            isInside{iStacks} = false(nROIs, 1);
            isSoma{iStacks} = false(nROIs, 1);
            for iROI = 1:nROIs
                % Identify if this ROI is inside the area
                maskROI = any(ll == iROI, 3);
                maskOverlap = BoundaryMask & maskROI;
                isInside{iStacks}(iROI) = any(maskOverlap(:));
                % look for overlap between somata and processes
                maskOverlap2= SomaMask & maskROI;
                isSoma{iStacks}(iROI) = any(maskOverlap2(:));
                
                % Create a 2p5d mask to make it cleaner to measure the ROI
                mask3D_temp = false(size(roiMask));
                mask3D_temp(cc.PixelIdxList{iROI}) = true;
                mask2p5D(:,:,iROI) = any(mask3D_temp, 3);
            end
            
            % Create the configs for only the processes
            config2p5D = ConfigCellScan(ConfigFindROIsDummy(...
                'roiMask', mask2p5D, 'roiNames', ...
                FLIKA_3D(iStacks).calcFindROIs.data.roiNames), ...
                ConfigMeasureROIsDummy(), ...
                ConfigDetectSigsDummy());
            FLIKA_2p5D(iStacks) = CellScan('', ImgArray(iStacks), config2p5D);
        end
        FLIKA_2p5D = FLIKA_2p5D.process();
        % FLIKA_2p5D(1).plot()
        % FLIKA_2p5D(1).plot('traces')
        % FLIKA_2p5D(1).plot('video')
        % figure, imagesc(FLIKA_2p5D(1).calcMeasureROIs.data.time, ...
        %     1:size(FLIKA_2p5D(1).calcMeasureROIs.data.tracesNorm, 2),...
        %     FLIKA_2p5D(1).calcMeasureROIs.data.tracesNorm');
        
        
        
        %% Output data
        
        for iStack = 1:nStacks
            % Prepare some temporary tables for 3D FLIKA and Field of View Summary
            tblTempRaw = table();
            tblTempSummary = table();
            
            % Extract the basic parameters from the 3D FLIKA CellScan
            nROIs = numel(FLIKA_3D(iStack).calcFindROIs.data.roiNames);
            tblTempRaw.animalname=repmat({CurrentAnimal}, nROIs, 1);
            tblTempRaw.Img = repmat({FLIKA_3D(iStack).rawImg.name}, nROIs, 1);
            tblTempRaw.trialname=repmat({strcat('trial', num2str(iStack))}, nROIs, 1);
            tblTempRaw.Spot=repmat({spotId}, nROIs, 1);
            tblTempRaw.celltype=repmat(CurrentCell(1,1), nROIs, 1);
            tblTempRaw.depth=repmat(CurrentDepth(1,1), nROIs, 1);
            
            tblTempRaw.ROI = FLIKA_3D(iStack).calcFindROIs.data.roiNames;
            tblTempRaw.volume = FLIKA_3D(iStack).calcFindROIs.data.volume;
            tblTempRaw.area = FLIKA_3D(iStack).calcFindROIs.data.area;
            tblTempRaw.duration = FLIKA_3D(iStack).calcFindROIs.data.duration;
            tblTempRaw.onset = FLIKA_3D(iStack).calcFindROIs.data.onset;
            
            
            % Extract some more parameters from the 2.5D traces, but using the
            % timing information from the 3D traces
            amplitude = zeros(nROIs, 1);
            auc = amplitude;
            for jROI = 1:nROIs
                traceExists = ...
                    FLIKA_3D(iStack).calcMeasureROIs.data.tracesExist(:,jROI);
                traceExtract = FLIKA_2p5D(iStack).calcMeasureROIs.data.tracesNorm(...
                    traceExists, jROI);
                timeExtract = FLIKA_2p5D(iStack).calcMeasureROIs.data.time(traceExists);
                amplitude(jROI) = max(traceExtract); %
                auc(jROI) = trapz(timeExtract, traceExtract);
                
            end
            tblTempRaw.Max_amplitude = amplitude;
            tblTempRaw.auc = auc;
            
            % Specify whether the ROI is a process or soma
            tblTempRaw.is_soma = isSoma{iStack};
            
            % delete the rows where the ROIs are outside the field of interest
            tblTempRaw(~isInside{iStack},:)=[];
            
            % Add the data from this image to the table
            tblRaw = [tblRaw; tblTempRaw];
            
            % Extract some summary data
            nCats = 2;
            tblTempSummary.animalname=repmat({CurrentAnimal}, nCats, 1);
            tblTempSummary.Img = repmat({FLIKA_3D(iStack).rawImg.name}, nCats, 1);
            tblTempSummary.trialname=repmat({strcat('trial', num2str(iStack))}, nCats, 1);
            tblTempSummary.Spot=repmat({spotId}, nCats, 1);
            tblTempSummary.celltype=repmat(CurrentCell(1,1), nCats, 1);
            tblTempSummary.depth=repmat(CurrentDepth(1,1), nCats, 1);
            
            tblTempSummary.fov_area = repmat(...
                (FLIKA_3D(iStack).rawImg.metadata.nPixelsPerLine.* ...
                FLIKA_3D(iStack).rawImg.metadata.pixelSize).*2, nCats, 1);
            tblTempSummary.img_duration = repmat(...
                FLIKA_3D(iStack).rawImg.metadata.nFrames./ ...
                FLIKA_3D(iStack).rawImg.metadata.frameRate, nCats, 1);
            tblTempSummary.num_somas = repmat(nSomata-1, nCats, 1);
            tblTempSummary.is_soma = [true; false];
            tblTempSummary.num_signals = [sum(isSoma{iStack}); sum(~isSoma{iStack})];
            tblTempSummary.frequency = (tblTempSummary.num_signals)./ ...
                ((tblTempSummary.img_duration./60).* ...
                (tblTempSummary.fov_area./100));
            tblSummary = [tblSummary; tblTempSummary];
            
            
            
            %% extract traces from ROIs for correlations etc.
            traces= FLIKA_2p5D(1,iStack).calcMeasureROIs.data.tracesNorm;
            %preallocate
            Trace_data=cell(size(traces,2),1);
            for iROI = 1:size(traces,2)
                Trace_data{iROI,1}= FLIKA_2p5D(1,iStack).calcFindROIs.data.roiNames{iROI,1};
                Trace_data{iROI,2}= FLIKA_3D(iStack).rawImg.name;
                Trace_data{iROI,3}= strcat('trial', num2str(iStack));
                Trace_data{iROI,4}= spotId;
                Trace_data{iROI,5}= CurrentAnimal;
                Trace_data{iROI,6}= CurrentCell;
                Trace_data{iROI,7} = CurrentDepth(1,1);
                Trace_data{iROI,8} = traces(:,iROI);
                Trace_data{iROI,9} = FLIKA_2p5D(1,iStack).calcFindROIs.data.roiMask(:,:,iROI);
                Trace_data{iROI,10} = FLIKA_2p5D(1,iStack).rawImg.metadata.pixelSize;
            end
            All_traces=vertcat(All_traces, Trace_data);
            clearvars Trace_data
            clearvars temp temp2
        end
        
        
        %% Pericyte signal propagation
        
        for itrial=1:length(FLIKA_3D)
            temp3=FLIKA_3D(1,itrial).calcFindROIs.data;
            
            %% Distance Calculations
            % find the minimium distance between the edges of the 3D
            % ROIs and somata
            
            % somata mask
            for iSoma=1:nSomata
                if isempty(strfind(somata(1).calcFindROIs.data.roiNames{iSoma,1}, 'B')) % look for border ROI
                    
                    SomaROI = double(somata(1).calcFindROIs.data.roiMask(:,:,iSoma));
                end
                
                for jROI = 1:nROIs
                    traceExists = ...
                        FLIKA_3D(iStack).calcMeasureROIs.data.tracesExist(:,jROI);
                    MaskExtract = Mask3D_temp(:,:,traceExist);
                    for iBlob = 1:size(MaskExtract)
                        ROIMask= double(MaskExtract(:,:,iBlob);
                        
                        % combine masks together
                        Mask=SomaROI + ROIMask;
                        Mask=im2bw(Mask);
                        
                        %Pythaogrean theorem method
                        
                        % Define object boundaries
                        boundaries = bwboundaries(Mask);
                        numberOfBoundaries = size(boundaries, 1);
                        if numberOfBoundaries==1
                            minDis_toSoma = 0;
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
                            minDis_toSoma = (min(minDistance)*FLIKA_3D(iStack).rawImg.metadata.pixelSize);
                        end
                        
                        % make a distance vector
                        DistanceVec{iBlob}=minDis_toSoma;
                        
                        clear boundaries boundary1 boundary2 boundary1x boundary1y boundary2x boundary2y minDistance indexofMin
                    end
                end
            end
            
            
            % create fields for trial, animal, spot, condition, etc.
            for iROI = 1:length(temp3.distance)
                temp4.roiNames{iROI,1}=temp3.roiNames{iROI,1};
                temp4.trialname{iROI,1}=strcat('trial', num2str(itrial));
                temp4.Spot{iROI,1}= spotId;
                temp4.animalname{iROI,1}= CurrentAnimal;
                temp4.celltype{iROI,1}=CurrentCell;
                temp4.depth{iROI,1} = CurrentDepth(1,1);
                temp4.propRate{iROI,1} = temp3.distance(iROI,1)/temp3.duration(iROI,1); % rate of propagation= distance of propagation/duration (um/s)
                temp4.pixelsize{iROI,1} = somata(1,itrial).rawImg.metadata.pixelSize;
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
    
    
    
    
    %% Save the data
    % Set up the filenames / directory
    dirData = 'F:\Data\2PLSM-temp\optic nerve\2016-10-06-N4005';
    fnData = '2016-10-06-N4005';
    fnFull = fullfile(dirData, fnData);
    % Save the data in 2 csv files, and also as a .mat file
    delim = '\t';
    writetable(tblRaw, [fnFull, '_raw.csv'], 'Delimiter', delim)
    writetable(tblSummary, [fnFull, '_summary.csv'], 'Delimiter', delim)
    save('-v7.3', [fnFull, '.mat'], 'cs*', 'ri*', 'tbl*')
