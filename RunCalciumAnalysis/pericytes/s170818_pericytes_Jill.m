close all; clear all;

All_traces= [];
tblRaw = table();
tblSummary = table();

%% Information about your images
Settings.MainDir = 'E:\Data\Pericyte_project\Two-photon-data\Calcium';

Settings.ScoreSheetNames = {
    'Control_CalciumFilesScoresheet.xlsx',...
    };

channel = struct('Ca_Memb_Astro', 1, 'blood_plasma', 2);

% final data file name
SaveFiles{1,1} = fullfile(Settings.MainDir, 'Results', 'ROI_table_18_08_2017.csv');% all ROI data
SaveFiles{1,2} = fullfile(Settings.MainDir, 'Results', 'Summary_table_18_08_2017.csv');% summary for field of view
SaveFiles{1,3} = fullfile(Settings.MainDir, 'Results', 'test_traces_18_08_2017.mat');% normalized 2.5D traces
SaveFiles{1,4} = fullfile(Settings.MainDir, 'Results', 'somataROIs_18_08_2017.mat');% CellScans, etc.
SaveFiles{1,5} = fullfile(Settings.MainDir, 'Results', '3DROIs_18_08_2017.mat');% CellScans, etc.

BorderROIName = 'B';

doPlot=1;

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
        
        if exist('file',SaveFiles{1,4})
            load(SaveFiles{1,4})
        else
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
        channelToUseMC = Settings.MotionCorrChannel(iSpot); % which channel to use
        refImg = squeeze(mean(ImgArray(1,1).rawdata(:,:,channelToUseMC, 5:10),4));
        ImgArray=ImgArray.motion_correct( 'refImg', refImg,'ch', channelToUseMC,'minCorr', 0.4);
        
        %ImgArray.plot();
        
        % only use part of the data
%         [test2, ~] = split1(ImgArray(1,1), 4, [200 size(ImgArray(1,1).rawdata, 4) - 200]);
%         ImgArray=test2;
        
        
        %% Configs for Finding ROIs
        % Pericyte calcium
        
        % automated selection of ROIs
        findConf{1} = ConfigFindROIsFLIKA_3D.from_preset('ca_memb_astro','baselineFrames',1:15,...
            'freqPassBand',2,'sigmaXY', 1,...
            'sigmaT', 0.1,'thresholdPuff', 7,... 'threshold2D', 0.2,...
            'minRiseTime',0.16, 'maxRiseTime', 4,'minROIArea', 2.5,...'maxROIArea'
            'minROITime', 0.425,'dilateXY', 0, 'dilateT', 0,'erodeXY', 0, 'erodeT', 0,...
            'discardBorderROIs',true);
        
        % hand selected ROIs
        x_pix= size(ImgArray(1,1).rawdata,2); y_pix= size(ImgArray(1,1).rawdata,1);
        scaleF = 1;
        % zipPath = 'D:/....';
        zipPath= fullfile(testRoot,'RoiSet.zip');
        findConf{2} = ConfigFindROIsDummy.from_ImageJ(zipPath, x_pix, y_pix, scaleF);
        
        
        %% Configs for measuring ROIs
        % AWAKE astrocyte membrane calcium
        detectConf = ConfigDetectSigsClsfy('propagateNaNs', false, 'excludeNaNs', false,...
            'lpWindowTime', 5, 'spFilterOrder', 2,'spPassBandMin',0.025, 'spPassBandMax', 1,...
            'thresholdLP', 7,'thresholdSP', 3);
        
        % for calculating AUC for each trace
        measureConf = ConfigMeasureROIsDummy();
        
        % Combine the configs into a CellScan config
        configCS{1,1} = ConfigCellScan(findConf{1,1}, measureConf,detectConf); % pericyte FLIKA, 3D
        configCS{1,2} = ConfigCellScan(findConf{1,2}, measureConf,detectConf); % pericyte hand selected
        
        %% Create CellScan objects
        somata = CellScan(fnList, ImgArray, configCS{1,2}, 1); % peaks from hand clicked ROIs
        FLIKA_3D = CellScan(fnList, ImgArray, configCS{1,1}, 1); % ROIs from 3D FLIKA
        
        %% Process the images
        somata =somata.process();
        FLIKA_3D =FLIKA_3D.process();
       
        
        %% Make the debugging plots
        if doPlot
            somata(1,1).plot();
            
            FLIKA_3D(1,1).plot();
            %         FLIKA_3D(1,1).plot('video');
            
            % for working out find peaks parameters
            %         %somata(iStack).opt_config()
            %         FLIKA_3D(1).opt_config()
        end
        
        %% Exclude 'noise' ROIs that are outside the pericyte then sort somata and processes
        
        nStacks = numel(ImgArray);
        isInside = cell(1, nStacks);
        for iStacks = 1:nStacks
            % boundary ROI mask (this is the area to be analyzed- ROIs must be
            % inside this area)
            BoundaryMaskRaw = any(somata(iStacks).calcFindROIs.data.roiMask, 3);
            BoundaryMask = imresize(BoundaryMaskRaw, size(ImgArray(1).rawdata(:,:,1,1)));
            
            % somata ROI mask
            nSomata=length(somata(iStacks).calcFindROIs.data.roiNames);
            AllSomaMask = false([size(BoundaryMask(:,:)), 1]);
            for iSoma=1:nSomata
                maskROI = somata(iStacks).calcFindROIs.data.roiMask(:,:,iSoma);
                % is it the border ROI based on name?
                if isempty(strfind(somata(iStacks).calcFindROIs.data.roiNames{iSoma,1}, BorderROIName)) % look for border ROI
                    AllSomaMask=AllSomaMask + maskROI;
                end
            end
            
            
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
                maskOverlap2= AllSomaMask & maskROI;
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
        
        for iStacks = 1:nStacks
            % Prepare some temporary tables for 3D FLIKA and Field of View Summary
            tblTempRaw = table();
            tblTempSummary = table();
            
            % Extract the basic parameters from the 3D FLIKA CellScan
            nROIs = numel(FLIKA_3D(iStacks).calcFindROIs.data.roiNames);
            tblTempRaw.animalname=repmat({CurrentAnimal}, nROIs, 1);
            tblTempRaw.Img = repmat({FLIKA_3D(iStacks).rawImg.name}, nROIs, 1);
            tblTempRaw.trialname=repmat({strcat('trial', num2str(iStacks))}, nROIs, 1);
            tblTempRaw.Spot=repmat({spotId}, nROIs, 1);
            tblTempRaw.Drug=repmat(CurrentDrug(1,1), nROIs, 1);
            tblTempRaw.celltype=repmat(CurrentCell(1,1), nROIs, 1);
            tblTempRaw.depth=repmat(CurrentDepth(1,1), nROIs, 1);
            
            tblTempRaw.ROI = FLIKA_3D(iStacks).calcFindROIs.data.roiNames;
            tblTempRaw.volume = FLIKA_3D(iStacks).calcFindROIs.data.volume;
            tblTempRaw.area = FLIKA_3D(iStacks).calcFindROIs.data.area;
            tblTempRaw.duration = FLIKA_3D(iStacks).calcFindROIs.data.duration;
            tblTempRaw.onset = FLIKA_3D(iStacks).calcFindROIs.data.onset;
            
            
            % Extract some more parameters from the 2.5D traces, but using the
            % timing information from the 3D traces
            amplitude = zeros(nROIs, 1);
            auc = amplitude;
            for jROI = 1:nROIs
                traceExists = ...
                    FLIKA_3D(iStacks).calcMeasureROIs.data.tracesExist(:,jROI);
                traceExtract = FLIKA_2p5D(iStacks).calcMeasureROIs.data.tracesNorm(...
                    traceExists, jROI);
                timeExtract = FLIKA_2p5D(iStacks).calcMeasureROIs.data.time(traceExists);
                amplitude(jROI) = max(traceExtract); %
                auc(jROI) = trapz(timeExtract, traceExtract);
                
            end
            tblTempRaw.Max_amplitude = amplitude;
            tblTempRaw.auc = auc;
            
            %distance calculations
            tblTempRaw.centroidDis_Traveled = FLIKA_3D(iStacks).calcFindROIs.data.distance; %the ROI centroid has moved
            % propagation rate of the centroid
            tblTempRaw.centroidProp_Rate = (FLIKA_3D(iStacks).calcFindROIs.data.distance)./ ...
                FLIKA_3D(iStacks).calcFindROIs.data.duration;

            % minimum distance between 3D ROI edge to soma (hand-clicked)
            
            % 2.5D FLIKA mask
            Dis_to_soma = zeros(nROIs, 1);
            for jROI = 1:nROIs
                FLIKA_ROI= double(FLIKA_2p5D(iStacks).calcFindROIs.data.roiMask(:,:,jROI));
                
                % somata mask
                for iSoma=1:nSomata
                    
                    SomaMask = double(somata(iStacks).calcFindROIs.data.roiMask(:,:,iSoma));
                    somaDistance(iSoma)=minDistance(FLIKA_ROI,SomaMask); % # of pixels
                    
                    
                    % exclude border ROI
                    if ~isempty(strfind(somata(iStacks).calcFindROIs.data.roiNames{iSoma,1}, BorderROIName))
                        somaDistance(iSoma) = NaN;
                    end
                end
                Dis_to_soma(jROI,1) = min(somaDistance)*FLIKA_3D(iStacks).rawImg.metadata.pixelSize; % in microns
            end
            tblTempRaw.Dis_to_Soma = Dis_to_soma;
            
            
            % Specify whether the ROI is a process or soma
            tblTempRaw.is_soma = isSoma{iStacks};
            
            % delete the rows where the ROIs are outside the field of interest
            tblTempRaw(~isInside{iStacks},:)=[];
            
            % Add the data from this image to the table
            tblRaw = [tblRaw; tblTempRaw];
            
            % Extract some summary data
            nCats = 2;
            tblTempSummary.animalname=repmat({CurrentAnimal}, nCats, 1);
            tblTempSummary.Img = repmat({FLIKA_3D(iStacks).rawImg.name}, nCats, 1);
            tblTempSummary.trialname=repmat({strcat('trial', num2str(iStacks))}, nCats, 1);
            tblTempSummary.Spot=repmat({spotId}, nCats, 1);
            tblTempSummary.celltype=repmat(CurrentCell(1,1), nCats, 1);
            tblTempSummary.Drug=repmat(CurrentDrug(1,1), nCats, 1);
            tblTempSummary.depth=repmat(CurrentDepth(1,1), nCats, 1);
            
            tblTempSummary.fov_area = repmat(...
                (FLIKA_3D(iStacks).rawImg.metadata.nPixelsPerLine.* ...
                FLIKA_3D(iStacks).rawImg.metadata.pixelSize).*2, nCats, 1);
            tblTempSummary.img_duration = repmat(...
                FLIKA_3D(iStacks).rawImg.metadata.nFrames./ ...
                FLIKA_3D(iStacks).rawImg.metadata.frameRate, nCats, 1);
            tblTempSummary.num_somas = repmat(nSomata-1, nCats, 1);
            tblTempSummary.is_soma = [true; false];
            tblTempSummary.num_signals = [sum(isSoma{iStacks}); sum(~isSoma{iStacks})];
            tblTempSummary.frequency = (tblTempSummary.num_signals)./ ...
                ((tblTempSummary.img_duration./60).* ...
                (tblTempSummary.num_somas./100));
            tblSummary = [tblSummary; tblTempSummary];
            
            
            
            %% extract traces and masks from ROIs.
            traces= FLIKA_2p5D(1,iStacks).calcMeasureROIs.data.tracesNorm;

            %preallocate
            Trace_data=cell(size(traces,2),1);
            for iROI = 1:size(traces,2)
                Trace_data{iROI,1}= FLIKA_2p5D(1,iStacks).calcFindROIs.data.roiNames{iROI,1};
                Trace_data{iROI,2}= FLIKA_3D(iStacks).rawImg.name;
                Trace_data{iROI,3}= strcat('trial', num2str(iStacks));
                Trace_data{iROI,4}= spotId;
                Trace_data{iROI,5}= CurrentAnimal;
                Trace_data{iROI,6}= CurrentCell;
                Trace_data{iROI,7}= CurrentDrug;
                Trace_data{iROI,8} = CurrentDepth(1,1);
                Trace_data{iROI,9} = isSoma{iStacks}(iROI);
                Trace_data{iROI,10} = traces(:,iROI);
                Trace_data{iROI,11} = FLIKA_2p5D(1,iStacks).calcFindROIs.data.roiMask(:,:,iROI);
                Trace_data{iROI,12} = FLIKA_2p5D(1,iStacks).rawImg.metadata.pixelSize;
            end
            All_traces=vertcat(All_traces, Trace_data);
        end
        
        
    end
    
end

%% Save all data for R analysis
cd(fullfile(Settings.MainDir, 'Results'));

% Save the data in 2 csv files, and also as a .mat file
delim = '\t';
writetable(tblRaw, SaveFiles{1,1}, 'Delimiter', delim)  % ROI table
writetable(tblSummary, SaveFiles{1,2}, 'Delimiter', delim) % summary table
save(SaveFiles{1,3}, 'All_traces','-v7.3');  % all traces
save('-v7.3', SaveFiles{1,4}, 'somata') % cell scans etc.
save('-v7.3', SaveFiles{1,5}, '3DFLIKA') % cell scans etc.

