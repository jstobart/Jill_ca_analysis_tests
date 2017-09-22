% analysis as of Sept. 2017

clearvars

tblRaw = table();
%% Information about your images
Settings.MainDir = 'E:\Data\Pericyte_project\Two-photon-data\Ret_ret_Mice\Baseline_BloodFlow';

Settings.AnimalNames = {
    'Prt7',...
    'Prt6',...
    'Prt5',...
    'Prt4',...
    'Prt3',...
    'Prt2',...
    'Prt1',...
    'Prr6',...
    'Prr5',...
    'Prr4',...
    'Prr3',...
    'Prr2',...
    'Nrr1',...
    'Prr7',...
    'Prr8',...
    'Prr9',...
    };
Settings.ScoreSheetNames = {
    'Prt7_DiamVelFilesScoresheet',...
    'Prt6_DiamVelFilesScoresheet',...
    'Prt5_DiamVelFilesScoresheet',...
    'Prt4_DiamVelFilesScoresheet',...
    'Prt3_DiamVelFilesScoresheet',...
    'Prt2_DiamVelFilesScoresheet',...
    'Prt1_DiamVelFilesScoresheet',...
    'Prr6_DiamVelFilesScoresheet',...
    'Prr5_DiamVelFilesScoresheet',...
    'Prr4_DiamVelFilesScoresheet',...
    'Prr3_DiamVelFilesScoresheet',...
    'Prr2_DiamVelFilesScoresheet',...
    'Nrr1_DiamVelFilesScoresheet',...
    'Prr7_DiamVelFilesScoresheet',...
    'Prr8_DiamVelFilesScoresheet',...
    'Prr9_DiamVelFilesScoresheet',...
    };

% final data file name
SaveFiles{1,1} = fullfile(Settings.MainDir, 'Results', 'Ret+_diam_velocity_22_09_2017.csv');%'Control_Peaks_3Conds.csv'); % all data
SaveFiles{1,2} = 'Velocity.mat';%velocity
SaveFiles{1,3} = 'Diameter.mat';%diameter

%% Load calibration file
calibration ='E:\matlab\CalibrationFiles\calibration_20x.mat';
CalFile = CalibrationPixelSize.load(calibration);


%% load scoresheet and loop through animal, vessel, etc.

Settings.ScoreSheetPath = fullfile(Settings.MainDir,'Scoresheets',Settings.ScoreSheetNames);

numAnimals = length(Settings.AnimalNames);
for iAnimal = 1:numAnimals
    CurrentAnimal = Settings.AnimalNames{iAnimal};
    disp(CurrentAnimal);
    %read Scoresheet of current animal
    CurrentSheet = Settings.ScoreSheetPath{iAnimal};
    Settings = readScoresheet2(CurrentSheet, Settings);
    
    % Get SpotID
    spots = unique(Settings.SpotIDs);
    
    for iSpot = 1:length(spots)
        
        % Extract spot name
        spotId = spots{iSpot};
        disp(spotId);
        % Find the idx of paths matching this spot (all Conditions)
        matchingSpotsIdx = find(~cellfun(@isempty, regexp(Settings.SpotIDs, spotId)));
        testRoot = Settings.LowresPath{matchingSpotsIdx};
        CurrentDepth = Settings.Depth(matchingSpotsIdx);
        CurrentBranch = Settings.BranchOrder(matchingSpotsIdx);
        channel = struct('blood_plasma',Settings.Channel(matchingSpotsIdx));
        
        if exist(fullfile(testRoot,SaveFiles{1,2}),'file')
            load(fullfile(testRoot,SaveFiles{1,2}))
        else
            
            % load velocity data
            Velocityexpfiles = dir(fullfile(testRoot,'v*'));
            V_fnTempList = {Velocityexpfiles(:).name};
            V_fnList = fullfile(testRoot, V_fnTempList);
            
            % Create an array of ScanImage Tiffs
            V_ImgArray =  SCIM_Tif(V_fnList, channel, CalFile);
            
            %% Configs for Flux calculations from velocity
            VelConfig = ConfigVelocityRadon ('thresholdProm',Settings.VelocityThresholdProm(matchingSpotsIdx),...
                'windowTime',Settings.VelocityWindowTime(matchingSpotsIdx),...
                'ThresholdSNR',Settings.ThresholdSNR(matchingSpotsIdx));
            
            ColstoUseVel=[32 96];
            Columns=ColstoUseVel*Settings.ScaleFactor(matchingSpotsIdx);
            
            %% Create LineScanVel objects
            
            Flux = LineScanVel(V_fnList, V_ImgArray,VelConfig, 1, Columns);
            
            Flux =Flux.process();
            %Flux(1).plot();
            Flux(1).opt_config();
            
            % SAVE OBJECT
            save('-v7.3', fullfile(testRoot,SaveFiles{1,2}), 'Flux')
        end
        
        %% load diameter data
        if ~exist(fullfile(testRoot,SaveFiles{1,3}),'file')
            
            Diameter_expfiles = dir(fullfile(testRoot,'d*'));
            D_fnTempList = {Diameter_expfiles(:).name};
            D_fnList = fullfile(testRoot, D_fnTempList);
            
            % Create an array of ScanImage Tiffs
            D_ImgArray =  SCIM_Tif(D_fnList, channel, CalFile);
            
            ColstoUseDiam=[2 122]; % 115
            Columns2=ColstoUseDiam*Settings.ScaleFactor(matchingSpotsIdx);
            
            % only use part of the data
            %          [D_ImgArray2,~]=split1(D_ImgArray,1,[7500 size(D_ImgArray.rawdata,1)-7500]);
            %         ;
            %% Create Diameter objects
            if Settings.ScaleFactor==1
                DiamConfig= ConfigDiameterFWHM('maxRate',5);
            else
                DiamConfig=ConfigDiameterFWHM('maxRate',10);
            end
            
            Diameter= LineScanDiam(D_fnList, D_ImgArray,DiamConfig,Columns2, Settings.Channel(matchingSpotsIdx),[]);
            Diameter.process();
            %Diameter(1).opt_config()
            Diameter(1).plot()
            
            save('-v7.3', fullfile(testRoot,SaveFiles{1,3}), 'Diameter')
        else
            load(fullfile(testRoot,SaveFiles{1,3}))
        end
        
        
        close all
        %% Reprocess Velocity to make sure all children are the same
        % this is needed to make sure the pulsatility is accurate
        
        %check that Flux is the same as the spreadsheet (for ALL
        % Children!)
        % reprocess it is not and output the data
        ThresholdProm=Settings.VelocityThresholdProm(matchingSpotsIdx);
        WindowTime=Settings.VelocityWindowTime(matchingSpotsIdx);
        ThresholdSNR=Settings.ThresholdSNR(matchingSpotsIdx);
        
        for iChildren = 1:length(Flux)
            if Flux(1,iChildren).calcVelocity.config.windowTime~=WindowTime || Flux(1,iChildren).calcVelocity.config.thresholdProm~=ThresholdProm || Flux(1,iChildren).calcVelocity.config.thresholdSNR~=ThresholdSNR
                Flux(1,iChildren).calcVelocity.config.windowTime=WindowTime;
                Flux(1,iChildren).calcVelocity.config.thresholdProm=ThresholdProm;
                Flux(1,iChildren).calcVelocity.config.thresholdSNR=ThresholdSNR;
                Flux(1,iChildren).process();
            end
        end
        save('-v7.3', fullfile(testRoot,SaveFiles{1,2}), 'Flux')
        
        
        
        for iChildren = 1:length(Diameter)
            if Settings.ScaleFactor==1
                if Diameter(1,iChildren).calcDiameter.config.maxRate~=5
                    Diameter(1,iChildren).calcDiameter.config.maxRate=5;
                    Diameter(1,iChildren).process();
                end
            end
            if Diameter(1,1).calcDiameter.config.lev50~=Diameter(1,3).calcDiameter.config.lev50
                Diameter(1,iChildren).calcDiameter.config.lev50=Diameter(1,1).calcDiameter.config.lev50;
                Diameter(1,iChildren).process();
            end
        end
        save('-v7.3', fullfile(testRoot,SaveFiles{1,3}), 'Diameter')
        
        
        
        
        % Prepare some temporary tables data
        tblTempRaw = table();
        MeanDiam=zeros(1,length(Diameter));
        SD_Diam=zeros(1,length(Diameter));
        MeanVel=zeros(1,length(Flux));
        SD_Vel=zeros(1,length(Flux));
        MeanLD=zeros(1,length(Flux));
        SD_LD=zeros(1,length(Flux));
        MeanFlux=zeros(1,length(Flux));
        SD_flux=zeros(1,length(Flux));
        PI=zeros(1,length(Flux));
        
        % Extract the basic parameters from the 3D FLIKA CellScan
        tblTempRaw.Genotype=Settings.Genotype(matchingSpotsIdx);
        tblTempRaw.AnimalName={CurrentAnimal};
        tblTempRaw.Spot={spotId};
        tblTempRaw.Depth=CurrentDepth;
        tblTempRaw.BranchOrder=CurrentBranch;
        tblTempRaw.NoFlow=Settings.NoFlow(matchingSpotsIdx);
        
        %mean diameter and velocity
        % calculate mean of multiple trials (if they exist)
        for iChildren = 1:length(Diameter)
            MeanDiam(iChildren)=Diameter(1,iChildren).calcDiameter.data.means.raw.diameter;
            SD_Diam(iChildren)=Diameter(1,iChildren).calcDiameter.data.stdevs.raw.diameter;
        end
        tblTempRaw.Diameter=mean(MeanDiam);
        tblTempRaw.SD_Diam=mean(SD_Diam);
        
        for iChildren = 1:length(Flux)
            MeanVel(iChildren)=Flux(1,iChildren).calcVelocity.data.means.raw.velocity;
            SD_Vel(iChildren)=Flux(1,iChildren).calcVelocity.data.stdevs.raw.velocity;
            
            if Settings.NoFlow(matchingSpotsIdx)
                MeanVel(iChildren)=NaN;
                SD_Vel(iChildren)=NaN;
            end
            if Settings.FluxGood(matchingSpotsIdx)
                MeanFlux(iChildren)=Flux(1,iChildren).calcVelocity.data.means.raw.flux;
                MeanLD(iChildren)=Flux(1,iChildren).calcVelocity.data.means.raw.linearDensity;
                
                SD_flux(iChildren)=Flux(1,iChildren).calcVelocity.data.stdevs.raw.flux;
                SD_LD(iChildren)=Flux(1,iChildren).calcVelocity.data.stdevs.raw.linearDensity;
            else
                MeanFlux(iChildren)=NaN;
                MeanLD(iChildren)=NaN;
                
                SD_flux(iChildren)=NaN;
                SD_LD(iChildren)=NaN;
            end
        end
        tblTempRaw.Velocity=abs(mean(MeanVel));
        tblTempRaw.SD_Vel=abs(mean(SD_Vel));
        
        tblTempRaw.Flux=mean(MeanFlux);
        tblTempRaw.SD_Flux=mean(SD_flux);
        
        tblTempRaw.linearDensity=mean(MeanLD);
        tblTempRaw.SD_LD=mean(SD_LD);
        
        % hematocrit (volume of RBC/ total blood volume)
        %Hct= (#RBCs x Volume)/ (Pi x (radius)^2 x length)  RBC volume = 45um3 (Hawkey et al. Haemtology 1991).
        %Volume_radius= RBC volume / Pi*r^2
        Volume_radius = (45/(pi*((mean(MeanDiam))/2)^2))*0.001; %convert to mm
        
        % hematocrit= linear density* RBC volume radius above
        tblTempRaw.Hematocrit= (mean(MeanLD))*Volume_radius;
        
        
        %% reprocess velocity with same parameters for ALL vessels for pulsatility index
        % output data
        
        AllVelocities=[];
        for iChildren = 1:length(Flux)
            % reprocess
            if Flux(1,iChildren).calcVelocity.config.windowTime~=1000
                Flux(1,iChildren).calcVelocity.config.windowTime=1000;
                Flux(1,iChildren).process();
            end
            AllVelocities=vertcat(AllVelocities,Flux(1,iChildren).calcVelocity.data.velocity);
        end
        % pulsatility index= (Vmax-Vmin)/Vmean (see Bouvy et al. 2014)
        %                 %use the max 95th percentile and min 5th percentile
        Vmax = prctile(abs(AllVelocities),95);
        Vmin = prctile(abs(AllVelocities),5);
                
        tblTempRaw.VDiff=Vmax-Vmin;
        tblTempRaw.PulsatilityIndex= (Vmax-Vmin)/mean(abs(AllVelocities));
        
        if Settings.NoFlow(matchingSpotsIdx)
            tblTempRaw.VDiff=NaN;
            tblTempRaw.PulsatilityIndex=NaN;
        end
        if ~Settings.NoFlow(matchingSpotsIdx) && mean(abs(AllVelocities))== Inf
            break
        end

        
        % Add the data from this image to the table
        if Settings.DataGood(matchingSpotsIdx) || Settings.NoFlow(matchingSpotsIdx)
            
            tblRaw = [tblRaw; tblTempRaw];
        end
        
    end
end

%% Save all data for R analysis
cd(fullfile(Settings.MainDir, 'Results'));

% Save the data in 2 csv files, and also as a .mat file
delim = '\t';
writetable(tblRaw, SaveFiles{1,1}, 'Delimiter', delim)  % data table
