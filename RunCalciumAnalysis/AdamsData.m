
%% Adam's data

clearvars

AllData=[];
AllDataFinal=[];
%% load data

channel = struct('cellular_signal',1,'Ca_Cyto_Astro',2);
BL_frames = 7;

% Get image paths
FileNames1 ={'G:\AdamData\Whisker-1102\Whisker-1102.xml',...
    'G:\AdamData\Whisker-1102\Whisker-1102.xml'};

testRoot ='G:\AdamData\Whisker-1102\';


% run a loop through multiple files

for iFile=1:length(FilesNames1)
    
    CurrentFile=FileNames1(1,iFile);
    
    %% Open file
    ImgArray =  BioFormats(CurrentFile,channel,[]);
    
    
    
    %% Run motion correction
    
    %             HighRes = BioFormats(RefImgName,channel,[]);  % use a high
    %             res image as the reference
    
    % Extract a reference image from the first 3 frames
    refImg = mean(ImgArray.rawdata(:,:,2,1:3),4);
    
    % plot reference image
    figure('Name','ref Image for motion correction'), imagesc(refImg), axis image, axis off, colormap(gray)
    
    ImgArray = ImgArray.motion_correct('refImg',refImg, 'ch',2,'maxShift', 10,'minCorr', 0.4);%,'doPlot',true);
    
    % plot the motion corrected movie
    ImgArray.plot();
    
    
    %% Configs for Finding ROIs
    %ASTROCYTES
    % 2D automated selection for peaks
    AC_findConf{1} = ConfigFindROIsFLIKA_2D.from_preset('ca_cyto_astro', 'baselineFrames',...
        BL_frames,'thresholdPuff', 10, 'threshold2D',0.3,...
        'minROIarea', 49,'maxROIarea', 200000,...
        'dilateXY', 10,'dilateT', 2.51,'erodeXY', 15,'erodeT', 2.51);
    
    % hand selected- peaks from cellular structures
    x_pix= 512; y_pix= 512;
    AC_findConf{2} = ConfigFindROIsDummy.from_ImageJ(fullfile(testRoot,'RoiSet.zip'), x_pix, y_pix,1);
    
    
    % Configuration for measuring ROIs
    detectConf{1} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...'normMethod', 'z-score','zIters', 100,...
        'propagateNaNs', false, 'excludeNaNs', false, 'lpWindowTime', 5, 'spFilterOrder', 2,...
        'spPassBandMin',0.025, 'spPassBandMax', 0.075, 'thresholdLP', 5,'thresholdSP', 7);
    
    % for calculating AUC for each trace
    measureConf = ConfigMeasureROIsDummy('baselineFrames', BL_frames);
    
    %%
    % Combine the configs into a CellScan config
    configCS{1,1} = ConfigCellScan(AC_findConf{1,1}, measureConf,detectConf{1,1}); % astrocyte FLIKA, peaks
    configCS{1,2} = ConfigCellScan(AC_findConf{1,2}, measureConf,detectConf{1,1}); % astrocyte hand, peaks
    
    %% Create CellScan objects
    Astrocyte_FLIKA = CellScan(CurrentFile, ImgArray, configCS{1,1}, 2);
    
    Astrocyte_Hand = CellScan(CurrentFile, ImgArray, configCS{1,2}, 2);
    
    
    %% Process the images
    Astrocyte_FLIKA =Astrocyte_FLIKA.process();
    Astrocyte_Hand =Astrocyte_Hand.process();
    
    
    %% Plots
    Astrocyte_FLIKA.plot();
    Astrocyte_Hand.plot();
    
    %% Make the debugging plots
    
    % for working out find peaks parameters
    %Astrocyte_FLIKA.opt_config()
    
    
    %% Output data
    
       
    % make a giant data table
    listFields = {'amplitude', 'fullWidth', 'halfWidth', ...
        'numPeaks', 'peakTime', 'peakStart', 'peakStartHalf', ...
        'peakType', 'prominence', 'roiName', 'peakAUC'};
        
        % peak output from FLIKA
        temp=Astrocyte_FLIKA.calcDetectSigs.data;
        temp2.animalname = {};
        temp2.channel = {};
        temp2.Spot = {};

        
        % extract fields from Class
        for jField = 1:numel(listFields)
            isFirst = (iFile == 1 );
            if isFirst
                data.(listFields{jField}) = {};
            end
            data.(listFields{jField}) = [data.(listFields{jField}); ...
                temp.(listFields{jField})];
        end
        
        % create fields for trial, animal, spot, condition, etc.
        for iPeak = 1:length(temp.amplitude)
            temp2.channel{iPeak,1}= 'GCaMP';
            temp2.Spot{iPeak,1}={};
            temp2.animalname{iPeak,1}= {};
                    
        isFirst = (iFile == 1 );
        if isFirst
            data.Animal = {};
            data.Channel = {};
            data.Spot = {};
        end
        data.Animal= [data.Animal; temp2.animalname];
        data.Channel= [data.Channel; temp2.channel];
        data.Spot= [data.Spot; temp2.Spot];
 
    clearvars temp temp2
    
    % Astrocyte Hand clicked ROIs
        
        temp=Astrocyte_Hand.calcDetectSigs.data;
        temp2.animalname = {};
        temp2.channel = {};
        temp2.Spot = {};
       
        % extract fields from Class
        for jField = 1:numel(listFields)
            data.(listFields{jField}) = [data.(listFields{jField}); ...
                temp.(listFields{jField})];
        end
        
        % create fields for trial, animal, spot, condition, etc.
        for iPeak = 1:length(temp.amplitude)
            temp2.channel{iPeak,1}= 'GCaMP';
            temp2.Spot{iPeak,1}= {};
            temp2.animalname{iPeak,1}= {};
        end
        
        data.Animal= [data.Animal; temp2.animalname];
        data.Channel= [data.Channel; temp2.channel];
        data.Spot= [data.Spot; temp2.Spot];
        
    
    dataNames=fieldnames(data);
    data2= struct2cell(data);
    data3= [data2{:}];
    
    AllData=vertcat(AllData, data3);
    
    clearvars data data3 temp temp2
    
end


% %% Save all data for R analysis
AllDataFinal= [dataNames';AllData];

cd(fullfile(testRoot));

% write date to created file
cell2csv(SaveFiles{1,1}, AllDataFinal);


