
%% Adam's data

clearvars



%% load data

channel = struct('cellular_signal',1,'Ca_Cyto_Astro',2);
BL_frames = 7;

% Get image paths
testRoot1 ='G:\AdamData\Whisker-1102\Whisker-1102.xml';

testRoot ='G:\AdamData\Whisker-1102\';

% Open file
ImgArray =  BioFormats(testRoot1,channel,[]);



%% Run motion correction

%             HighRes = BioFormats(RefImgName,channel,[]);  % use a high
%             res image as the reference

% Extract a reference image from the first 3 frames
refImg = mean(ImgArray.rawdata(:,:,2,1:3),4);

% plot reference image
figure, imagesc(refImg), axis image, axis off, colormap(gray)

ImgArray = ImgArray.motion_correct('refImg',refImg, 'ch',2,'maxShift', 10,'minCorr', 0.4);%,'doPlot',true);
% fill in bad data?
if plotMotion==1
    ImgArray.plot();
end

%% Configs for Finding ROIs
%ASTROCYTES
% 2D automated selection for peaks
AC_findConf{1} = ConfigFindROIsFLIKA_2D.from_preset('ca_cyto_astro', 'baselineFrames',...
    BL_frames,'thresholdPuff', 10, 'threshold2D',0.3,...
    'minROIarea', 49,'maxROIarea', 200000,...
    'dilateXY', 10,'dilateT', 2.51,'erodeXY', 15,'erodeT', 2.51);

% hand selected- peaks from cellular structures
%x_pix= 512; y_pix= 512;
%AC_findConf{2} = ConfigFindROIsDummy.from_ImageJ(fullfile(testRoot,'Astrocytes.zip'), x_pix, y_pix,1);



%% Configuration for measuring ROIs
detectConf{1} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...'normMethod', 'z-score','zIters', 100,...
    'propagateNaNs', false, 'excludeNaNs', false, 'lpWindowTime', 5, 'spFilterOrder', 2,...
    'spPassBandMin',0.025, 'spPassBandMax', 0.075, 'thresholdLP', 5,'thresholdSP', 7);


% for calculating AUC for each trace
measureConf = ConfigMeasureROIsDummy('baselineFrames', BL_frames);

%%
% Combine the configs into a CellScan config
configCS{1,1} = ConfigCellScan(AC_findConf{1,1}, measureConf,detectConf{1,1}); % astrocyte FLIKA, peaks
%configCS{1,2} = ConfigCellScan(AC_findConf{1,2}, measureConf,detectConf{1,1}); % astrocyte hand, peaks

%% Create CellScan objects
Astrocyte_FLIKA = CellScan([], [], configCS{1,1}, 2);

%Astrocyte_Hand = CellScan([], [], configCS{1,2}, 1);


%% Process the images
Astrocyte_FLIKA =Astrocyte_FLIKA.process();
% Astrocyte_Hand =Astrocyte_Hand.process();


%% Make the debugging plots

% for working out find peaks parameters
Astrocyte_FLIKA.opt_config()



