function configCS= loadConfigCS(BL_frames)
% find ROIs config

% For frame rate = 11.84Hz
   % CYTO GCaMP
   AC_findConf{1} = ConfigFindROIsFLIKA_2D.from_preset('ca_cyto_astro', 'sigmaXY', 2,...
       'sigmaT', 0.5,'threshold_std', 5,...
       'min_rise_time',0.1689, 'max_rise_time', 8,...
       'dilateXY', 4,...
       'dilateT', 0.5,'erodeXY', 2);

% LCK GCaMP
AC_findConf{2} = ConfigFindROIsFLIKA_2D.from_preset('ca_memb_astro', 'baselineFrames',...
    BL_frames,'freqPassBand',1,'sigmaXY', 2,...
    'sigmaT', 0.1,'threshold_std', 7, 'threshold_2D', 0.2,...
    'min_rise_time',0.0845, 'max_rise_time', 1,'minPuffArea', 10,...
    'dilateXY', 5, 'dilateT', 0.3,'erodeXY', 1, 'erodeT', 0.1,...
    'discardBorderROIs',true);

% hand selected- peaks from cellular structures
x_pix= 128; y_pix= 128;
AC_findConf{3} = ConfigFindROIsDummy.from_ImageJ([], x_pix, y_pix,1);

   AC_findConf{4} = ConfigFindROIsFLIKA_3D.from_preset('ca_cyto_astro', 'sigmaXY', 2,...
       'sigmaT', 0.5,'threshold_std', 5,...
       'min_rise_time',0.1689, 'max_rise_time', 8,...
       'dilateXY', 4,...
       'dilateT', 0.5,'erodeXY', 2);

% 2D FLIKA selected for peaks from "dendrites"
Neur_findConf{1} = ConfigFindROIsFLIKA_2D.from_preset('ca_neuron', 'baselineFrames',...
    BL_frames,'freqPassBand',1,'sigmaXY', 2,...
    'sigmaT', 0.1,'threshold_std', 7, 'threshold_2D', 0.2,...
    'min_rise_time',0.0845, 'max_rise_time', 2,'minPuffArea', 10,...
    'dilateXY', 3, 'dilateT', 0.5,'erodeXY', 2, 'erodeT', 0.3,...
    'discardBorderROIs',true);

% hand selected for peaks from somata
Neur_findConf{2} = ConfigFindROIsDummy.from_ImageJ([], x_pix, y_pix,1);

Neur_findConf{3} = ConfigFindROIsFLIKA_3D.from_preset('ca_neuron', 'baselineFrames',...
    BL_frames,'freqPassBand',1,'sigmaXY', 2,...
    'sigmaT', 0.1,'threshold_std', 7, 'threshold_2D', 0.2,...
    'min_rise_time',0.0845, 'max_rise_time', 2,'minPuffArea', 10,...
    'dilateXY', 3, 'dilateT', 0.5,'erodeXY', 2, 'erodeT', 0.3,...
    'discardBorderROIs',true);


% AWAKE astrocyte cyto calcium
detectConf{1} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...'normMethod', 'z-score','zIters', 100,...
    'propagateNaNs', false, 'excludeNaNs', false, 'lpWindowTime', 3, 'spFilterOrder', 2,...
    'spPassBandMin',0.05, 'spPassBandMax', 0.2, 'thresholdSD_low', 3,'thresholdSD_band', 5);

% AWAKE neuron calcium
detectConf{2} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,... 'normMethod','z-score',... 'zIters', 10000,...
    'propagateNaNs', false,'excludeNaNs', false, 'lpWindowTime', 2, 'spFilterOrder', 2,...
    'spPassBandMin',0.1, 'spPassBandMax', 1, 'thresholdSD_low', 3,'thresholdSD_band', 5);


% for 3D FLIKA
detectConf{3} = ConfigDetectSigsDummy();

% for calculating AUC for each trace
measureConf = ConfigMeasureROIsDummy('baselineFrames', BL_frames);

%
% Combine the configs into a CellScan config
configCS{1,1} = ConfigCellScan(AC_findConf{1,1}, measureConf,detectConf{1,1}); % astrocyte cyto FLIKA, peaks
configCS{1,2} = ConfigCellScan(AC_findConf{1,2}, measureConf,detectConf{1,1}); % astrocyte membrane peaks
configCS{1,3} = ConfigCellScan(AC_findConf{1,3}, measureConf,detectConf{1,1}); % astrocyte hand
configCS{1,4} = ConfigCellScan(Neur_findConf{1,1}, measureConf,detectConf{1,2}); % neuronal FLIKA, peaks
configCS{1,5} = ConfigCellScan(Neur_findConf{1,2}, measureConf,detectConf{1,2}); % neuronal hand peaks
configCS{1,6} = ConfigCellScan(AC_findConf{1,4}, measureConf,detectConf{1,3}); % astrocyte hand
configCS{1,7} = ConfigCellScan(Neur_findConf{1,3}, measureConf,detectConf{1,3}); % astrocyte hand

end