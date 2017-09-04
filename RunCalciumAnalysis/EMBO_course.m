
%% make a mixing matrix for spectral unmixing of RCaMP and GCaMP

mixedImg = SCIM_Tif();

[unmixedImg, RCaMP_GCaMP_Matrix] = unmix_chs(mixedImg); %returns the mixing matrix to be used for all imaging



%%

%load data images
ImgArray =  SCIM_Tif();


%spectral unmix
ImgArray= ImgArray.unmix_chs(false, [], cell2mat(RCaMP_GCaMP_Matrix));
%ImgArray.plot()

%% motion correct 
refImg = squeeze(mean(ImgArray(1,1).rawdata(:,:,1, 1:5),4));

%2D convolution
ImgArray_MC1=ImgArray.motion_correct( 'refImg', refImg,'ch', 1,'minCorr', 0.4);
%ImgArray_MC1.plot()
% 2D motion corrected movie
[img1ch_MC1, ~] = ImgArray_MC1.split1(3, [1, 1]);
img1ch_MC1.plot('CAxis', [0, 3000])


%HMM
% ImgArray2 =  SCIM_Tif();
% refImg2 = squeeze(mean(ImgArray2(1,1).rawdata(:,:,2, 1:5),4));
% ImgArray_MC3=ImgArray2.motion_correct( 'refImg', refImg2,'ch', 2, 'method', 'hmm');
% 
% % HMM motion corrected movie
% [img1ch_MC3, img2ch_MC3] = ImgArray_MC3.split1(3, [1, 1]);
% img2ch_MC3.plot
% 
% 
% ImgArray_MC2=ImgArray.motion_correct( 'refImg', refImg,'ch', 1, 'method', 'hmm');
% 
% % HMM motion corrected movie
% [img1ch_MC2, ~] = ImgArray_MC2.split1(3, [1, 1]);
% img1ch_MC2.plot('CAxis', [0, 3000])

% % denoising with block matching 3D
 %ImgArray.denoise();
% 
% [img1ch, ~] = ImgArray.split1(3, [1, 1]);
% img1ch.plot


FLIKA3D = CellScan([], ImgArray, [], 1);
FLIKA3D.process();
FLIKA3D.plot();
FLIKA3D.plot('video');

FLIKA3D.opt_config()


FLIKA2D = CellScan([], ImgArray, [], 1);
FLIKA2D.process();
FLIKA2D.plot('plotROIs',[6,7,8,9,10,12,13,14,15,16,17,18,19], 'spacingFactor',2)
FLIKA2D.opt_config();

%figure, FLIKA2D.calcMeasureROIs.plot(FLIKA2D, 'plotROIs',[6,7,8,9,10,12,13,14,15,16,17,18,19], 'spacingFactor',2)

Hand = CellScan([], ImgArray, [], 1);
Hand.process();
Hand.plot()%('ROIs',[5:10]);

figure, Hand.calcMeasureROIs.plot(Hand)



%%
% BL_frames=50;
% % find ROIs config
% 
% % For frame rate = 11.84Hz
%    % CYTO GCaMP
%    AC_findConf{1} = ConfigFindROIsFLIKA_2D.from_preset('ca_cyto_astro', 'sigmaXY', 2,...
%        'sigmaT', 0.5,'threshold_std', 5,...
%        'min_rise_time',0.1689, 'max_rise_time', 8,...
%        'dilateXY', 4,...
%        'dilateT', 0.5,'erodeXY', 2);
% 
% % LCK GCaMP
% AC_findConf{2} = ConfigFindROIsFLIKA_2D.from_preset('ca_memb_astro', 'baselineFrames',...
%     BL_frames,'freqPassBand',1,'sigmaXY', 2,...
%     'sigmaT', 0.1,'threshold_std', 7, 'threshold_2D', 0.2,...
%     'min_rise_time',0.0845, 'max_rise_time', 1,'minPuffArea', 10,...
%     'dilateXY', 5, 'dilateT', 0.3,'erodeXY', 1, 'erodeT', 0.1,...
%     'discardBorderROIs',true);
% 
% % hand selected- peaks from cellular structures
% x_pix= 128; y_pix= 128;
% AC_findConf{3} = ConfigFindROIsDummy.from_ImageJ([], x_pix, y_pix,1);
% 
%    AC_findConf{4} = ConfigFindROIsFLIKA_3D.from_preset('ca_cyto_astro', 'sigmaXY', 2,...
%        'sigmaT', 0.5,'threshold_std', 5,...
%        'min_rise_time',0.1689, 'max_rise_time', 8,...
%        'dilateXY', 4,...
%        'dilateT', 0.5,'erodeXY', 2);
% 
% % 2D FLIKA selected for peaks from "dendrites"
% Neur_findConf{1} = ConfigFindROIsFLIKA_2D.from_preset('ca_neuron', 'baselineFrames',...
%     BL_frames,'freqPassBand',1,'sigmaXY', 2,...
%     'sigmaT', 0.1,'threshold_std', 7, 'threshold_2D', 0.2,...
%     'min_rise_time',0.0845, 'max_rise_time', 2,'minPuffArea', 10,...
%     'dilateXY', 3, 'dilateT', 0.5,'erodeXY', 2, 'erodeT', 0.3,...
%     'discardBorderROIs',true);
% 
% % hand selected for peaks from somata
% Neur_findConf{2} = ConfigFindROIsDummy.from_ImageJ([], x_pix, y_pix,1);
% 
% Neur_findConf{3} = ConfigFindROIsFLIKA_3D.from_preset('ca_neuron', 'baselineFrames',...
%     BL_frames,'freqPassBand',1,'sigmaXY', 2,...
%     'sigmaT', 0.1,'threshold_std', 7, 'threshold_2D', 0.2,...
%     'min_rise_time',0.0845, 'max_rise_time', 2,'minPuffArea', 10,...
%     'dilateXY', 3, 'dilateT', 0.5,'erodeXY', 2, 'erodeT', 0.3,...
%     'discardBorderROIs',true);
% 
% 
% % AWAKE astrocyte cyto calcium
% detectConf{1} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...'normMethod', 'z-score','zIters', 100,...
%     'propagateNaNs', false, 'excludeNaNs', false, 'lpWindowTime', 3, 'spFilterOrder', 2,...
%     'spPassBandMin',0.05, 'spPassBandMax', 0.2, 'thresholdSD_low', 3,'thresholdSD_band', 5);
% 
% % AWAKE neuron calcium
% detectConf{2} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,... 'normMethod','z-score',... 'zIters', 10000,...
%     'propagateNaNs', false,'excludeNaNs', false, 'lpWindowTime', 2, 'spFilterOrder', 2,...
%     'spPassBandMin',0.1, 'spPassBandMax', 1, 'thresholdSD_low', 3,'thresholdSD_band', 5);
% 
% 
% % for 3D FLIKA
% detectConf{3} = ConfigDetectSigsDummy();
% 
% % for calculating AUC for each trace
% measureConf = ConfigMeasureROIsDummy('baselineFrames', BL_frames);
% 
% %
% % Combine the configs into a CellScan config
% configCS{1,1} = ConfigCellScan(AC_findConf{1,1}, measureConf,detectConf{1,1}); % astrocyte cyto FLIKA, peaks
% configCS{1,2} = ConfigCellScan(AC_findConf{1,2}, measureConf,detectConf{1,1}); % astrocyte membrane peaks
% configCS{1,3} = ConfigCellScan(AC_findConf{1,3}, measureConf,detectConf{1,1}); % astrocyte hand
% configCS{1,4} = ConfigCellScan(Neur_findConf{1,1}, measureConf,detectConf{1,2}); % neuronal FLIKA, peaks
% configCS{1,5} = ConfigCellScan(Neur_findConf{1,2}, measureConf,detectConf{1,2}); % neuronal hand peaks
% configCS{1,6} = ConfigCellScan(AC_findConf{1,4}, measureConf,detectConf{1,3}); % astrocyte hand
% configCS{1,7} = ConfigCellScan(Neur_findConf{1,3}, measureConf,detectConf{1,3}); % astrocyte hand
% 
% 
% %%
% % FLIKA
% FLIKA_2D = CellScan([], ImgArray, configCS{1,1}, 1);  % astrocyte 2D FLIKA
% FLIKA_2D.process();
% FLIKA_2D.plot();
% 
% FLIKA_2DN = CellScan([], ImgArray, configCS{1,4}, 2);  % neuronal 2D FLIKA
% FLIKA_2DN.process();
% FLIKA_2DN.plot();
% 
% 
% Hand_N = CellScan([], ImgArray, configCS{1,5}, 2);  % neuronal hand selected
% Hand_N.process();
% Hand_N.plot();
% 
% 
% FLIKA3D = CellScan([], ImgArray, configCS{1,6}, 1);
% FLIKA3D.process();
% FLIKA3D.plot();
% FLIKA3D.plot('video');
% 
% 
% 
% FLIKA3D.opt_config()

%             %% Create CellScan objects
%             CSArray_Ch1_FLIKA = CellScan(fnList, ImgArray, configCS{1,1}, 1);
%
%             CSArray_Ch1_Hand = CellScan(fnList, ImgArray, configCS{1,2}, 1);
%
%             CSArray_Ch2_FLIKA = CellScan(fnList, ImgArray, configCS{1,4}, 2);
%
%             CSArray_Ch2_Hand = CellScan(fnList, ImgArray, configCS{1,5}, 2);
%
%
%
%             %% Process the images
%             CSArray_Ch1_FLIKA =CSArray_Ch1_FLIKA.process();
%             CSArray_Ch1_Hand =CSArray_Ch1_Hand.process();
%             CSArray_Ch2_Hand =CSArray_Ch2_Hand.process();
%             CSArray_Ch2_FLIKA =CSArray_Ch2_FLIKA.process();
%
%
%             %% Output data
%
%             % make a giant data table
%             listFields = {'amplitude', 'fullWidth', 'halfWidth', ...
%                 'numPeaks', 'peakTime', 'peakStart', 'peakStartHalf', ...
%                 'peakType', 'prominence', 'ROIname', 'peakAUC'};
%
%             % Astrocyte FLIKA
%             for itrial=1:length(CSArray_Ch1_FLIKA)
%
%                 % peak output
%                 temp=CSArray_Ch1_FLIKA(1,itrial).calcDetectSigs.data;
%                 temp2.trialname ={};
%                 temp2.animalname = {};
%                 temp2.channel = {};
%                 temp2.Spot = {};
%                 temp2.Cond = {};
%                 temp2.depth = {};
%                 temp2.overlap = {};
%                 temp2.area = {};
%                 temp2.pixelsize = {};
%
%                 % extract fields from Class
%                 for jField = 1:numel(listFields)
%                     isFirst = (itrial == 1 );
%                     if isFirst
%                         data.(listFields{jField}) = {};
%                     end
%                     data.(listFields{jField}) = [data.(listFields{jField}); ...
%                         temp.(listFields{jField})];
%                 end
%
%                 % create fields for trial, animal, spot, condition, etc.
%                 for iPeak = 1:length(temp.amplitude)
%                     temp2.trialname{iPeak,1}=strcat('trial', num2str(itrial,'%02d'));
%                     temp2.channel{iPeak,1}= 'GCaMP';
%                     temp2.Spot{iPeak,1}= spotId;
%                     temp2.animalname{iPeak,1}= CurrentAnimal;
%                     temp2.Cond{iPeak,1} = CurrentCondition;
%                     temp2.depth{iPeak,1} = CurrentDepth(1,1);
%                     temp2.pixelsize{iPeak,1} = CSArray_Ch1_FLIKA(1,itrial).rawImg.metadata.pixelSize;
%
%                     % get the indices  and area for a particular ROI
%                     jROIname = temp.ROIname{iPeak};
%                     ROIindex= strcmp(CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.roiNames,jROIname);
%                     temp2.area{iPeak,1} = CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.area(ROIindex);
%
%                     ROIpuffIdx = CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.puffIdxs(ROIindex);
%                     x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
%                     [jy,jx,~] = ind2sub([y_pix,x_pix],cell2mat(ROIpuffIdx));
%                     jx = round(mean(jx));
%                     jy = round(mean(jy));
%                     xclose = round(x_pix*0.03); % 3 percent of pixels
%                     yclose = round(y_pix*0.03); % 3 percent of pixels
%
%                     % find astrocyte process and soma ROIs that have
%                     % similar centroids
%                     for kROI= 1:length(CSArray_Ch1_Hand(1,itrial).calcFindROIs.data.roiNames)
%                         % get the indices of the handclicked ROIs
%                         ROIMask{kROI}= CSArray_Ch1_Hand(1, 1).calcFindROIs.data.roiMask(:,:,kROI);
%                         ROI_Idx{kROI} = find(ROIMask{kROI});
%                         % mean pixels for soma ROI
%                         [ky,kx,~] = ind2sub([y_pix,x_pix],ROI_Idx{kROI});
%                         kx = round(mean(kx));
%                         ky = round(mean(ky));
%
%                         % check if they are too close (actually same region)
%                         if jx<=(kx+xclose) && jx>=(kx-xclose) && jy<=(ky+yclose) && jy>=(ky-yclose) % only %3 of pixels apart
%                             spatialcorr(kROI)  = 1;
%                         end
%                     end
%
%                     if exist('spatialcorr','var')
%                         indx = find(spatialcorr>0);
%                         temp2.overlap{iPeak,1}= CSArray_Ch1_Hand(1,itrial).calcFindROIs.data.roiNames{indx};
%                     else
%                         temp2.overlap{iPeak,1} = 0;
%                     end
%                     clear spatialcorr
%                 end
%
%                 isFirst = (itrial == 1 );
%                 if isFirst
%                     data.Trial = {};
%                     data.Animal = {};
%                     data.Channel = {};
%                     data.Spot = {};
%                     data.Condition = {};
%                     data.Depth = {};
%                     data.area = {};
%                     data.overlap = {};
%                     data.pixelsize={};
%                 end
%                 data.Trial= [data.Trial; temp2.trialname];
%                 data.Animal= [data.Animal; temp2.animalname];
%                 data.Channel= [data.Channel; temp2.channel];
%                 data.Spot= [data.Spot; temp2.Spot];
%                 data.Condition= [data.Condition; temp2.Cond];
%                 data.Depth= [data.Depth; temp2.depth];
%                 data.area= [data.area; temp2.area];
%                 data.overlap= [data.overlap; temp2.overlap];
%                 data.pixelsize= [data.pixelsize; temp2.pixelsize];
%
%                 %traces output processes
%                 if strcmp(CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.roiNames{1,1}, 'none')
%                     continue
%                 else
%                     traces= CSArray_Ch1_FLIKA(1,itrial).calcMeasureROIs.data.tracesNorm;
%                     %preallocate
%                     Trace_data=cell(size(traces,2),10);
%                     for iROI = 1:size(traces,2)
%                         Trace_data{iROI,1}= CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.roiNames{iROI,1};
%                         Trace_data{iROI,2}= strcat('trial', num2str(itrial,'%02d'));
%                         Trace_data{iROI,3}= 'GCaMP';
%                         Trace_data{iROI,4}= spotId;
%                         Trace_data{iROI,5}= CurrentAnimal;
%                         Trace_data{iROI,6}= CurrentCondition;
%                         Trace_data{iROI,7} = CurrentDepth(1,1);
%                         Trace_data{iROI,8} = traces(:,iROI);
%                         Trace_data{iROI,9} = CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.centroid{iROI,1};
%                         Trace_data{iROI,10} = CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.puffIdxs{iROI,1};
%                         Trace_data{iROI,11} = CSArray_Ch1_FLIKA(1,itrial).rawImg.metadata.pixelSize;
%
%                         % get the indices  and area for a particular ROI
%                         ROIpuffIdx = CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.puffIdxs{iROI,1};
%                         x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
%                         [jy,jx,~] = ind2sub([y_pix,x_pix],ROIpuffIdx);
%                         jx = round(mean(jx));
%                         jy = round(mean(jy));
%                         xclose = round(x_pix*0.03); % 3 percent of pixels
%                         yclose = round(y_pix*0.03); % 3 percent of pixels
%
%                         % find astrocyte process and soma ROIs that have
%                         % similar centroids
%                         for kROI= 1:length(CSArray_Ch1_Hand(1,itrial).calcFindROIs.data.roiNames)
%                             % get the indices of the handclicked ROIs
%                             ROIMask{kROI}= CSArray_Ch1_Hand(1, 1).calcFindROIs.data.roiMask(:,:,kROI);
%                             ROI_Idx{kROI} = find(ROIMask{kROI});
%                             % mean pixels for soma ROI
%                             [ky,kx,~] = ind2sub([y_pix,x_pix],ROI_Idx{kROI});
%                             kx = round(mean(kx));
%                             ky = round(mean(ky));
%
%                             % check if they are too close (actually same region)
%                             if jx<=(kx+xclose) && jx>=(kx-xclose) && jy<=(ky+yclose) && jy>=(ky-yclose) % only %3 of pixels apart
%                                 spatialcorr(kROI)  = 1;
%                             end
%                         end
%
%                         if exist('spatialcorr','var')
%                             indx = find(spatialcorr>0);
%                             Trace_data{iROI,12}= CSArray_Ch1_Hand(1,itrial).calcFindROIs.data.roiNames{indx};
%                         else
%                             Trace_data{iROI,12} = 0;
%                         end
%                         clear spatialcorr
%
%
%                     end
%                     All_traces=vertcat(All_traces, Trace_data);
%                     clearvars Trace_data
%                 end
%             end
%             clearvars temp temp2
%
%             % Astrocyte Hand clicked ROIs
%             for itrial=1:length(CSArray_Ch1_Hand)
%
%                 temp=CSArray_Ch1_Hand(1,itrial).calcDetectSigs.data;
%                 temp2.trialname ={};
%                 temp2.animalname = {};
%                 temp2.channel = {};
%                 temp2.Spot = {};
%                 temp2.Cond = {};
%                 temp2.depth = {};
%                 temp2.overlap = {};
%                 temp2.area = {};
%                 temp2.pixelsize = {};
%
%                 % extract fields from Class
%                 for jField = 1:numel(listFields)
%                     data.(listFields{jField}) = [data.(listFields{jField}); ...
%                         temp.(listFields{jField})];
%                 end
%
%                 % create fields for trial, animal, spot, condition, etc.
%                 for iPeak = 1:length(temp.amplitude)
%                     temp2.trialname{iPeak,1}=strcat('trial', num2str(itrial,'%02d'));
%                     temp2.channel{iPeak,1}= 'GCaMP';
%                     temp2.Spot{iPeak,1}= spotId;
%                     temp2.animalname{iPeak,1}= CurrentAnimal;
%                     temp2.Cond{iPeak,1} = CurrentCondition;
%                     temp2.depth{iPeak,1} = CurrentDepth(1,1);
%                     temp2.area{iPeak,1} = 0;
%                     temp2.overlap{iPeak,1} = 0;
%                     temp2.pixelsize{iPeak,1} = CSArray_Ch1_FLIKA(1,itrial).rawImg.metadata.pixelSize;
%
%                 end
%
%                 data.Trial= [data.Trial; temp2.trialname];
%                 data.Animal= [data.Animal; temp2.animalname];
%                 data.Channel= [data.Channel; temp2.channel];
%                 data.Spot= [data.Spot; temp2.Spot];
%                 data.Condition= [data.Condition; temp2.Cond];
%                 data.Depth= [data.Depth; temp2.depth];
%                 data.area= [data.area; temp2.area];
%                 data.overlap= [data.overlap; temp2.overlap];
%                 data.pixelsize= [data.pixelsize; temp2.pixelsize];
%
%                 %traces output
%                 traces= CSArray_Ch1_Hand(1,itrial).calcMeasureROIs.data.tracesNorm;
%                 %preallocate
%                 Trace_data=cell(size(traces,2),10);
%                 for iROI = 1:size(traces,2)
%                     Trace_data{iROI,1}= CSArray_Ch1_Hand(1,itrial).calcFindROIs.data.roiNames{iROI,1};
%                     Trace_data{iROI,2}= strcat('trial', num2str(itrial,'%02d'));
%                     Trace_data{iROI,3}= 'GCaMP';
%                     Trace_data{iROI,4}= spotId;
%                     Trace_data{iROI,5}= CurrentAnimal;
%                     Trace_data{iROI,6}= CurrentCondition;
%                     Trace_data{iROI,7} = CurrentDepth(1,1);
%                     Trace_data{iROI,8} = traces(:,iROI);
%                     Trace_data{iROI,9} = 0;
%                     Trace_data{iROI,10} = CSArray_Ch1_Hand(1,itrial).calcFindROIs.data.roiMask(:,:,iROI);
%                     Trace_data{iROI,12} = 0;
%                     Trace_data{iROI,11} = CSArray_Ch1_Hand(1,itrial).rawImg.metadata.pixelSize;
%
%                 end
%                 All_traces=vertcat(All_traces, Trace_data);
%                 clearvars Trace_data
%             end
%             clearvars temp temp2
%
%
%             % Neuronal 2D FLIKA
%             for itrial=1:length(CSArray_Ch2_FLIKA)
%
%                 % peak output
%                 temp=CSArray_Ch2_FLIKA(1,itrial).calcDetectSigs.data;
%                 temp2.trialname ={};
%                 temp2.animalname = {};
%                 temp2.channel = {};
%                 temp2.Spot = {};
%                 temp2.Cond = {};
%                 temp2.depth = {};
%                 temp2.overlap = {};
%                 temp2.area = {};
%                 temp2.pixelsize = {};
%
%                 % extract fields from Class
%                 for jField = 1:numel(listFields)
%                     data.(listFields{jField}) = [data.(listFields{jField}); ...
%                         temp.(listFields{jField})];
%                 end
%
%                 % create fields for trial, animal, spot, condition, etc.
%                 for iPeak = 1:length(temp.amplitude)
%                     temp2.trialname{iPeak,1}=strcat('trial', num2str(itrial,'%02d'));
%                     temp2.channel{iPeak,1}= 'RCaMP';
%                     temp2.Spot{iPeak,1}= spotId;
%                     temp2.animalname{iPeak,1}= CurrentAnimal;
%                     temp2.Cond{iPeak,1} = CurrentCondition;
%                     temp2.depth{iPeak,1} = CurrentDepth(1,1);
%                     temp2.pixelsize{iPeak,1} = CSArray_Ch2_FLIKA(1,itrial).rawImg.metadata.pixelSize;
%
%                     % get the indices  and area for a particular ROI
%                     jROIname = temp.ROIname{iPeak};
%                     ROIindex= strcmp(CSArray_Ch2_FLIKA(1,itrial).calcFindROIs.data.roiNames,jROIname);
%                     temp2.area{iPeak,1} = CSArray_Ch2_FLIKA(1,itrial).calcFindROIs.data.area(ROIindex);
%
%                     ROIpuffIdx = CSArray_Ch2_FLIKA(1,itrial).calcFindROIs.data.puffIdxs(ROIindex);
%                     x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
%                     [jy,jx,~] = ind2sub([y_pix,x_pix],cell2mat(ROIpuffIdx));
%                     jx = round(mean(jx));
%                     jy = round(mean(jy));
%                     xclose = round(x_pix*0.03); % 3 percent of pixels
%                     yclose = round(y_pix*0.03); % 3 percent of pixels
%
%                     % find astrocyte process and soma ROIs that have
%                     % similar centroids
%                     for kROI= 1:length(CSArray_Ch2_Hand(1,itrial).calcFindROIs.data.roiNames)
%                         % get the indices of the handclicked ROIs
%                         ROIMask{kROI}= CSArray_Ch2_Hand(1, 1).calcFindROIs.data.roiMask(:,:,kROI);
%                         ROI_Idx{kROI} = find(ROIMask{kROI});
%                         % mean pixels for soma ROI
%                         [ky,kx,~] = ind2sub([y_pix,x_pix],ROI_Idx{kROI});
%                         kx = round(mean(kx));
%                         ky = round(mean(ky));
%
%                         % check if they are too close (actually same region)
%                         if jx<=(kx+xclose) && jx>=(kx-xclose) && jy<=(ky+yclose) && jy>=(ky-yclose) % only %3 of pixels apart
%                             spatialcorr(kROI)  = 1;
%                         end
%                     end
%
%                     if exist('spatialcorr','var')
%                         indx = find(spatialcorr>0);
%                         temp2.overlap{iPeak,1}= CSArray_Ch2_Hand(1,itrial).calcFindROIs.data.roiNames{indx};
%                     else
%                         temp2.overlap{iPeak,1} = 0;
%                     end
%                     clear spatialcorr
%                 end
%
%                 data.Trial= [data.Trial; temp2.trialname];
%                 data.Animal= [data.Animal; temp2.animalname];
%                 data.Channel= [data.Channel; temp2.channel];
%                 data.Spot= [data.Spot; temp2.Spot];
%                 data.Condition= [data.Condition; temp2.Cond];
%                 data.Depth= [data.Depth; temp2.depth];
%                 data.area= [data.area; temp2.area];
%                 data.overlap= [data.overlap; temp2.overlap];
%                 data.pixelsize= [data.pixelsize; temp2.pixelsize];
%
%                 %traces output processes
%                 if strcmp(CSArray_Ch2_FLIKA(1,itrial).calcFindROIs.data.roiNames{1,1}, 'none')
%                     continue
%                 else
%                     traces= CSArray_Ch2_FLIKA(1,itrial).calcMeasureROIs.data.tracesNorm;
%                     %preallocate
%                     Trace_data=cell(size(traces,2),10);
%                     for iROI = 1:size(traces,2)
%                         Trace_data{iROI,1}= CSArray_Ch2_FLIKA(1,itrial).calcFindROIs.data.roiNames{iROI,1};
%                         Trace_data{iROI,2}= strcat('trial', num2str(itrial,'%02d'));
%                         Trace_data{iROI,3}= 'RCaMP';
%                         Trace_data{iROI,4}= spotId;
%                         Trace_data{iROI,5}= CurrentAnimal;
%                         Trace_data{iROI,6}= CurrentCondition;
%                         Trace_data{iROI,7} = CurrentDepth(1,1);
%                         Trace_data{iROI,8} = traces(:,iROI);
%                         Trace_data{iROI,9} = CSArray_Ch2_FLIKA(1,itrial).calcFindROIs.data.centroid{iROI,1};
%                         Trace_data{iROI,10} = CSArray_Ch2_FLIKA(1,itrial).calcFindROIs.data.puffIdxs{iROI,1};
%                         Trace_data{iROI,11} = CSArray_Ch2_FLIKA(1,itrial).rawImg.metadata.pixelSize;
%
%                         % get the indices  and area for a particular ROI
%                         ROIpuffIdx = CSArray_Ch2_FLIKA(1,itrial).calcFindROIs.data.puffIdxs{iROI,1};
%                         x_pix= Settings.Xres(1,1); y_pix= Settings.Yres(1,1);
%                         [jy,jx,~] = ind2sub([y_pix,x_pix],ROIpuffIdx);
%                         jx = round(mean(jx));
%                         jy = round(mean(jy));
%                         xclose = round(x_pix*0.03); % 3 percent of pixels
%                         yclose = round(y_pix*0.03); % 3 percent of pixels
%
%                         % find astrocyte process and soma ROIs that have
%                         % similar centroids
%                         for kROI= 1:length(CSArray_Ch2_Hand(1,itrial).calcFindROIs.data.roiNames)
%                             % get the indices of the handclicked ROIs
%                             ROIMask{kROI}= CSArray_Ch2_Hand(1, 1).calcFindROIs.data.roiMask(:,:,kROI);
%                             ROI_Idx{kROI} = find(ROIMask{kROI});
%                             % mean pixels for soma ROI
%                             [ky,kx,~] = ind2sub([y_pix,x_pix],ROI_Idx{kROI});
%                             kx = round(mean(kx));
%                             ky = round(mean(ky));
%
%                             % check if they are too close (actually same region)
%                             if jx<=(kx+xclose) && jx>=(kx-xclose) && jy<=(ky+yclose) && jy>=(ky-yclose) % only %3 of pixels apart
%                                 spatialcorr(kROI)  = 1;
%                             end
%                         end
%
%                         if exist('spatialcorr','var')
%                             indx = find(spatialcorr>0);
%                             Trace_data{iROI,12}= CSArray_Ch2_Hand(1,itrial).calcFindROIs.data.roiNames{indx};
%                         else
%                             Trace_data{iROI,12} = 0;
%                         end
%                         clear spatialcorr
%
%
%                     end
%                     All_traces=vertcat(All_traces, Trace_data);
%                     clearvars Trace_data
%                 end
%             end
%             clearvars temp temp2
%
%
%             % Neuronal Hand clicked ROIs
%             for itrial=1:length(CSArray_Ch2_Hand)
%
%                 temp=CSArray_Ch2_Hand(1,itrial).calcDetectSigs.data;
%                 temp2.trialname ={};
%                 temp2.animalname = {};
%                 temp2.channel = {};
%                 temp2.Spot = {};
%                 temp2.Cond = {};
%                 temp2.depth = {};
%                 temp2.overlap = {};
%                 temp2.area = {};
%                 temp2.pixelsize = {};
%
%                 % extract fields from Class
%                 for jField = 1:numel(listFields)
%                     data.(listFields{jField}) = [data.(listFields{jField}); ...
%                         temp.(listFields{jField})];
%                 end
%
%                 % create fields for trial, animal, spot, condition, etc.
%                 for iPeak = 1:length(temp.amplitude)
%                     temp2.trialname{iPeak,1}=strcat('trial', num2str(itrial,'%02d'));
%                     temp2.channel{iPeak,1}= 'RCaMP';
%                     temp2.Spot{iPeak,1}= spotId;
%                     temp2.animalname{iPeak,1}= CurrentAnimal;
%                     temp2.Cond{iPeak,1} = CurrentCondition;
%                     temp2.depth{iPeak,1} = CurrentDepth(1,1);
%                     temp2.area{iPeak,1} = 0;
%                     temp2.overlap{iPeak,1} = 0;
%                     temp2.pixelsize{iPeak,1} = CSArray_Ch1_FLIKA(1,itrial).rawImg.metadata.pixelSize;
%
%                 end
%
%                 data.Trial= [data.Trial; temp2.trialname];
%                 data.Animal= [data.Animal; temp2.animalname];
%                 data.Channel= [data.Channel; temp2.channel];
%                 data.Spot= [data.Spot; temp2.Spot];
%                 data.Condition= [data.Condition; temp2.Cond];
%                 data.Depth= [data.Depth; temp2.depth];
%                 data.area= [data.area; temp2.area];
%                 data.overlap= [data.overlap; temp2.overlap];
%                 data.pixelsize= [data.pixelsize; temp2.pixelsize];
%
%                 %traces output
%                 traces= CSArray_Ch2_Hand(1,itrial).calcMeasureROIs.data.tracesNorm;
%                 %preallocate
%                 Trace_data=cell(size(traces,2),10);
%                 for iROI = 1:size(traces,2)
%                     Trace_data{iROI,1}= CSArray_Ch2_Hand(1,itrial).calcFindROIs.data.roiNames{iROI,1};
%                     Trace_data{iROI,2}= strcat('trial', num2str(itrial,'%02d'));
%                     Trace_data{iROI,3}= 'RCaMP';
%                     Trace_data{iROI,4}= spotId;
%                     Trace_data{iROI,5}= CurrentAnimal;
%                     Trace_data{iROI,6}= CurrentCondition;
%                     Trace_data{iROI,7} = CurrentDepth(1,1);
%                     Trace_data{iROI,8} = traces(:,iROI);
%                     Trace_data{iROI,9} = 0;
%                     Trace_data{iROI,10} = CSArray_Ch2_Hand(1,itrial).calcFindROIs.data.roiMask(:,:,iROI);
%                     Trace_data{iROI,12} = 0;
%                     Trace_data{iROI,11} = CSArray_Ch2_Hand(1,itrial).rawImg.metadata.pixelSize;
%
%                 end
%                 All_traces=vertcat(All_traces, Trace_data);
%                 clearvars Trace_data
%
%             end
%             dataNames=fieldnames(data);
%             data2= struct2cell(data);
%             data3= [data2{:}];
%
%             AllData=vertcat(AllData, data3);
%
%
%             clearvars data data3 temp temp2
%         end
%     end
% end
%
% % %% Save all data for R analysis
% AllData2= [dataNames';AllData];
%
% cd(fullfile(Settings.MainDir, 'Results'));
% % write date to created file
% cell2csv(SaveFiles{1,1}, AllData2);
% save(SaveFiles{1,3}, 'All_traces','-v7.3');
% save(SaveFiles{1,2}, 'AllData2','-v7.3');
