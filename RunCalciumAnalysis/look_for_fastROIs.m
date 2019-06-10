%% 1 lck no stim
clearvars
%addpath(genpath('C:\Program Files\2p-img-analysis'));
%addpath(genpath('C:\Program Files\CHIPS'));
%addpath(genpath('E:\Jill\matlab'));
%AllData= [];
All_traces= [];


%% Information about your images
channel = struct('Ca_Memb_Astro',2,'Ca_Neuron',1);

calibration ='E:\matlab\CalibrationFiles\calibration_20x.mat';
CalibrationFile = CalibrationPixelSize.load(calibration);

FastTime=1;
%% make a mixing matrix for spectral unmixing of RCaMP and GCaMP

load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\RCaMP_mGCaMP_Matrix_4Ch.mat');



%% load images
      
            
            CurrentCondition = 'Stim';
            BL_frames = 59*2;
            
            
            % Create an array of ScanImage Tiffs
                        ImgArray =  BioFormats([], channel, CalibrationFile);
           
            
            % Spectral Unmixing of GCaMP and RCaMP
            ImgArray= ImgArray.unmix_chs(false, [], cell2mat(RCaMP_mGCaMP_Matrix));
            
            % Run motion correction          
            
            % Extract a reference image
            refImg = mean(HighRes.rawdata(:,:,2,:),4);
            %figure, imagesc(refImg), axis image, axis off, colormap(gray)
            
            ImgArray = ImgArray.motion_correct('refImg',refImg, 'ch',2,'maxShift', 10,'minCorr', 0.4);%,'doPlot',true);
            % fill in bad data?
         
               % ImgArray(1,3).plot();

            
            %% Configs for Finding ROIs
            %ASTROCYTES
            % 2D automated selection for peaks
            AC_findConf{1} = ConfigFindROIsFLIKA_2D.from_preset('ca_memb_astro', 'baselineFrames',...
                BL_frames,'freqPassBand',1,'sigmaXY', 2,...
                'sigmaT', 0.1,'thresholdPuff', 7, 'threshold2D', 0.2,...
                'minRiseTime',0.0845, 'maxRiseTime', 1,'minPuffArea', 10,...
                'dilateXY', 5, 'dilateT', 0.3,'erodeXY', 1, 'erodeT', 0.1,...
                'discardBorderROIs',false);
            
            % 3D automated selection for time and space estimations
            AC_findConf{3} = ConfigFindROIsFLIKA_3D.from_preset('ca_memb_astro', 'baselineFrames',...
                BL_frames,'freqPassBand',1,'sigmaXY', 2,...
                'sigmaT', 0.1,'thresholdPuff', 7, 'threshold2D', 0.2,...
                'minRiseTime',0.0845, 'maxRiseTime', 1,'minPuffArea', 10,...
                'dilateXY', 5, 'dilateT', 0.3,'erodeXY', 1, 'erodeT', 0.1,...
                'discardBorderROIs',false);
            
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
            configCS{1,3} = ConfigCellScan(AC_findConf{1,3}, measureConf,detectConf{1,3}); % astrocyte 3D FLIKA
            
            %% Create CellScan objects
            CSArray_Ch1_FLIKA = CellScan([], ImgArray, configCS{1,1}, 1);
                  
            
            
            %% Process the images
            CSArray_Ch1_FLIKA =CSArray_Ch1_FLIKA.process();

           % CSArray_Ch1_FLIKA.plot();
            %CSArray_Ch1_FLIKA.opt_config();
            %% Output data
            
            % Astrocyte FLIKA
            for itrial=1:length(CSArray_Ch1_FLIKA)

                %traces output processes
                if strcmp(CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.roiNames{1,1}, 'none')
                    continue
                else
                    traces= CSArray_Ch1_FLIKA(1,itrial).calcMeasureROIs.data.tracesNorm;
                    %preallocate
                    Trace_data=cell(size(traces,2),10);
                    for iROI = 1:size(traces,2)
                        Trace_data{iROI,1}= CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.roiNames{iROI,1};
                        Trace_data{iROI,2}= strcat('trial', num2str(itrial,'%02d'));
                        Trace_data{iROI,3}= 'GCaMP';
                        Trace_data{iROI,4}= NaN;
                        Trace_data{iROI,5}= NaN;
                        Trace_data{iROI,6}= CurrentCondition;
                        Trace_data{iROI,7} = NaN;
                        Trace_data{iROI,8} = traces(:,iROI);
                        Trace_data{iROI,9} = CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.centroidX(iROI,1);
                        Trace_data{iROI,10} = CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.centroidY(iROI,1);
                        Trace_data{iROI,11} = CSArray_Ch1_FLIKA(1,itrial).calcFindROIs.data.roiIdxs{iROI,1};
                        Trace_data{iROI,12} = CSArray_Ch1_FLIKA(1,itrial).rawImg.metadata.pixelSize;
                        
                    end
                    All_traces=vertcat(All_traces, Trace_data);
                    clearvars Trace_data
                end
            end

            
            
 %% find onset times
 
FrameRate=CSArray_Ch1_FLIKA(1,1).rawImg.metadata.frameRate;
nframes=CSArray_Ch1_FLIKA(1,1).rawImg.metadata.nFrames;
TimeX(1:nframes) = (1:nframes)/FrameRate;

baselineCorrectedTime=TimeX-(BL_frames/FrameRate);

% peak onsets and AUC in the first second after stim for each ROI
for iROI= 1:size(All_traces,1)
    trace=All_traces{iROI,8};
    % onset time
    if size(trace,1)>590
        Onsets=find_first_onset_time(baselineCorrectedTime(10:end), trace(10:592),2.5,2);
        if isempty(Onsets)
            Onsets=nan(1,1);
        end
        All_traces{iROI, 16}= Onsets;

    else
        All_traces{iROI,16}=NaN;
    end    
end


% plot onset times histogram
% X=cell2mat(All_traces(:,16));
% X(isnan(X)) = [];
% figure
% histc(X, 'BinWidth',0.5);

%% show fast ROIs
for iROI=1:size(All_traces, 1)
  
fastIdx(iROI)=~isempty(find(All_traces{iROI,16}<FastTime));
delayedIdx(iROI)=~isempty(find(All_traces{iROI,16}<12 && All_traces{iROI,16}>FastTime));
end

fastROIs=All_traces(fastIdx,:);
slowROIs=All_traces(delayedIdx,:);

ProcMap1=zeros(127,128);
for iFast=1:size(fastROIs,1)
    Image1=zeros(127,128);
    Image1(fastROIs{iFast,11})=1;
    ProcMap1=ProcMap1+Image1;
end
ProcMap2=zeros(1,128);
ProcMaps=[ProcMap1;ProcMap2];


ACMask=im2bw(ProcMaps);
AC_B=bwboundaries(ACMask);

figure('name','fastROI mask');
%imshow(mean(CSArray_Ch1_FLIKA(1,1).rawImg.rawdata(:,:,2,:),4));%((128,128)); 
imshow(zeros(128,128)); 
hold on
for k=1:length(AC_B)
    border=AC_B{k};
    plot(border(:,2),border(:,1),'g','linewidth',1.5);
end

clearvars fastIdx TimeX


%% Plot fastROIs
nframes=round(30*FrameRate);

figure ('name', 'fast ROI traces')
hold on
axis off
ROIs=39:45;
for ii=1:size(fastROIs,1)
    
    tempY1=smooth(fastROIs{ii,8},3);
    tempY1=tempY1(1:nframes);
    plot(baselineCorrectedTime(1:nframes),tempY1'+(3*(ii-1)),'g')%'LineWidth',1);
    plot([0 0],[0 50], 'k','LineWidth', 1)
end

figure ('name', 'fast ROI traces- 4')
hold on
axis off
ROIs=[10,4];
for ii=1:size(ROIs,2)
    ROInum=ROIs(ii);
    tempY1=smooth(fastROIs{ROInum,8},5);
    tempY1=tempY1(1:nframes);    
    plot(baselineCorrectedTime(1:nframes),tempY1'+(3*(ii-1)),'b')%'LineWidth',1);
end
plot([0 0],[0 6], 'k','LineWidth', 1)
plot([-10 -10],[0 1], 'k','LineWidth', 1)




figure ('name', 'slow ROIs trial 4')
hold on
axis off
ROIs=39:45;
for ii=1:size(ROIs,2)
    ROInum=ROIs(ii);
    tempY1=smooth(slowROIs{ROInum,8},3);
    tempY1=tempY1(1:nframes);    
    plot(baselineCorrectedTime(1:nframes),tempY1'+(3*(ii-1)),'g')%'LineWidth',1);
end
plot([0 0],[0 40], 'k','LineWidth', 1)

figure ('name', 'slow ROIs trial 8')
hold on
axis off
ROIs=94:107;
for ii=1:size(ROIs,2)
    ROInum=ROIs(ii);
    tempY1=smooth(slowROIs{ROInum,8},3);
    tempY1=tempY1(1:nframes);    
    plot(baselineCorrectedTime(1:nframes),tempY1'+(3*(ii-1)),'g')%'LineWidth',1);
end
plot([0 0],[0 40], 'k','LineWidth', 1)

% PLOT THE ROI MASK

figure ('name', 'slow ROIs trial 1')
hold on
axis off
ROIs=1:7;
for ii=1:size(ROIs,2)
    ROInum=ROIs(ii);
    tempY1=smooth(slowROIs{ROInum,8},3);
    tempY1=tempY1(1:nframes);    
    plot(baselineCorrectedTime(1:nframes),tempY1'+(3*(ii-1)),'g')%'LineWidth',1);
end
plot([0 0],[0 40], 'k','LineWidth', 1)


figure ('name', 'fast and slow ROI traces- trial 4')
hold on
axis off
ROIs=[10,4];
for ii=1:size(ROIs,2)
    ROInum=ROIs(ii);
    tempY1=smooth(fastROIs{ROInum,8},3);
    tempY1=tempY1(1:nframes);    
    plot(baselineCorrectedTime(1:nframes),tempY1'+(3*(ii-1)),'b')%'LineWidth',1);
end
ROIs=[41,42];
for ii=1:size(ROIs,2)
    ROInum=ROIs(ii);
    tempY1=smooth(slowROIs{ROInum,8},3);
    tempY1=tempY1(1:nframes);    
    plot(baselineCorrectedTime(1:nframes),tempY1'+(5*ii+2),'g')%'LineWidth',1);
end
rectangle('Position', [0 -1 8 20])
plot([-10 -10],[0 1], 'k','LineWidth', 1)




%% extract ROI masks


% fast ROIs
Mask1=zeros(127,128);
Mask1(fastROIs{10,11})=1;
Mask1(fastROIs{4,11})=1;
Mask1(slowROIs{41,11})=1;
Mask1(slowROIs{42,11})=1;


ProcMap2=zeros(1,128);
ProcMaps=[Mask1;ProcMap2];

ACMask=im2bw(ProcMaps);
AC_B=bwboundaries(ACMask);

figure('name','example ROI mask');
imshow(zeros(128,128)); 
hold on
for k=1:length(AC_B)
    border=AC_B{k};
    plot(border(:,2),border(:,1),'g','linewidth',1.5);
end

%% extract ROI masks


% fast ROIs
Mask1=zeros(127,128);
Mask1(fastROIs{10,11})=1;
Mask1(fastROIs{4,11})=1;
Mask1(slowROIs{41,11})=1;
Mask1(slowROIs{42,11})=1;


ProcMap2=zeros(1,128);
ProcMaps=[Mask1;ProcMap2];

ACMask=im2bw(ProcMaps);
AC_B=bwboundaries(ACMask);

figure('name','example ROI mask');
imshow(zeros(128,128)); 
hold on
for k=1:length(AC_B)
    border=AC_B{k};
    plot(border(:,2),border(:,1),'g','linewidth',1.5);
end
%% cytosolic

nframes=round(30*FrameRate);

figure ('name', 'fast ROI traces')
hold on
axis off
ROIs=99:118;
for ii=1:size(ROIs,2)    
    ROInum=ROIs(ii);
    tempY1=smooth(All_traces{ii,10},3);
    tempY1=tempY1(1:nframes);
    plot(baselineCorrectedTime(1:nframes),tempY1'+(3*(ii-1)),'g')%'LineWidth',1);
   % plot([0 0],[0 50], 'k','LineWidth', 1)
end
rectangle('Position', [0 -1 8 50])

figure ('name', 'slow ROIs trial 4')
hold on
axis off
ROIs=[106, 115];
for ii=1:size(ROIs,2)
    ROInum=ROIs(ii);
    tempY1=smooth(All_traces{ROInum,10},5);
    tempY1=tempY1(1:nframes);    
    plot(baselineCorrectedTime(1:nframes),tempY1'+(3*(ii-1)),'g')%'LineWidth',1);
end
plot([0 0],[0 50], 'k','LineWidth', 1)

figure ('name', 'neuron traces')
hold on
axis off
ROIs=[119:155];
for ii=1:size(ROIs,2)    
    ROInum=ROIs(ii);
    tempY1=smooth(All_traces{ii,10},3);
    tempY1=tempY1(1:nframes);
    plot(baselineCorrectedTime(1:nframes),tempY1'+(3*(ii-1)),'r')%'LineWidth',1);
   % plot([0 0],[0 50], 'k','LineWidth', 1)
end
rectangle('Position', [0 -1 8 50])


