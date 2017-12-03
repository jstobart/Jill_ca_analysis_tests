%% 1 lck no stim
clearvars

AllData= [];
All_traces= [];


%% Information about your images
channel = struct('Ca_Memb_Astro',1,'Ca_Neuron',2);


%How fast should we consider

FastTime=2;
%% make a mixing matrix for spectral unmixing of RCaMP and GCaMP

load(fullfile(Settings.MainDir, 'Results','RCaMP_mGCaMP_Matrix.mat'));



%% load images
      
            
            CurrentCondition = 'Stim';
            BL_frames = 59;
            
            
            % Create an array of ScanImage Tiffs
            ImgArray =  SCIM_Tif([], channel, []);
            
            % Spectral Unmixing of GCaMP and RCaMP
            ImgArray= ImgArray.unmix_chs(false, [], cell2mat(RCaMP_mGCaMP_Matrix));
            
            % Run motion correction          
            HighRes = SCIM_Tif([],channel, []);
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
                'sigmaT', 0.1,'threshold_std', 7, 'threshold_2D', 0.2,...
                'min_rise_time',0.0845, 'max_rise_time', 1,'minPuffArea', 10,...
                'dilateXY', 5, 'dilateT', 0.3,'erodeXY', 1, 'erodeT', 0.1,...
                'discardBorderROIs',true);
            
            % 3D automated selection for time and space estimations
            AC_findConf{3} = ConfigFindROIsFLIKA_3D.from_preset('ca_memb_astro', 'baselineFrames',...
                BL_frames,'freqPassBand',1,'sigmaXY', 2,...
                'sigmaT', 0.1,'threshold_std', 7, 'threshold_2D', 0.2,...
                'min_rise_time',0.0845, 'max_rise_time', 1,'minPuffArea', 10,...
                'dilateXY', 5, 'dilateT', 0.3,'erodeXY', 1, 'erodeT', 0.1,...
                'discardBorderROIs',true);
            
            %% Configuration for measuring ROIs
            % AWAKE astrocyte membrane calcium
            detectConf{1} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,...'normMethod', 'z-score','zIters', 100,...
                'propagateNaNs', false, 'excludeNaNs', false, 'lpWindowTime', 1.5, 'spFilterOrder', 2,...
                'spPassBandMin',0.05, 'spPassBandMax', 0.5, 'thresholdSD_low', 3,'thresholdSD_band', 5);
            
            % AWAKE neuron calcium
            detectConf{2} = ConfigDetectSigsClsfy('baselineFrames', BL_frames,... 'normMethod','z-score',... 'zIters', 10000,...
                'propagateNaNs', false,'excludeNaNs', false, 'lpWindowTime', 2, 'spFilterOrder', 2,...
                'spPassBandMin',0.1, 'spPassBandMax', 1, 'thresholdSD_low', 3,'thresholdSD_band', 5);
            
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

            CSArray_Ch1_FLIKA.plot();
            
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

baselineCorrectedTime=TimeX-5;

% peak onsets and AUC in the first second after stim for each ROI
for iROI= 1:length(All_traces)
    trace=All_traces{iROI,8};
    %first 1 sec after stim onset
    x1=round(FrameRate*5);
    x2=round(FrameRate*6);
    x3= round(FrameRate*10);
    x4= round(FrameRate*15);
    % onset time
    if size(trace,1)>590
        Onsets=find_first_onset_time(baselineCorrectedTime(10:end), trace(10:592),2.5,2);
        if isempty(Onsets)
            Onsets=nan(1,1);
        end
        All_traces{iROI, 16}= Onsets;
        % AUC
        All_traces{iROI,17}=trapz(trace(x1:x2));
        All_traces{iROI,18}=trapz(trace(x3:x4));
    else
        All_traces{iROI,16}=NaN;
        All_traces{iROI,17}=NaN;
        All_traces{iROI,18}=NaN;
    end    
end


% plot onset times histogram

figure
histogram(cell2mat(All_traces(:,16)), 'BinWidth',0.5);

%% show fast ROIs
for iROI=1:size(All_traces, 1)
  
fastIdx(iROI)=~isempty(find(All_traces{iROI,16}<FastTime));

end

fastROIs=All_traces(fastIdx,:);

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

figure();
imshow(zeros(128,128)); hold on
for k=1:length(AC_B)
    border=AC_B{k};
    plot(border(:,2),border(:,1),'g','linewidth',1.5);
end





