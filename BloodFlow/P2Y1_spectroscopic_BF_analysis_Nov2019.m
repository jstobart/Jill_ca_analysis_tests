close all; clear variables;

All_traces=[]; AllData=[]; All_traces2=[]; AllData2=[]; LckData=[]; LckFieldData=[];
%% Information about your images

% folder where data should be saved for each animal
Settings.ResultsFolder = 'D:\Data\GCaMP_RCaMP\NR1_KD\Results';
CalibrationFile= CalibrationPixelSize([1 2 3 4], [10 5 2 1], 240, '4x', '2019', 'Jill', 'Jill');  %((ZOOM, PXSIZE, IMGSIZE, OBJECTIVE, 
        %       DATE, NAME, PERSON, FUNRAW) 

% File names for saving
SaveFiles{1,1} = fullfile(Settings.ResultsFolder,'FilesforMatlab', '74_peaks_longtrials3.mat');
SaveFiles{1,2} = fullfile(Settings.ResultsFolder,'FilesforR','74_peaks_longtrials3.csv');


Settings.FileNames = {  % folder names where images and roiSet.zip is found
    'G:\ZurichData\LaserSpeckle_MWSpec_Data\MWSpec\P2YRG7\baseline_04102017\Trial_1',...
    'G:\ZurichData\LaserSpeckle_MWSpec_Data\MWSpec\P2YRG7\baseline_04102017\Trial_2',...
    };

Settings.Baseline = 5; % time (s) before the whisker stimulator starts,  2s for short trials (10s long), 5s for long trials
Settings.Animal = 'P2YRG7';

channel = struct('intrinsic_optical',[1,2,3,4,5,6]);%,'intrinsic_optical',2,'intrinsic_optical',3,'intrinsic_optical',4,'intrinsic_optical',5,'intrinsic_optical',6);  % what's on each channel


%% Loop through each file and make Cell Scans for LckGCaMP and RCaMP

for Iteration= 1:length(Settings.FileNames)
    
    % loop through each spot
    SpotRoot= Settings.FileNames{Iteration};
    SpotId=SpotRoot(end-24:end);
    
    % Get a list of all dat files in this folder.
    % different stimulation conditions
  
    dataFiles= dir(fullfile(SpotRoot,'*.dat'));
    fnTempList = {dataFiles(:).name};
    dataPath = fullfile(SpotRoot, fnTempList);
            
    % find folders for each condition (stim vs nostim, different
    % frequencies)
    for iTrial= 1:length(dataPath)
        DatName=dataPath(iTrial);
        
           SpectroscopicImageObj = SpectroscopicImage(DatName, channel, CalibrationFile);
            %% Create Spectroscopic Object
            BloodFlow = SpectroscopicOI(DatName,[],[]); % peaks from automated ROIs
            
            % Process the images
            BloodFlow =BloodFlow.process();
            
            % Make the debugging plots
            BloodFlow.plot();

            
    end
end

%             %% Extract/Calculate data we want and create table to output data
% %             if bad
% %                 continue
% %             else
%             % make a giant data table for the signal peaks
%             listFields = {'amplitude', 'fullWidth', 'halfWidth', ...
%                 'numPeaks', 'peakTime', 'peakStart', 'peakStartHalf', ...
%                 'peakType', 'prominence', 'roiName', 'peakAUC'};
%             
%             CellScans=vertcat(BloodFlow, neurons1, neurons2);
%             
%             
%             % loop through cellscans and pull out the data
%             for iScan=1:size(CellScans,1)
%                 
%                 % peak output
%                 temp=CellScans(iScan).calcDetectSigs.data;
%                 temp2.trialname ={};
%                 temp2.animalname = {};
%                 temp2.channel = {};
%                 temp2.Spot = {};
%                 temp2.Cond = {};
%                 temp2.area = {};
%                 temp2.pixelsize = {};
%                 temp2.overlap ={};
%                 
%                 % extract fields from Class
%                 for jField = 1:numel(listFields)
%                     isFirst = (iTrial == 1 && iScan == 1);
%                     if isFirst
%                         data.(listFields{jField}) = {};
%                     end
%                     data.(listFields{jField}) = [data.(listFields{jField}); ...
%                         temp.(listFields{jField})];
%                 end
%                 
%                 % create fields for trial, animal, spot, condition, etc.
%                 for iPeak = 1:length(temp.amplitude)
%                     temp2.trialname{iPeak,1}=strcat('trial', num2str(iTrial,'%02d'));
%                     if iScan==1
%                         temp2.channel{iPeak,1}= 'GCaMP';
%                     else
%                         temp2.channel{iPeak,1}= 'RCaMP';
%                     end
%                     temp2.Spot{iPeak,1}= SpotId;
%                     temp2.animalname{iPeak,1}= Settings.Animal;
%                     temp2.Cond{iPeak,1} = subDirsNames{1,iTrial};
%                     temp2.pixelsize{iPeak,1} = CellScans(iScan).rawImg.metadata.pixelSize;
%                     
%                     % get the indices  and area for a particular ROI
%                     if iScan==3
%                         jROIname = temp.roiName{iPeak};
%                         ROIindex= strcmp(CellScans(iScan,1).calcFindROIs.data.roiNames,jROIname);
%                         temp2.area{iPeak,1} = CellScans(iScan,1).calcFindROIs.data.area(ROIindex);
%                         
%                         jx = CellScans(iScan,1).calcFindROIs.data.centroidX(ROIindex,1);
%                         jy = CellScans(iScan,1).calcFindROIs.data.centroidY(ROIindex,1);
%                         xclose = round(x_pix*0.03); % 3 percent of pixels
%                         yclose = round(y_pix*0.03); % 3 percent of pixels
%                         
%                         if isempty(jx)
%                             temp2.overlap{iPeak,1} = 0;
%                         else
%                             % find NEURONAL FLIKA ROIs and hand selected ROIS that have
%                             % similar centroids
%                             
%                             for kROI= 1:size(CellScans(2,1).calcFindROIs.data.roiNames,1)
%                                 % get the centroids of the handclicked ROIs
%                                 kx = CellScans(2,1).calcFindROIs.data.centroidX(kROI,1);
%                                 ky = CellScans(2,1).calcFindROIs.data.centroidY(kROI,1);
%                                 
%                                 % check if they are too close (actually same region)
%                                 if jx<=(kx+xclose) && jx>=(kx-xclose) && jy<=(ky+yclose) && jy>=(ky-yclose) % only %3 of pixels apart
%                                     spatialcorr(kROI)  = 1;
%                                 end
%                             end
%                             
%                             if exist('spatialcorr','var')
%                                 indx = find(spatialcorr>0);
%                                 temp2.overlap{iPeak,1}= CellScans(2,1).calcFindROIs.data.roiNames{indx};
%                             else
%                                 temp2.overlap{iPeak,1} = 0;
%                             end
%                             
%                             clear spatialcorr
%                         end
%                     else
%                         temp2.area{iPeak,1} = 0;
%                         temp2.overlap{iPeak,1} = 0;
%                     end
%                     
%                 end
%                 
%                 
%                 
%                 %%
%                 isFirst = (iTrial == 1 && iScan == 1);
%                 if isFirst
%                     data.Trial = {};
%                     data.Animal = {};
%                     data.Channel = {};
%                     data.Spot = {};
%                     data.Condition = {};
%                     data.area = {};
%                     data.pixelsize={};
%                     data.overlap ={};
%                 end
%                 data.Trial= [data.Trial; temp2.trialname];
%                 data.Animal= [data.Animal; temp2.animalname];
%                 data.Channel= [data.Channel; temp2.channel];
%                 data.Spot= [data.Spot; temp2.Spot];
%                 data.Condition= [data.Condition; temp2.Cond];
%                 data.area= [data.area; temp2.area];
%                 data.pixelsize= [data.pixelsize; temp2.pixelsize];
%                 data.overlap= [data.overlap; temp2.overlap];
%                 
%                 clearvars temp temp2
%                 
%                 
%                 %% make a table of ROI info and traces
%                 
%                 %traces output processes
%                 if strcmp(CellScans(iScan).calcFindROIs.data.roiNames{1,1}, 'none')
%                     continue
%                 else
%                     traces= CellScans(iScan).calcMeasureROIs.data.tracesNorm;
%                     
%                     %preallocate
%                     Trace_data=cell(size(traces,2),10);
%                     mask = zeros(x_pix, y_pix);
%                     
%                     for iROI = 1:size(traces,2)
%                         Trace_data{iROI,1}= CellScans(iScan).calcFindROIs.data.roiNames{iROI,1};
%                         Trace_data{iROI,2}= strcat('trial', num2str(iTrial,'%02d'));
%                         if iScan==1
%                             Trace_data{iROI,3}= 'GCaMP';
%                         else
%                             Trace_data{iROI,3}= 'RCaMP';
%                         end
%                         
%                         Trace_data{iROI,4}= SpotId;
%                         Trace_data{iROI,5}= Settings.Animal;
%                         Trace_data{iROI,6}= subDirsNames{1,iTrial};
%                         Trace_data{iROI,7} = Settings.Baseline;
%                         Trace_data{iROI,8} = traces(:,iROI);
%                         if iScan==1 || iScan==3
%                             Trace_data{iROI,9} = CellScans(iScan).calcFindROIs.data.roiIdxs{iROI,1};
%                             Trace_data{iROI,10} = CellScans(iScan).rawImg.metadata.pixelSize;
%                         else
%                             
%                             Trace_data{iROI,9} = CellScans(iScan).calcFindROIs.data.roiMask(:,:,iROI);
%                             Trace_data{iROI,10} = CellScans(iScan).rawImg.metadata.pixelSize;
%                         end
%                         
%                         FrameRate= CellScans(1).rawImg.metadata.frameRate;
%                         Trace_data{iROI,11} = FrameRate; % frameRate
%                         
%                         nFrames=length(traces(:,iROI));
%                         TimeX(1:nFrames) = (1:nFrames)/FrameRate;
%                         
%                         % Calculate the first peak onset time and AUC after stim
%                         BL_time=Settings.Baseline;  % number of s for baseline
%                         baselineCorrectedTime=TimeX-BL_time;
%                         
%                         % onset time  % 2.5SD from baseline and
%                         % smoothing trace at 11 points (5 each side
%                         % of middle)
%                         Onsets=find_first_onset_time(baselineCorrectedTime(10:end), traces(10:end,iROI),2.5,2);
%                         if isempty(Onsets)
%                             Onsets=nan(1,1);
%                         end
%                         Trace_data{iROI,12}= Onsets;
%                         
%                         % trace AUC in the 8 s follow stimulation
%                         x2=round(FrameRate*(BL_time+8));
%                         Trace_data{iROI,13}=trapz(traces(BL_frames:x2,iROI));
%                         
%                         % create a ROI mask of all the astrocyte ROIs
%                         if iScan==1
%                             roiIdx = CellScans(1).calcFindROIs.data.roiIdxs{iROI,1};
%                             mask(roiIdx) = true;
%                         end
%                         
%                         % get the indices  and area for a particular ROI
%                         if iScan==3
%                             jx = CellScans(iScan,1).calcFindROIs.data.centroidX(iROI,1);
%                             jy = CellScans(iScan,1).calcFindROIs.data.centroidY(iROI,1);
%                             xclose = round(x_pix*0.03); % 3 percent of pixels
%                             yclose = round(y_pix*0.03); % 3 percent of pixels
%                             
%                             if isempty(jx)
%                                 Trace_data{iROI,14} = 0;
%                             else
%                                 % find FLIKA ROIs and hand selected ROIS that have
%                                 % similar centroids
%                                 for kROI= 1:length(CellScans(2,1).calcFindROIs.data.roiNames)
%                                     % get the centroids of the handclicked ROIs
%                                     kx = CellScans(2,1).calcFindROIs.data.centroidX(kROI,1);
%                                     ky = CellScans(2,1).calcFindROIs.data.centroidY(kROI,1);
%                                     
%                                     % check if they are too close (actually same region)
%                                     if jx<=(kx+xclose) && jx>=(kx-xclose) && jy<=(ky+yclose) && jy>=(ky-yclose) % only %3 of pixels apart
%                                         spatialcorr(kROI)  = 1;
%                                     end
%                                 end
%                                 
%                                 
%                                 if exist('spatialcorr','var')
%                                     indx = find(spatialcorr>0);
%                                     Trace_data{iROI,14}= CellScans(2,1).calcFindROIs.data.roiNames{indx};
%                                 else
%                                     Trace_data{iROI,14} = 0;
%                                 end
%                                 clear spatialcorr
%                             end
%                         else
%                             Trace_data{iROI,14} = 0;
%                         end
%                     end
%                     
%                 end
%                 
%                 All_traces=vertcat(All_traces, Trace_data);
%                 clearvars Trace_data
%                 
%                 %% Lck field of view specific data
%                 if iScan==1
%                     Lck.trialname{iTrial,1} =strcat('trial', num2str(iTrial,'%02d'));
%                     Lck.Spot{iTrial,1}= SpotId;
%                     Lck.animalname{iTrial,1}= Settings.Animal;
%                     Lck.Cond{iTrial,1} = subDirsNames{1,iTrial};
%                     Lck.pixelsize{iTrial,1} = CellScans(iScan).rawImg.metadata.pixelSize;
%                     
%                     % get fraction of active MD area
%                     neuroMask = CellScans(2).calcFindROIs.data.roiMask;
%                     if ndims(neuroMask) == 3
%                         neuroMask = max(neuroMask, [], 3);
%                     end
%                     [nFluoPix, nActivePix, nTotalPix] = ...
%                         getFracActive(CellScans(1), 'nanMask', ...
%                         neuroMask);
%                     Lck.nFluoPix{iTrial,1} = nFluoPix;
%                     Lck.nActivePix{iTrial,1} = nActivePix;
%                     Lck.nTotalPix{iTrial,1} = nTotalPix;
%                     
%                     % ROI masks output
%                     Lck.Trial_ROIMask{iTrial,1} = mask;
%                     
%                     % all astrocyte ROI masks together
%                     trialmasks(:,:,iTrial) = mask;
%                     
%                 end
%             end
%             
% %             end
%         end
%         
%         % append peak data
%         dataNames=fieldnames(data);
%         data2= struct2cell(data);
%         data3= [data2{:}];
%         
%         AllData=vertcat(AllData, data3);
%         
%         clearvars data data2 data3
%         
%         % Mask of Active Pixels (proportional to number of trials)
%         sumImg =sum(trialmasks,3);
%         fracImg = sum(trialmasks,3)./numel(TrialNames);
%         
%         % calculate the score for responses (normalized to the "threshold"
%         thresh = 1/numel(TrialNames); % or 1/numel(trials)
%         activePxIdx = fracImg > 0;
%         activePx = sum(activePxIdx(:));
%         fracImg(fracImg <= thresh) = NaN;
%         score = nansum(nansum(fracImg)) / activePx;
%         
%         % save the mask of active pixels
%         FigFileName1 = fullfile(DatName,'Sum_Pixel_Mask.tif');
%         FigFileName2 = fullfile(DatName,'FracActive_Pixel_Mask.tif');
%         
%         f = figure('visible', 'off');
%         imagesc(sumImg);
%         axis square
%         caxis([0 1])
%         colormap('jet')
%         colorbar
%         saveas(gcf,FigFileName1{1,1})
%         close(f)
%         
%         
%         f = figure('visible', 'off');
%         imagesc(fracImg);
%         axis square
%         caxis([0 1])
%         colormap('jet')
%         colorbar
%         saveas(gcf,FigFileName2{1,1})
%         close(f)
%         
%         % store Lck field of view information
%         
%         for iLck= 1:length(TrialNames)
%             Lck.Total_ROIMask{iLck,1} = fracImg;
%             Lck.Response_Score{iLck,1}=score;
%         end
%         
%         LckNames=fieldnames(Lck);
%         Lckdata2= struct2cell(Lck);
%         Lckdata3= [Lckdata2{:}];
%         
%         LckData=vertcat(LckData, Lckdata3);
%         
%         clearvars Lck Lckdata2 Lckdata3
%         end
% 
% end
% 
% 
% % %% Save all data for R analysis
% % peak data
% AllData2= [dataNames';AllData];
% 
% %Field of View data
% LckFieldData= [LckNames';LckData];
% LckFieldData2=LckFieldData;
% LckFieldData(:,8:9)=[];
% 
% %onsetTimeTable (ROI table)
% names={'ROI','Trial','Channel','Spot','Animal', 'Condition','baseline',...
%     'trace','ROIIdx','PixelSize','FrameRate','OnsetTime','TraceAUC1','TraceAUC10'};
% All_traces2=vertcat(names, All_traces);
% All_traces2(:,8:9)=[];
% 
% 
% 
% cd(fullfile(Settings.ResultsFolder));
% 
% % write date to created file
% save(SaveFiles{1,1}, 'AllData2','-v7.3');
% cell2csv(SaveFiles{1,2}, AllData2);
% save(SaveFiles{1,3}, 'All_traces','-v7.3');
% cell2csv(SaveFiles{1,4}, All_traces2);
% save(SaveFiles{1,5}, 'LckFieldData2','-v7.3');
% cell2csv(SaveFiles{1,6}, LckFieldData);
% 
% close all; clear variables;
% 
