%% figure- no stimulation vs. stimulation- exact ROI pairs


%function [AC_ROIs,AP_ROIs,BR_ROIs,EF_ROIs] = Raw_LineGraphs(DirName)
% [AC_ROIs,AP_ROIs,BR_ROIs,EF_ROIs] = Raw_LineGraphs(DirName)
DirName = 'E:\Data\Two_Photon_Data\GCaMP6\GC605\spots\Results\2013_6_11_spot1\GC605_2013_06_11_16_29\bhv2PROI.mat';
load(DirName)


%%
 %calculate the sum of each channel for each trial and make a cell array 
for itrial = 1:length(bhvROI.ROIFname)  %trial number
    for iROI = 1:length(bhvROI.ROIFname{1,1}) % ROI number   
    sumtrials{iROI}{itrial} = bhvROI.ROIch{1,itrial}{1,iROI} + bhvROI.ROIch{1,itrial}{2,iROI}; %sum
    end
end
%% Determine the baseline length (# of frames) and imaging length (# of frames) of each trial
   
%Determine the baseline length (# of frames) and imaging length (# of frames) of each trial
%20 SEC BASELINE!!!!
Baseline= (20* bhvROI.ExcelFrameRate(1,1));
Baseline_frames= round(Baseline*10^0)/10^0;
Imaging_frames= length(bhvROI.ROIch{1,1}{1,1}) - Baseline_frames;

%Calculate the Baseline change in fluorescence with average and SD
for itrial = 1:length(bhvROI.ROIFname)  %trial number
    for iROI = 1:length(bhvROI.ROIFname{1,1}) % ROI number   
    %average of baseline values
    F0{itrial}{iROI} = mean(sumtrials{iROI}{itrial}(1:Baseline_frames));
       
    %Imaging Values
    F1{itrial}{iROI} = sumtrials{iROI}{itrial}((Baseline_frames+1):(Baseline_frames+Imaging_frames));
    
     % calculate the change in fluorescence (deltaF) of the peak_values
      deltaF{iROI}{itrial} = ((F1{itrial}{iROI}(1:Imaging_frames) - F0{itrial}{iROI})/F0{itrial}{iROI})* 100;
           
    end        
end
%%
for itrial = 1:length(bhvROI.ROIFname)  %trial number
    for iROI = 1:length(bhvROI.ROIFname{1,1}) % ROI number 
       xtrial = length(bhvROI.ROIFname);
%calculate the average of all trials at each time point and add to structure
reshape_deltaF= (cellfun(@cell2mat,deltaF,'UniformOutput',0)); 
d{iROI} = reshape(reshape_deltaF{iROI}',Imaging_frames,xtrial);
timeseriesAVG{iROI}= mean(d{iROI}');

    end
end
%% Plots of each trial and average
for iROI = 15:20
 figure('name','AP ROIs')
 %No stim   
 %all 10 trials in grey
 for itrial = 1:length(bhvROI.ROIFname)
     hold on
     trial_deltaF = deltaF{iROI}{itrial};
     tvec= (1:1:length(deltaF{iROI}{itrial}))/bhvROI.ExcelFrameRate(1,1); 
     grey = [0.4,0.4,0.4];
     plot (tvec,trial_deltaF, 'color', grey)
     
   %average of trials in red
   plot(tvec,timeseriesAVG{iROI}, 'r', 'LineWidth', 2)
   set(gca,'Xlim',[0 200],'YLim',[-50 200],'XTick', 0:60:200, 'YTick', -50:50:200)
   xlabel('Time (s)')
   ylabel('% DF/Fo')
   
 end
end
 
 

