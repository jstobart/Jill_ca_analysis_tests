%% figure- no stimulation vs. stimulation- exact ROI pairs
% AFTER BARIUM INTRACORTICAL INJECTION

%function [AC_ROIs,AP_ROIs,BR_ROIs,EF_ROIs] = Raw_LineGraphs(DirName)
% [AC_ROIs,AP_ROIs,BR_ROIs,EF_ROIs] = Raw_LineGraphs(DirName)
DirName = 'E:\Data\Two_Photon_Data\GCaMP6\Intracortical\spots\Results\2013_10_10_spot1\GC603_2013_10_10_10_13\bhv2PROI.mat';
load(DirName)


%%
 %calculate the sum of each channel for each trial and make a cell array 

    for iROI = 1:length(bhvROI.ROIFname{1,1}) % ROI number   
    sumtrials{iROI} = bhvROI.ROIch{1,1}{1,iROI} + bhvROI.ROIch{1,1}{2,iROI}; %sum
    end

%% Determine the baseline length (# of frames) and imaging length (# of frames) of each trial
   
%Determine the baseline length (# of frames) and imaging length (# of frames) of each trial
%100 frame BASELINE!!!!
Baseline_frames= 100;
Imaging_frames= length(bhvROI.ROIch{1,1}{1,1}) - Baseline_frames;

%Calculate the Baseline change in fluorescence with average and SD
%trial number= 1
    for iROI = 1:length(bhvROI.ROIFname{1,1}) % ROI number   
    %average of baseline values
    F0{iROI} = mean(sumtrials{iROI}(1:Baseline_frames));
       
    %Imaging Values
    F1{iROI} = sumtrials{iROI}((Baseline_frames+1):(Baseline_frames+Imaging_frames));
    
     % calculate the change in fluorescence (deltaF) of the peak_values
      deltaF{iROI} = ((F1{iROI}(1:Imaging_frames) - F0{iROI})/F0{iROI})* 100;
           
    end        


%% Plots of each ROI 

 figure('name','Trial3_BariumInjection_55PSI_1min_ACROIs','numbertitle','off')
  %trace over all trials
  
    tvec= (1:1:length(deltaF{1,1}))/bhvROI.ExcelFrameRate(1,1); 
    
    hold on
   
    plot(tvec,deltaF{1,1}, 'k','LineWidth', 2)
    plot(tvec,(deltaF{1,2}+100), 'g','LineWidth', 2)
    plot(tvec,deltaF{1,3}+200, 'b','LineWidth', 2)
    plot(tvec,deltaF{1,4}+300, 'r','LineWidth', 2)
    plot(tvec,deltaF{1,5}+400, 'y','LineWidth', 2)
    plot(tvec,deltaF{1,6}+500, 'm','LineWidth', 2)
    plot(tvec,deltaF{1,7}+600, 'c','LineWidth', 2)
    
    plot(tvec,deltaF{1,8}+700, 'k','LineWidth', 2) 
    plot(tvec,deltaF{1,9}+800, 'g','LineWidth', 2)
    plot(tvec,deltaF{1,10}+900, 'b','LineWidth', 2)
    plot(tvec,deltaF{1,11}+1000, 'r','LineWidth', 2)
    plot(tvec,deltaF{1,12}+1100, 'y','LineWidth', 2)
    plot(tvec,deltaF{1,13}+1200, 'm','LineWidth', 2)
    plot(tvec,deltaF{1,14}+1300, 'c','LineWidth', 2)
    plot(tvec,deltaF{1,15}+1400, 'k','LineWidth', 2)
    plot(tvec,deltaF{1,16}+1500, 'g','LineWidth', 2)
    plot(tvec,deltaF{1,17}+1600, 'b','LineWidth', 2)
    plot(tvec,deltaF{1,18}+1700, 'r','LineWidth', 2)
   
    plot([-30 -30],[0 50],'LineWidth', 2)
    %Injection line- 1 min
        hold on
        plot([0 88],[-50 -50], 'm','LineWidth', 4)

    set(gca,'XTick',0:20:140)
    set(gca, 'YTick', [])
    xlabel('Time (sec)')
    grid on
    
  
  legend('AC1','AC2','AC3','AC4','AC5','AC6','AC7','AC8','AC9','AC10','AC11','AC12','AC13','AC14','AC15','AC16','AC17','AC18')
 
%% AP ROIS
 figure('name','Trial3_BariumInjection_55PSI_1min_APROIs','numbertitle','off')
  %trace over all trials
  
    tvec= (1:1:length(deltaF{1,1}))/bhvROI.ExcelFrameRate(1,1); 
    
    hold on
   
    plot(tvec,deltaF{1,19}, 'k','LineWidth', 2)
    plot(tvec,(deltaF{1,20}+100), 'g','LineWidth', 2)
    plot(tvec,deltaF{1,21}+200, 'b','LineWidth', 2)
    plot(tvec,deltaF{1,22}+300, 'r','LineWidth', 2)
    plot(tvec,deltaF{1,23}+400, 'y','LineWidth', 2)
    plot(tvec,deltaF{1,24}+500, 'm','LineWidth', 2)
    plot(tvec,deltaF{1,25}+600, 'c','LineWidth', 2)
    
    plot(tvec,deltaF{1,26}+700, 'k','LineWidth', 2) 
    plot(tvec,deltaF{1,27}+800, 'g','LineWidth', 2)
    plot(tvec,deltaF{1,28}+900, 'b','LineWidth', 2)
    plot(tvec,deltaF{1,29}+1000, 'r','LineWidth', 2)
    plot(tvec,deltaF{1,30}+1100, 'y','LineWidth', 2)
    plot(tvec,deltaF{1,31}+1200, 'm','LineWidth', 2)
    plot(tvec,deltaF{1,32}+1300, 'c','LineWidth', 2)
    plot(tvec,deltaF{1,33}+1400, 'k','LineWidth', 2)
    
   
    plot([-30 -30],[0 50],'LineWidth', 2)
    %Injection line- 1 min
        hold on
        plot([0 88],[-50 -50], 'm','LineWidth', 4)

    set(gca,'XTick',0:20:140)
    set(gca, 'YTick', [])
    xlabel('Time (sec)')
    grid on
    
  
  legend('AP1','AP2','AP3','AP4','AP5','AP6','AP7','AP8','AP9','AP10','AP11','AP12','AP13','AP14','AP15')
 
%% EF ROIS
 figure('name','Trial3_BariumInjection_55PSI_1min_EF&LargeROIs','numbertitle','off')
  %trace over all trials
  
    tvec= (1:1:length(deltaF{1,1}))/bhvROI.ExcelFrameRate(1,1); 
    
    hold on
   
    plot(tvec,deltaF{1,36}, 'k','LineWidth', 2)
    plot(tvec,(deltaF{1,34}+100), 'g','LineWidth', 2)
    plot(tvec,deltaF{1,35}+200, 'b','LineWidth', 2)
   
    plot([-30 -30],[0 50],'LineWidth', 2)
    %Injection line- 1 min
        hold on
        plot([0 88],[-50 -50], 'm','LineWidth', 4)

    set(gca,'XTick',0:20:140)
    set(gca, 'YTick', [])
    xlabel('Time (sec)')
    grid on
    
  
  legend('G','EF1','EF2')
 
 

