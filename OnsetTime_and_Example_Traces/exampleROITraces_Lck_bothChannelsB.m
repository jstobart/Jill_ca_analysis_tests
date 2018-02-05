clearvars
close all

% time windows based on stimulation
NOnsetWindow= 2; %1 for short stim % neuronal onset times
AOnsetWindow= 12; % 5 for short stim % astrocyte onset times
Fast_OnsetWindow=1.09;

NPTWindow= 9; % one second longer than stimulation for peak times
APTWindow= 15; % astrocyte longer than stimulation for peak times

stimwindow=30; % 5 s baseline, 15 s imaging
nframes=round(30*11.8371);

%% load control trace data

load('D:\Data\GCaMP_RCaMP\Revision\Lck_GCaMP\FilesforMatlab\Control_untreated\Traces_2ndCohort_Lck_nostim_vs_longstim_01_2018.mat');

All_traces(:,7)=[];
All_traces(:,14)=[];
All_traces(:,14)=[];
All_traces(:,14)=[];

Shortstim=All_traces;

TimeX(1:nframes) = (1:nframes)/Shortstim{iROI,13};
BL_time=Shortstim{iROI, 8}/Shortstim{iROI,13};
baselineCorrectedTime=TimeX-BL_time;

% Extract info for the required mice
for iROI=1:length(Shortstim)
    str4(iROI)= ~isempty(strfind(Shortstim{iROI,5},'WT_LR4'));
end

WT_LR4=Shortstim(str4',:);

%% Control

for iROI=1:length(WT_LR4)
    spt3B_str(iROI)= ~isempty(strfind(WT_LR4{iROI,4},'spot3'));
end
WT_LR4_spot3=WT_LR4(spt3B_str',:);

for iROI=1:length(WT_LR4_spot3)
    rcamp_str(iROI)= ~isempty(strfind(WT_LR4_spot3{iROI,3},'RCaMP'));
    gcamp_str(iROI)= ~isempty(strfind(WT_LR4_spot3{iROI,3},'GCaMP'));
end
WT_LR4_spot3_RC=WT_LR4_spot3(rcamp_str',:);
WT_LR4_spot3_GC=WT_LR4_spot3(gcamp_str',:);

% trial 1
for iROI=1:length(WT_LR4_spot3_RC)
    trial_str(iROI)= ~isempty(strfind(WT_LR4_spot3_RC{iROI,2},'trial01'));
end
trial8=WT_LR4_spot3_RC(trial_str',:);

figure('name', 'responding RCaMP ROIs- WT_LR4, Spot3 Trial1')
hold on
axis off
for ii=1:length(trial8)
    %Responding(ii)=~isempty(find((trial8{ii,14}>0 && trial8{ii,14}<=Fast_OnsetWindow),1));
    tempY1=smooth(trial8{ii,9},3);
    tempY1=tempY1(1:nframes);
    %if Responding(ii)
        plot(baselineCorrectedTime(1:nframes),tempY1'+(3*(ii-1)),'r')%'LineWidth',1);
    %end
    
end

for iROI=1:length(WT_LR4_spot3_GC)
    trial_str(iROI)= ~isempty(strfind(WT_LR4_spot3_GC{iROI,2},'trial01'));
end
trial1_GC=WT_LR4_spot3_GC(trial_str',:);

figure('name', 'responding GCaMP ROIs- WT_LR4, Spot3 Trial1')
hold on
axis off
for ii=1:length(trial1_GC)
    %Responding(ii)=~isempty(find((trial8{ii,14}>0 && trial8{ii,14}<=Fast_OnsetWindow),1));
    tempY1=smooth(trial1_GC{ii,9},3);
    tempY1=tempY1(1:nframes);
    %if Responding(ii)
        plot(baselineCorrectedTime(1:nframes),tempY1'+(3*(ii-1)),'r')%'LineWidth',1);
    %end
    
end
%%
figure('name','GCaMP and RCaMP examples')
hold on
axis off
for ii=1:size(shortstim_ACtraces,2)
    tempY1=smooth(shortstim_ACtraces(:,ii),3);
    tempY1=tempY1(1:stimwindow);
plot(TimeX(1:stimwindow),tempY1'+(3*(ii-1)),'g')%'LineWidth',1);
   
end
for ii=1:size(shortstim_Ntraces,2)
    tempY2=smooth(shortstim_Ntraces(:,ii),3);
   tempY2=tempY2(1:stimwindow);

plot(TimeX(1:stimwindow),tempY2'+(3*(ii-1)),'r')
rectangle('Position', [5 -1 8 20])
%plot([-1 -1],[-1 1], 'k','LineWidth', 1)
%plot([-5 -1],[-1 -1], 'k','LineWidth', 1)
end