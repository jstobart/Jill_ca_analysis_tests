
clearvars
close all

% time windows based on stimulation
NOnsetWindow= 2; %1 for short stim % neuronal onset times
AOnsetWindow= 12; % 5 for short stim % astrocyte onset times
Fast_OnsetWindow=1.09;

NPTWindow= 9; % one second longer than stimulation for peak times
APTWindow= 15; % astrocyte longer than stimulation for peak times

stimwindow=30; % 5 s baseline, 15 s imaging

% %% Load peak data
%
% %load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\FilesforMatlab\Peaks_1stCohort_Lck_nostim_vs_longstim_12_2017.mat');
% load('D:\Data\GCaMP_RCaMP\Revision\Lck_GCaMP\FilesforMatlab\Control_untreated\Peaks_1stCohort_Lck_nostim_vs_longstim_12_2017.mat');
%
% peaks1=AllData2(2:end,:);
%
% load('D:\Data\GCaMP_RCaMP\Revision\Lck_GCaMP\FilesforMatlab\Control_untreated\Peaks_2ndCohort_Lck_nostim_vs_longstim_01_2018.mat');
%
% peaks2=AllData2(2:end,:);
%
% Lck_tracesPeaks=vertcat(peaks1, peaks2);
%
%
% % Peak info- ROiTypes, etc.
%
% for iROI=1:length(Lck_tracesPeaks)
%     %find ROITypes
%     N_str= strfind(Lck_tracesPeaks{iROI, 1},'N');
%     D_str= strfind(Lck_tracesPeaks{iROI, 1},'D'); %hand selected dendrite
%     r_str= strfind(Lck_tracesPeaks{iROI, 1},'r'); %FLIKA ROIs
%     EF_str= strfind(Lck_tracesPeaks{iROI, 1},'E');
%
%     if ~isempty(N_str)
%         Lck_tracesPeaks{iROI,23}='Neuron';
%     elseif ~isempty(D_str)
%         Lck_tracesPeaks{iROI,23}='Dendrite';
%     elseif ~isempty(r_str)
%         P_str= strfind(Lck_tracesPeaks{iROI, 3},'GCaMP');
%         if ~isempty(P_str)
%             Lck_tracesPeaks{iROI,23}='Process';
%         else
%             Lck_tracesPeaks{iROI,23}='Dendrite';
%         end
%     elseif ~isempty(EF_str)
%         Lck_tracesPeaks{iROI,23}='Endfeet';
%     end
%
%     % make new unique trial names
%     Lck_tracesPeaks{iROI,24}=strcat(Lck_tracesPeaks{iROI,13},'_',Lck_tracesPeaks{iROI,15},'_',Lck_tracesPeaks{iROI,12});
%     % unique ROI names
%     Lck_tracesPeaks{iROI,25}=strcat(Lck_tracesPeaks{iROI,13},'_',Lck_tracesPeaks{iROI,15},'_',Lck_tracesPeaks{iROI,12}, '_',Lck_tracesPeaks{iROI,10});
% end
%
% for iROI=1:length(Lck_tracesPeaks)
%     %if a process ROI is overlapping with a soma or process, exclude it
%     %from data
%     nonOverlapIdx(iROI)=~ischar(Lck_tracesPeaks{iROI,19});
% end
%
% % remove overlapping processes
% Lck_tracesPeaks = Lck_tracesPeaks(nonOverlapIdx',:);


%% load trace data

%AstrocyteExcelFile='E:\Data\Two_Photon_Data\GCaMP_RCaMP\Manuscript\Figures\DataForTraceOverlay_Heatmaps\cytovslck_Nostimvsstim_comparison\cyto-AstrocyteTraces_Stim.xlsx';
%NeuronalExcelFile ='E:\Data\Two_Photon_Data\GCaMP_RCaMP\Manuscript\Figures\DataForTraceOverlay_Heatmaps\cytovslck_Nostimvsstim_comparison\cyto-NeuronalTraces_Stim.xlsx';

AstrocyteExcelFile='D:\Data\GCaMP_RCaMP\Manuscript\Neuron_2017\Revision\Figures\Lck-AstrocyteTraces_Stim.xlsx';
NeuronalExcelFile ='D:\Data\GCaMP_RCaMP\Manuscript\Neuron_2017\Revision\Figures\Lck-NeuronalTraces_Stim.xlsx';

% for the first cohort

% Load trace data and run onset time calculation
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\\FilesforMatlab\Traces_1stCohort_Lck_nostim_vs_longstim_12_2017.mat');
load('D:\Data\GCaMP_RCaMP\Revision\Lck_GCaMP\FilesforMatlab\Control_untreated\Traces_1stCohort_Lck_nostim_vs_longstim_12_2017.mat');

cohort1=All_traces;

load('D:\Data\GCaMP_RCaMP\Revision\Lck_GCaMP\FilesforMatlab\Control_untreated\Traces_2ndCohort_Lck_nostim_vs_longstim_01_2018.mat');

All_traces(:,7)=[];
All_traces(:,14)=[];
All_traces(:,14)=[];
All_traces(:,14)=[];

Shortstim=vertcat(cohort1, All_traces);


% trace info- ROIType etc.
% get info for plotting

% ONLY CONSIDER STIM DATA
% get rid of astrocyte neuropil traces
for xROI=1:size(Shortstim,1)
    NostimIdx(xROI)= strcmp(Shortstim{xROI,6},'Nostim');
end
Shortstim=Shortstim(~NostimIdx',:);


% get rid of handclicked neuropil ROIs because they are not relevant if
% FLIKA ROIs are included for the neuron channel

% get rid of astrocyte neuropil traces
for xROI=1:size(Shortstim,1)
    GCNP2(xROI)= strcmp(Shortstim{xROI,1},'np');
end
Shortstim=Shortstim(~GCNP2',:);

for iROI=1:length(Shortstim)
    %find ROITypes
    N_str= strfind(Shortstim{iROI, 1},'N');
    D_str= strfind(Shortstim{iROI, 1},'D'); %hand selected dendrite
    r_str= strfind(Shortstim{iROI, 1},'r'); %FLIKA ROIs
    EF_str= strfind(Shortstim{iROI, 1},'E');
    
    if ~isempty(N_str)
        Shortstim{iROI,14}='Neuron';
    elseif ~isempty(D_str)
        Shortstim{iROI,14}='Dendrite';
    elseif ~isempty(r_str)
        P_str= strfind(Shortstim{iROI, 3},'GCaMP');
        if ~isempty(P_str)
            Shortstim{iROI,14}='Process';
        else
            Shortstim{iROI,14}='Dendrite';
        end
    elseif ~isempty(EF_str)
        Shortstim{iROI,14}='Endfeet';
    end
    
    % make new unique trial names
    Shortstim{iROI,15}=strcat(Shortstim{iROI,5},'_',Shortstim{iROI,4},'_',Shortstim{iROI,2});
    % unique ROI names
    Shortstim{iROI,16}=strcat(Shortstim{iROI,5},'_',Shortstim{iROI,4},'_',Shortstim{iROI,2}, '_',Shortstim{iROI,1});
end

for iROI=1:length(Shortstim)
    %if a process ROI is overlapping with a soma or process, exclude it
    %from data
    nonOverlapIdx2(iROI)=~ischar(Shortstim{iROI,12});
end

% remove overlapping processes
Shortstim = Shortstim(nonOverlapIdx2',:);



% Calculate the first peak onset time and AUC after stim

% peak onsets and AUC in the first second after stim for each ROI
for iROI= 1:length(Shortstim)
    
    trace=Shortstim{iROI,9};
    FrameRate=Shortstim{iROI,13};
    nframes=length(trace);
    TimeX(1:nframes) = (1:nframes)/FrameRate;
    BL_time=Shortstim{iROI, 8}/Shortstim{iROI,13};
    baselineCorrectedTime=TimeX-BL_time;
    %first 1 sec after stim onset
    x1=Shortstim{iROI, 8};
    x2=round((BL_time+1)*FrameRate);
    x3= round(FrameRate*(BL_time+10));
    
    % onset time
    Onsets=find_first_onset_time(baselineCorrectedTime(10:end), trace(10:end),2.5,2);
    if isempty(Onsets)
        Onsets=nan(1,1);
    end
    Shortstim{iROI, 17}= Onsets;
    % AUC
    Shortstim{iROI,18}=trapz(trace(x1:x2));
    Shortstim{iROI,19}=trapz(trace(x1:x3));
    
end


%% separate out IP3 mice from wildtypes

for xROI=1:size(Shortstim,1)
    IP3Match(xROI)= ~isempty(strfind(Shortstim{xROI,5},'IP'));
end
IP3=Shortstim(IP3Match',:);
Lck_traces=Shortstim(~IP3Match',:);


% genotypes in IP3
for iROI=1:length(IP3)
    %find controls
    str2= strfind(IP3{iROI, 5},'2');
    str3= strfind(IP3{iROI, 5},'3');
    str6= strfind(IP3{iROI, 5},'6');
    
    if ~isempty(str2)
        IP3{iROI,20}='WT';
    elseif ~isempty(str3)
        IP3{iROI,20}='WT';
    elseif ~isempty(str6)
        IP3{iROI,20}='WT';
    else
        IP3{iROI,20}='KO';
    end    
end


%% Mean Traces- Lck Control Data

% ROI with a response to stimulation

nframes2=round(stimwindow*FrameRate);

for iROI=1:length(Lck_traces)
    rc_str(iROI)= ~isempty(strfind(Lck_traces{iROI,3},'RCaMP'));
    gc_str(iROI)= ~isempty(strfind(Lck_traces{iROI,3},'GCaMP'));
end

RCaMP=Lck_traces(rc_str',:);
GCaMP=Lck_traces(gc_str',:);

OT_RCaMP_traces=[];
for iROI=1:length(RCaMP)
    respOTIdx1(iROI)=~isempty(find((RCaMP{iROI,17}>0 && RCaMP{iROI,17}<=Fast_OnsetWindow),1));
    tempY= RCaMP{iROI,9};
    BL_time=RCaMP{iROI,8}/RCaMP{iROI,13};
    if length(tempY)>(nframes2+59)
        if BL_time<6
            tempY2=tempY(1:nframes2); %(stimwindow*FrameRate));
        else
            tempY2=tempY(59:(nframes2+58));
        end
        if respOTIdx1(iROI)
            OT_RCaMP_traces= horzcat(OT_RCaMP_traces, tempY2);
        end
    end
end
RrespOT=RCaMP(respOTIdx1',:); % responding neurons
OT_RCaMP_mean= nanmean(OT_RCaMP_traces,2);

fastAC_traces=[];
slowAC_traces=[];
for iROI=1:length(GCaMP)
    fast_respOTIdx2(iROI)=~isempty(find((GCaMP{iROI,17}>0 && GCaMP{iROI,17}<=Fast_OnsetWindow),1));
    delayed_respOTIdx2(iROI)=~isempty(find((GCaMP{iROI,17}>Fast_OnsetWindow && GCaMP{iROI,17}<=AOnsetWindow),1));
    tempY= GCaMP{iROI,9};
    BL_time=GCaMP{iROI, 8}/GCaMP{iROI,13};
    if length(tempY)>(nframes2+59)
        if BL_time<6
            tempY2=tempY(1:nframes2); %(stimwindow*FrameRate));
        else
            tempY2=tempY(59:(nframes2+58));
        end
        if fast_respOTIdx2(iROI)
            fastAC_traces= horzcat(fastAC_traces, tempY2);
        end
        if delayed_respOTIdx2(iROI)
            slowAC_traces= horzcat(slowAC_traces, tempY2);
        end
    end
end
fastAC=GCaMP(fast_respOTIdx2',:); % fast responding astrocytes
delayedAC=GCaMP(delayed_respOTIdx2',:); % delayed responding astrocytes

GrespOT =vertcat(fastAC,delayedAC);
% mean trace of astrocytes
fastAC_mean= nanmean(fastAC_traces,2);
slowAC_mean= nanmean(slowAC_traces,2);



%% shaded error bar with
green=[(27/255) (120/255) (55/255)];
purple=[(123/255) (50/255) (148/255)];
blue= [(0/255) (114/255) (178/255)];

TimeX2(1:nframes2) = (1:nframes2)/FrameRate;

% SEM calculations
fastAC_SDTrace = nanstd(fastAC_traces');
fastAC_SEM=fastAC_SDTrace/sqrt(size(fastAC_traces,2));

slowAC_SDTrace = nanstd(slowAC_traces');
slowAC_SEM=slowAC_SDTrace/sqrt(size(slowAC_traces,2));

RC_SDTrace = nanstd(OT_RCaMP_traces');
RC_SEM=RC_SDTrace/sqrt(size(OT_RCaMP_traces,2));

figure('name', 'Lck means- plus SEM')
hold on
axis off
xlim([-1 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

%ylim([-0.2 3]);
lineProps.col = {blue};
mseb(TimeX2,slowAC_mean',slowAC_SEM,lineProps)
lineProps.col = {green};
mseb(TimeX2,(fastAC_mean'+0.75),fastAC_SEM,lineProps)
lineProps.col = {purple};
mseb(TimeX2,(OT_RCaMP_mean'+2),RC_SEM,lineProps)

% H2=shadedErrorBar(TimeX,(fastAC_mean+0.75),fastAC_SEM)%,green,1);
% H3=shadedErrorBar(TimeX,(OT_RCaMP_mean+2),RC_SEM)%,purple,1);
rectangle('Position', [5 -0.3 8 4])
plot([-0.1 -0.1],[0 1], 'k','LineWidth', 1)


% figure('name', 'Lck short stim all means- plus SD')
% hold on
% axis off
% xlim([-1 25]);
% lineProps.width = 1;
% lineProps.edgestyle = ':';
% %ylim([-0.2 3]);
%
% lineProps.col = {blue};
% mseb(TimeX2,slowAC_mean',slowAC_SDTrace,lineProps)
% lineProps.col = {green};
% mseb(TimeX2,(fastAC_mean'+3),fastAC_SDTrace,lineProps)
% lineProps.col = {purple};
% mseb(TimeX2,(OT_RCaMP_mean'+7),RC_SDTrace,lineProps)
% rectangle('Position', [5 -0.3 8 10])
% plot([-0.1 -0.1],[0 1], 'k','LineWidth', 1)


%% traces for fast endfeet and fast processes

% Group Responders by ROIType
Resp_EF_traces=[];
Resp_processes_traces=[];

%find ROITypes
for xROI= 1:length(fastAC)
    EF_str= strfind(fastAC{xROI, 14},'Endfeet');
    P_str= strfind(fastAC{xROI, 14},'Process');
    
    tempY= fastAC{xROI,9};
    BL_time=fastAC{xROI, 8}/fastAC{xROI,13};
    if length(tempY)>(nframes2+59)
        if BL_time<6
            tempY2=tempY(1:nframes2); %(stimwindow*FrameRate));
        else
            tempY2=tempY(59:(nframes2+58));
        end
        if ~isempty(EF_str)
            Resp_EF_traces= horzcat(Resp_EF_traces, tempY2);
        elseif ~isempty(P_str)
            Resp_processes_traces= horzcat(Resp_processes_traces, tempY2);
        end
    end
end

% means and SDs

%endfeet
RespEFmeanTrace = nanmean(Resp_EF_traces,2)';
RespEFSEMTrace = nanstd(Resp_EF_traces')/sqrt(size(Resp_EF_traces,2));

%processes
RespPmeanTrace = nanmean(Resp_processes_traces,2)';
RespPSEMTrace = nanstd(Resp_processes_traces')/sqrt(size(Resp_processes_traces,2));


figure('name', 'Lck FAST endfeet vs processes- plus SEM')
hold on
axis off
xlim([-1 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

%ylim([-0.2 3]);
lineProps.col = {blue};
mseb(TimeX2,RespEFmeanTrace,RespEFSEMTrace,lineProps)
lineProps.col = {green};
mseb(TimeX2,(RespPmeanTrace+2.5),RespPSEMTrace,lineProps)
rectangle('Position', [5 -1 8 5.5])
plot([-0.1 -0.1],[0 1], 'k','LineWidth', 1)


%% Data for heat maps

% sort astrocytes by onset time
[~, OT_GCidx] = sort([GrespOT{:,17}], 'ascend');
GrespOT_sort=GrespOT(OT_GCidx,:);

% sort neuronss by onset time
[~, OT_RCidx] = sort([RrespOT{:,17}], 'ascend');
RrespOT_sort=RrespOT(OT_RCidx,:);


nframes3=round(20*FrameRate);


% astrocyte table
AC_traces=[];
for xROI= 1:length(GrespOT_sort)
    tempOnset = GrespOT_sort{xROI,17}; % individual onset time
    tempY = GrespOT_sort{xROI,9}; % individual trace
    BL_time=GrespOT_sort{xROI, 8}/GrespOT_sort{xROI,13};
    if length(tempY)>(nframes3+59)
        if BL_time<6
            tempY2=tempY(1:nframes3); %(stimwindow*FrameRate));
        else
            tempY2=tempY(59:(nframes3+58));
        end
    end
    temp=vertcat(tempOnset,tempY2);
    AC_traces=vertcat(AC_traces,temp');
end

% neuronal table
N_traces=[];
for xROI= 1:length(RrespOT_sort)
    tempOnset = RrespOT_sort{xROI,17}; % individual onset time
    tempY = RrespOT_sort{xROI,9}; % individual trace
    
    BL_time=RrespOT_sort{xROI, 8}/RrespOT_sort{xROI,13};
    if length(tempY)>(nframes3+59)
        if BL_time<6
            tempY2=tempY(1:nframes3); %(stimwindow*FrameRate));
        else
            tempY2=tempY(59:(nframes3+58));
        end
    end
    temp=vertcat(tempOnset,tempY2);
    N_traces=vertcat(N_traces,temp');
end


xlswrite(AstrocyteExcelFile, AC_traces)
xlswrite(NeuronalExcelFile, N_traces)


%% Mean Traces- IP3 knockouts vs littermates

% ROI with a response to stimulation

nframes2=round(stimwindow*FrameRate);

for iROI=1:length(IP3)
    ko_str(iROI)= ~isempty(strfind(IP3{iROI,20},'KO'));
    wt_str(iROI)= ~isempty(strfind(IP3{iROI,20},'WT'));
end

IP3KO=IP3(ko_str',:);
IP3WT=IP3(wt_str',:);

% only consider GCaMP traces!
for iROI=1:length(IP3KO)
    gc_str3(iROI)= ~isempty(strfind(IP3KO{iROI,3},'GCaMP'));
end

IP3KO_GC=IP3KO(gc_str3',:);

for iROI=1:length(IP3WT)
    gc_str2(iROI)= ~isempty(strfind(IP3WT{iROI,3},'GCaMP'));
end

IP3WT_GC=IP3WT(gc_str2',:);

%knockouts
fastKO_traces=[];
slowKO_traces=[];
for iROI=1:length(IP3KO_GC)
    fast_KOIdx(iROI)=~isempty(find((IP3KO_GC{iROI,17}>0 && IP3KO_GC{iROI,17}<=Fast_OnsetWindow),1));
    slow_KOIdx(iROI)=~isempty(find((IP3KO_GC{iROI,17}>Fast_OnsetWindow && IP3KO_GC{iROI,17}<=AOnsetWindow),1));
    tempY= IP3KO_GC{iROI,9};
    BL_time=IP3KO_GC{iROI, 8}/IP3KO_GC{iROI,13};
    if length(tempY)>(nframes2+59)
        if BL_time<6
            tempY2=tempY(1:nframes2); %(stimwindow*FrameRate));
        else
            tempY2=tempY(59:(nframes2+58));
        end
        if fast_KOIdx(iROI)
            fastKO_traces= horzcat(fastKO_traces, tempY2);
        end
        if slow_KOIdx(iROI)
            slowKO_traces= horzcat(slowKO_traces, tempY2);
        end
    end
end
fastKO=IP3KO_GC(fast_KOIdx',:); % fast responding astrocytes
delayedKO=IP3KO_GC(slow_KOIdx',:); % delayed responding astrocytes

respKO =vertcat(fastKO,delayedKO);

% mean trace of KO
fastKO_mean= nanmean(fastKO_traces,2);
slowKO_mean= nanmean(slowKO_traces,2);


%wildtypes
fastWT_traces=[];
slowWT_traces=[];
for iROI=1:length(IP3WT_GC)
    fast_WTIdx(iROI)=~isempty(find((IP3WT_GC{iROI,17}>0 && IP3WT_GC{iROI,17}<=Fast_OnsetWindow),1));
    slow_WTIdx(iROI)=~isempty(find((IP3WT_GC{iROI,17}>Fast_OnsetWindow && IP3WT_GC{iROI,17}<=AOnsetWindow),1));
    tempY= IP3WT_GC{iROI,9};
    BL_time=IP3WT_GC{iROI, 8}/IP3WT_GC{iROI,13};
    if length(tempY)>(nframes2+59)
        if BL_time<6
            tempY2=tempY(1:nframes2); %(stimwindow*FrameRate));
        else
            tempY2=tempY(59:(nframes2+58));
        end
        if fast_WTIdx(iROI)
            fastWT_traces= horzcat(fastWT_traces, tempY2);
        end
        if slow_WTIdx(iROI)
            slowWT_traces= horzcat(slowWT_traces, tempY2);
        end
    end
end
fastWT=IP3WT_GC(fast_WTIdx',:); % fast responding astrocytes
delayedWT=IP3WT_GC(slow_WTIdx',:); % delayed responding astrocytes

respWT =vertcat(fastWT,delayedWT);

% mean trace of WT
fastWT_mean= nanmean(fastWT_traces,2);
slowWT_mean= nanmean(slowWT_traces,2);


%% shaded error bar with
green=[(27/255) (120/255) (55/255)];
blue= [(0/255) (114/255) (178/255)];

TimeX2(1:nframes2) = (1:nframes2)/FrameRate;

% SEM calculations
fastKO_SDTrace = nanstd(fastKO_traces');
fastKO_SEM=fastKO_SDTrace/sqrt(size(fastKO_traces,2));

slowKO_SDTrace = nanstd(slowKO_traces');
slowKO_SEM=slowKO_SDTrace/sqrt(size(slowKO_traces,2));

fastWT_SDTrace = nanstd(fastWT_traces');
fastWT_SEM=fastWT_SDTrace/sqrt(size(fastWT_traces,2));

slowWT_SDTrace = nanstd(slowWT_traces');
slowWT_SEM=slowWT_SDTrace/sqrt(size(slowWT_traces,2));


%1
figure('name', 'Knockouts fast vs delayed- plus SEM')
hold on
axis off
xlim([-1 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

%ylim([-0.2 3]);
lineProps.col = {blue};
mseb(TimeX2,slowKO_mean',slowKO_SEM,lineProps)
lineProps.col = {green};
mseb(TimeX2,(fastKO_mean'+0.75),fastKO_SEM,lineProps)

% H2=shadedErrorBar(TimeX,(fastAC_mean+0.75),fastAC_SEM)%,green,1);
% H3=shadedErrorBar(TimeX,(OT_RCaMP_mean+2),RC_SEM)%,purple,1);
rectangle('Position', [5 -0.3 8 4])
plot([-0.1 -0.1],[0 1], 'k','LineWidth', 1)


%2
figure('name', 'Wildtypes fast vs delayed- plus SEM')
hold on
axis off
xlim([-1 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

%ylim([-0.2 3]);
lineProps.col = {blue};
mseb(TimeX2,slowWT_mean',slowWT_SEM,lineProps)
lineProps.col = {green};
mseb(TimeX2,(fastWT_mean'+1),fastWT_SEM,lineProps)

% H2=shadedErrorBar(TimeX,(fastAC_mean+0.75),fastAC_SEM)%,green,1);
% H3=shadedErrorBar(TimeX,(OT_RCaMP_mean+2),RC_SEM)%,purple,1);
rectangle('Position', [5 -0.3 8 4])
plot([-0.1 -0.1],[0 1], 'k','LineWidth', 1)


%3
figure('name', 'Knockout fast vs wildtype fast')
hold on
axis off
xlim([-1 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

%ylim([-0.2 3]);
lineProps.col = {blue};
mseb(TimeX2,fastKO_mean',fastKO_SEM,lineProps)
lineProps.col = {green};
mseb(TimeX2,(fastWT_mean'+1),fastWT_SEM,lineProps)

rectangle('Position', [5 -0.3 8 4])
plot([-0.1 -0.1],[0 1], 'k','LineWidth', 1)

%4
figure('name', 'Knockout delayed vs wildtype delayed')
hold on
axis off
xlim([-1 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

%ylim([-0.2 3]);
lineProps.col = {blue};
mseb(TimeX2,slowKO_mean',slowKO_SEM,lineProps)
lineProps.col = {green};
mseb(TimeX2,(slowWT_mean'+1),slowWT_SEM,lineProps)

rectangle('Position', [5 -0.3 8 4])
plot([-0.1 -0.1],[0 1], 'k','LineWidth', 1)


%5
figure('name', 'Knockout vs wildtype, fast vs delayed')
hold on
axis off
xlim([-1 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

lineProps.col = {blue};
mseb(TimeX2,(fastKO_mean'+2),fastKO_SEM,lineProps)
lineProps.col = {green};
mseb(TimeX2,(fastWT_mean'+3),fastWT_SEM,lineProps)
lineProps.col = {blue};
mseb(TimeX2,slowKO_mean',slowKO_SEM,lineProps)
lineProps.col = {green};
mseb(TimeX2,(slowWT_mean'+1),slowWT_SEM,lineProps)

rectangle('Position', [5 -0.3 8 4])
plot([-0.1 -0.1],[0 1], 'k','LineWidth', 1)






