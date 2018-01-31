
clearvars
close all

% time windows based on stimulation
NOnsetWindow= 2; %1 for short stim % neuronal onset times
AOnsetWindow= 12; % 5 for short stim % astrocyte onset times
Fast_OnsetWindow=1.09;

NPTWindow= 9; % one second longer than stimulation for peak times
APTWindow= 15; % astrocyte longer than stimulation for peak times

stimwindow=30; % 5 s baseline, 15 s imaging

% timing for line scans (because line time in table is not correct)
FrameTime=592/11.8371;
LineScanLength=75184;
lineTime=1/(LineScanLength/FrameTime);

nLines=round(stimwindow/lineTime);
TimeX(1:nLines) = (1:nLines)*lineTime;

% Calculate the first peak onset time and AUC after stim
BL_time=118/11.8371;  % number of s for baseline
nBL_lines=BL_time/lineTime; % number of lines in baseline time
stimStart=TimeX(1,round(nBL_lines));  % exact time of stimulus start
baselineCorrectedTime=TimeX-stimStart;

%% load trace data

%AstrocyteExcelFile='E:\Data\Two_Photon_Data\GCaMP_RCaMP\Manuscript\Figures\DataForTraceOverlay_Heatmaps\cytovslck_Nostimvsstim_comparison\cyto-AstrocyteTraces_Stim.xlsx';
%NeuronalExcelFile ='E:\Data\Two_Photon_Data\GCaMP_RCaMP\Manuscript\Figures\DataForTraceOverlay_Heatmaps\cytovslck_Nostimvsstim_comparison\cyto-NeuronalTraces_Stim.xlsx';

%AstrocyteExcelFile='D:\Data\GCaMP_RCaMP\Manuscript\Neuron_2017\Revision\Figures\Lck-Linescan_AstrocyteTraces_Stim.xlsx';
%NeuronalExcelFile ='D:\Data\GCaMP_RCaMP\Manuscript\Neuron_2017\Revision\Figures\Lck-Linescan_NeuronalTraces_Stim.xlsx';

% for the first cohort

% Load trace data and run onset time calculation
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\\FilesforMatlab\Traces_1stCohort_Lck_nostim_vs_longstim_12_2017.mat');
load('D:\Data\GCaMP_RCaMP\Revision\Lck_GCaMP\FilesforMatlab\Linescans\EndfootLinescanTraces_AllMice_Lck_nostim_vs_longstim_01_2018.mat');

Shortstim=All_traces;

% trace info- ROIType etc.
% get info for plotting

% ONLY CONSIDER STIM DATA
% get rid of astrocyte neuropil traces
for xROI=1:size(Shortstim,1)
    NostimIdx(xROI)= strcmp(Shortstim{xROI,6},'Nostim');
end
Shortstim=Shortstim(~NostimIdx',:);


for iROI=1:length(Shortstim)
    % make new unique trial names
    Shortstim{iROI,16}=strcat(Shortstim{iROI,5},'_',Shortstim{iROI,4},'_',Shortstim{iROI,2});
    % unique ROI names
    Shortstim{iROI,17}=strcat(Shortstim{iROI,5},'_',Shortstim{iROI,4},'_',Shortstim{iROI,2}, '_',Shortstim{iROI,1});
end


%% separate out control data, IP3 mice and pharmacology

%all IP3 KO and WT
for xROI=1:size(Shortstim,1)
    IP3Match(xROI)= ~isempty(strfind(Shortstim{xROI,8},'IP'));
end
IP3=Shortstim(IP3Match',:);

% remove pharmacology from IP3 WT
for iROI=1:size(IP3,1)
    IP3PharmMatch(iROI)= ~isempty(strfind(IP3{iROI,7},'Control'));
end
IP3=IP3(IP3PharmMatch',:);


%Pharmacology data
for xROI=1:size(Shortstim,1)
    IPKO_Match(xROI)= ~isempty(strfind(Shortstim{xROI,8},'KO'));
end
Pharmacology=Shortstim(~IPKO_Match',:);

% control lines scans (for fast AC vs delayed AC vs Neuron comparisons)
for xROI=1:size(Pharmacology,1)
    ControlMatch(xROI)= ~isempty(strfind(Pharmacology{xROI,7},'Control'));
end
Lck_traces=Pharmacology(ControlMatch',:);


%% Mean Line scan Traces- Lck Control Data

for iROI=1:length(Lck_traces)
    rc_str(iROI)= ~isempty(strfind(Lck_traces{iROI,3},'RCaMP'));
    gc_str(iROI)= ~isempty(strfind(Lck_traces{iROI,3},'GCaMP'));
end

RCaMP=Lck_traces(rc_str',:);
GCaMP=Lck_traces(gc_str',:);


RCaMP_traces=[];
for iROI=1:length(RCaMP)
    respOTIdx1(iROI)=~isempty(find((RCaMP{iROI,15}>0 && RCaMP{iROI,15}<=Fast_OnsetWindow),1));
    tempY= RCaMP{iROI,11};
    tempY=tempY(1:nLines);
    if respOTIdx1(iROI)
        RCaMP_traces= horzcat(RCaMP_traces, tempY);
    end
end
RrespOT=RCaMP(respOTIdx1',:); % responding neurons
RCaMP_mean= nanmean(RCaMP_traces,2);

fastAC_traces=[];
slowAC_traces=[];
for iROI=1:length(GCaMP)
    fast_respOTIdx2(iROI)=~isempty(find((GCaMP{iROI,15}>0 && GCaMP{iROI,15}<=Fast_OnsetWindow),1));
    delayed_respOTIdx2(iROI)=~isempty(find((GCaMP{iROI,15}>Fast_OnsetWindow && GCaMP{iROI,15}<=AOnsetWindow),1));
    tempY= GCaMP{iROI,11};
    tempY=tempY(1:nLines);
    if fast_respOTIdx2(iROI)
        fastAC_traces= horzcat(fastAC_traces, tempY);
    end
    if delayed_respOTIdx2(iROI)
        slowAC_traces= horzcat(slowAC_traces, tempY);
    end
end
fastAC=GCaMP(fast_respOTIdx2',:); % fast responding astrocytes
delayedAC=GCaMP(delayed_respOTIdx2',:); % delayed responding astrocytes

GrespOT =vertcat(fastAC,delayedAC);
% mean trace of astrocytes
fastAC_mean= nanmean(fastAC_traces,2);
slowAC_mean= nanmean(slowAC_traces,2);


% %% low pass filter tests
%
% Hd2 = designfilt('lowpassfir','FilterOrder',100,'CutoffFrequency',20, ...
%        'DesignMethod','window','SampleRate',LineScanLength/FrameTime);
%
% y1 = filter(Hd2,fastAC_mean);
%
%    figure
%    hold on
%    plot(baselineCorrectedTime,fastAC_mean)
%    plot(baselineCorrectedTime,y1)
%
%    y2 = filter(Hd2,RCaMP_mean);
%
%    figure
%    hold on
%    plot(baselineCorrectedTime,RCaMP_mean)
%    plot(baselineCorrectedTime,y2)

%% shaded error bar with
green=[(27/255) (120/255) (55/255)];
purple=[(123/255) (50/255) (148/255)];
blue= [(0/255) (114/255) (178/255)];

%TimeX2(1:nframes2) = (1:nframes2)/FrameRate;

% SEM calculations
fastAC_SDTrace = nanstd(fastAC_traces');
fastAC_SEM=fastAC_SDTrace/sqrt(size(fastAC_traces,2));

slowAC_SDTrace = nanstd(slowAC_traces');
slowAC_SEM=slowAC_SDTrace/sqrt(size(slowAC_traces,2));

RC_SDTrace = nanstd(RCaMP_traces');
RC_SEM=RC_SDTrace/sqrt(size(RCaMP_traces,2));

figure('name', 'Lck means- plus SEM')
hold on
axis off
xlim([-3 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

%ylim([-0.2 3]);
lineProps.col = {blue};
mseb(baselineCorrectedTime,smooth(slowAC_mean',9),smooth(slowAC_SEM,9)',lineProps)
lineProps.col = {green};
mseb(baselineCorrectedTime,(smooth(fastAC_mean',9)+1.5),smooth(fastAC_SEM,9)',lineProps)
lineProps.col = {purple};
mseb(baselineCorrectedTime,(smooth(RCaMP_mean',9)+3),smooth(RC_SEM,9)',lineProps)
legend('Delayed EF','Fast EF','Neurons')
rectangle('Position', [0 -0.3 8 4])
plot([-2.1 -2.1],[0 1], 'k','LineWidth', 1)


