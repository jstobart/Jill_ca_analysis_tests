
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
load('D:\Data\GCaMP_RCaMP\Revision\Lck_GCaMP\FilesforMatlab\Linescans\LinescanTraces_AllMice_Lck_nostim_vs_longstim_01_2018.mat');

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
mseb(baselineCorrectedTime,(smooth(fastAC_mean',9)+0.75),smooth(fastAC_SEM,9)',lineProps)
lineProps.col = {purple};
mseb(baselineCorrectedTime,(smooth(RCaMP_mean',9)+2),smooth(RC_SEM,9)',lineProps)

% H2=shadedErrorBar(TimeX,(fastAC_mean+0.75),fastAC_SEM)%,green,1);
% H3=shadedErrorBar(TimeX,(OT_RCaMP_mean+2),RC_SEM)%,purple,1);
rectangle('Position', [0 -0.3 8 4])
plot([-2.1 -2.1],[0 1], 'k','LineWidth', 1)


%% Data for heat maps

% % sort astrocytes by onset time
% [~, OT_GCidx] = sort([GrespOT{:,15}], 'ascend');
% GrespOT_sort=GrespOT(OT_GCidx,:);
%
% % sort neuronss by onset time
% [~, OT_RCidx] = sort([RrespOT{:,15}], 'ascend');
% RrespOT_sort=RrespOT(OT_RCidx,:);
%
% % astrocyte table
% AC_traces=[];
% for xROI= 1:length(GrespOT_sort)
%     tempOnset = GrespOT_sort{xROI,15}; % individual onset time
%     tempY = GrespOT_sort{xROI,11}; % individual trace
%     tempY=smooth(tempY(1:nLines),9);
%     temp=vertcat(tempOnset,tempY);
%     AC_traces=vertcat(AC_traces,temp');
% end
%
% % neuronal table
% N_traces=[];
% for xROI= 1:length(RrespOT_sort)
%     tempOnset = RrespOT_sort{xROI,15}; % individual onset time
%     tempY = RrespOT_sort{xROI,11}; % individual trace
%     tempY=smooth(tempY(1:nLines),9);
%     temp=vertcat(tempOnset,tempY);
%     N_traces=vertcat(N_traces,temp');
% end
%
%
% xlswrite(AstrocyteExcelFile, AC_traces)
% xlswrite(NeuronalExcelFile, N_traces)


%% Mean Traces- IP3 knockouts vs littermates

% ROI with a response to stimulation

%nframes2=round(stimwindow*FrameRate);

for iROI=1:length(IP3)
    ko_str(iROI)= ~isempty(strfind(IP3{iROI,8},'KO'));
    wt_str(iROI)= ~isempty(strfind(IP3{iROI,8},'WT'));
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
    fast_KOIdx(iROI)=~isempty(find((IP3KO_GC{iROI,15}>0 && IP3KO_GC{iROI,15}<=Fast_OnsetWindow),1));
    slow_KOIdx(iROI)=~isempty(find((IP3KO_GC{iROI,15}>Fast_OnsetWindow && IP3KO_GC{iROI,15}<=AOnsetWindow),1));
    tempY= IP3KO_GC{iROI,11};
    tempY=tempY(1:nLines);
    if fast_KOIdx(iROI)
        fastKO_traces= horzcat(fastKO_traces, tempY);
    end
    if slow_KOIdx(iROI)
        slowKO_traces= horzcat(slowKO_traces, tempY);
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
    fast_WTIdx(iROI)=~isempty(find((IP3WT_GC{iROI,15}>0 && IP3WT_GC{iROI,15}<=Fast_OnsetWindow),1));
    slow_WTIdx(iROI)=~isempty(find((IP3WT_GC{iROI,15}>Fast_OnsetWindow && IP3WT_GC{iROI,15}<=AOnsetWindow),1));
    tempY= IP3WT_GC{iROI,11};
    tempY=tempY(1:nLines);
    if fast_WTIdx(iROI)
        fastWT_traces= horzcat(fastWT_traces, tempY);
    end
    if slow_WTIdx(iROI)
        slowWT_traces= horzcat(slowWT_traces, tempY);
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
red= [(255/255) (0/255) (0/255)];

%TimeX2(1:nframes2) = (1:nframes2)/FrameRate;

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
figure('name', 'Knockout vs wildtype, fast vs delayed')
hold on
axis off
xlim([-3 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

lineProps.col = {blue};
mseb(baselineCorrectedTime,(smooth(slowKO_mean',9)+1),smooth(slowKO_SEM,9)',lineProps)
mseb(baselineCorrectedTime,(smooth(fastKO_mean',9)+4.5),smooth(fastKO_SEM,9)',lineProps)
lineProps.col = {green};
mseb(baselineCorrectedTime,(smooth(fastWT_mean',9)+2.5),smooth(fastWT_SEM,9)',lineProps)
mseb(baselineCorrectedTime,smooth(slowWT_mean',9),smooth(slowWT_SEM,9)',lineProps)
rectangle('Position', [0 -0.3 8 4])
plot([-2.1 -2.1],[0 1], 'k','LineWidth', 1)

%% Pharmacology

% only consider Gcamp traces
for iROI=1:length(Pharmacology)
    gcamp_str(iROI)= ~isempty(strfind(Pharmacology{iROI,3},'GCaMP'));
end
Pharma_GC=Pharmacology(gcamp_str',:);


% make a table for each drug
for iROI=1:length(Pharma_GC)
    con_str(iROI)= ~isempty(strfind(Pharma_GC{iROI,7},'Control'));
    DSP_str(iROI)= ~isempty(strfind(Pharma_GC{iROI,7},'DSP'));
    atr_str(iROI)= ~isempty(strfind(Pharma_GC{iROI,7},'Atropine'));
    Met_str(iROI)= ~isempty(strfind(Pharma_GC{iROI,7},'Metergoline'));
    Traz_str(iROI)= ~isempty(strfind(Pharma_GC{iROI,7},'Trazodone'));
    Praz_str(iROI)= ~isempty(strfind(Pharma_GC{iROI,7},'Prazosin'));
end

Control=Pharma_GC(con_str',:);
DSP4=Pharma_GC(DSP_str',:);
Atropine=Pharma_GC(atr_str',:);
Metergoline=Pharma_GC(Met_str',:);
Trazodone=Pharma_GC(Traz_str',:);
Prazosin=Pharma_GC(Praz_str',:);


%% Control
fastCon_traces=[];
slowCon_traces=[];
for iROI=1:length(Control)
    fastIdx(iROI)=~isempty(find((Control{iROI,15}>0 && Control{iROI,15}<=Fast_OnsetWindow),1));
    slowIdx(iROI)=~isempty(find((Control{iROI,15}>Fast_OnsetWindow && Control{iROI,15}<=AOnsetWindow),1));
    tempY= Control{iROI,11};
    tempY=tempY(1:nLines);
    if fastIdx(iROI)
        fastCon_traces= horzcat(fastCon_traces, tempY);
    end
    if slowIdx(iROI)
        slowCon_traces= horzcat(slowCon_traces, tempY);
    end
end
fastControl=Control(fastIdx',:); % fast responding astrocytes
delayedControl=Control(slowIdx',:); % delayed responding astrocytes
clearvars fastIdx slowIdx
respControl=vertcat(fastControl,delayedControl);

AllControl_traces=horzcat(fastCon_traces,slowCon_traces);

% mean trace of KO
fastControl_mean= nanmean(fastCon_traces,2);
slowControl_mean= nanmean(slowCon_traces,2);
AllControl_mean= nanmean(AllControl_traces,2);

% SEM calculations
fastControl_SEM=nanstd(fastCon_traces')/sqrt(size(fastCon_traces,2));
slowControl_SEM=nanstd(slowCon_traces')/sqrt(size(slowCon_traces,2));

AllControl_SEM=nanstd(AllControl_traces')/sqrt(size(AllControl_traces,2));


% Plot
figure('name', 'Control, fast vs delayed')
hold on
axis off
xlim([-3 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

lineProps.col = {red};
mseb(baselineCorrectedTime,smooth(AllControl_mean',9),smooth(AllControl_SEM,9)',lineProps)
lineProps.col = {green};
mseb(baselineCorrectedTime,(smooth(slowControl_mean',9)+1),smooth(slowControl_SEM,9)',lineProps)
lineProps.col = {blue};
mseb(baselineCorrectedTime,(smooth(fastControl_mean',9)+2),smooth(fastControl_SEM,9)',lineProps)
rectangle('Position', [0 -0.3 8 4])
plot([-2.1 -2.1],[0 1], 'k','LineWidth', 1)
legend('All','Delayed','Fast')

%% DSP4
fastDSP_traces=[];
slowDSP_traces=[];
for iROI=1:length(DSP4)
    fastIdx(iROI)=~isempty(find((DSP4{iROI,15}>0 && DSP4{iROI,15}<=Fast_OnsetWindow),1));
    slowIdx(iROI)=~isempty(find((DSP4{iROI,15}>Fast_OnsetWindow && DSP4{iROI,15}<=AOnsetWindow),1));
    tempY= DSP4{iROI,11};
    tempY=tempY(1:nLines);
    if fastIdx(iROI)
        fastDSP_traces= horzcat(fastDSP_traces, tempY);
    end
    if slowIdx(iROI)
        slowDSP_traces= horzcat(slowDSP_traces, tempY);
    end
end
fastDSP=DSP4(fastIdx',:); % fast responding astrocytes
delayedDSP=DSP4(slowIdx',:); % delayed responding astrocytes
clearvars fastIdx slowIdx
respDSP4 =vertcat(fastDSP,delayedDSP);

AllDSP_traces=horzcat(fastDSP_traces,slowDSP_traces);

% mean trace of KO
fastDSP_mean= nanmean(fastDSP_traces,2);
slowDSP_mean= nanmean(slowDSP_traces,2);
AllDSP_mean= nanmean(AllDSP_traces,2);

% SEM calculations
%fastDSP_SEM=nanstd(fastDSP_traces,2)/sqrt(size(fastDSP_traces,2));
slowDSP_SEM=nanstd(slowDSP_traces')/sqrt(size(slowDSP_traces,2));

AllDSP_SEM=nanstd(AllDSP_traces')/sqrt(size(AllDSP_traces,2));


% Plot
figure('name', 'DSP4, fast vs delayed')
hold on
axis off
xlim([-3 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

lineProps.col = {red};
mseb(baselineCorrectedTime,smooth(AllDSP_mean',9),smooth(AllDSP_SEM,9)',lineProps)
lineProps.col = {green};
mseb(baselineCorrectedTime,(smooth(slowDSP_mean',9)+1),smooth(slowDSP_SEM,9)',lineProps)
lineProps.col = {blue};
plot(baselineCorrectedTime,(smooth(fastDSP_mean',9)+2))%,lineProps)
rectangle('Position', [0 -0.3 8 4])
plot([-2.1 -2.1],[0 1], 'k','LineWidth', 1)
legend('All','Delayed','Fast')

%% Atropine
fastAtropine_traces=[];
slowAtropine_traces=[];
for iROI=1:length(Atropine)
    fastIdx(iROI)=~isempty(find((Atropine{iROI,15}>0 && Atropine{iROI,15}<=Fast_OnsetWindow),1));
    slowIdx(iROI)=~isempty(find((Atropine{iROI,15}>Fast_OnsetWindow && Atropine{iROI,15}<=AOnsetWindow),1));
    tempY= Atropine{iROI,11};
    tempY=tempY(1:nLines);
    if fastIdx(iROI)
        fastAtropine_traces= horzcat(fastAtropine_traces, tempY);
    end
    if slowIdx(iROI)
        slowAtropine_traces= horzcat(slowAtropine_traces, tempY);
    end
end
fastAtropine=Atropine(fastIdx',:); % fast responding astrocytes
delayedAtropine=Atropine(slowIdx',:); % delayed responding astrocytes
clearvars fastIdx slowIdx
respAtropine=vertcat(fastAtropine,delayedAtropine);

AllAtropine_traces=horzcat(fastAtropine_traces,slowAtropine_traces);

% mean trace of KO
fastAtropine_mean= nanmean(fastAtropine_traces,2);
slowAtropine_mean= nanmean(slowAtropine_traces,2);
AllAtropine_mean= nanmean(AllAtropine_traces,2);

% SEM calculations
fastAtropine_SEM=nanstd(fastAtropine_traces')/sqrt(size(fastAtropine_traces,2));
slowAtropine_SEM=nanstd(slowAtropine_traces')/sqrt(size(slowAtropine_traces,2));

AllAtropine_SEM=nanstd(AllAtropine_traces')/sqrt(size(AllAtropine_traces,2));


% Plot
figure('name', 'Atropine, fast vs delayed')
hold on
axis off
xlim([-3 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

lineProps.col = {red};
mseb(baselineCorrectedTime,smooth(AllAtropine_mean',9),smooth(AllAtropine_SEM,9)',lineProps)
lineProps.col = {green};
mseb(baselineCorrectedTime,(smooth(slowAtropine_mean',9)+1),smooth(slowAtropine_SEM,9)',lineProps)
lineProps.col = {blue};
mseb(baselineCorrectedTime,(smooth(fastAtropine_mean',9)+2),smooth(fastAtropine_SEM,9)',lineProps)
rectangle('Position', [0 -0.3 8 4])
plot([-2.1 -2.1],[0 1], 'k','LineWidth', 1)
legend('All','Delayed','Fast')

%% Metergoline
fastMetergoline_traces=[];
slowMetergoline_traces=[];
for iROI=1:length(Metergoline)
    fastIdx(iROI)=~isempty(find((Metergoline{iROI,15}>0 && Metergoline{iROI,15}<=Fast_OnsetWindow),1));
    slowIdx(iROI)=~isempty(find((Metergoline{iROI,15}>Fast_OnsetWindow && Metergoline{iROI,15}<=AOnsetWindow),1));
    tempY= Metergoline{iROI,11};
    tempY=tempY(1:nLines);
    if fastIdx(iROI)
        fastMetergoline_traces= horzcat(fastMetergoline_traces, tempY);
    end
    if slowIdx(iROI)
        slowMetergoline_traces= horzcat(slowMetergoline_traces, tempY);
    end
end
fastMetergoline=Metergoline(fastIdx',:); % fast responding astrocytes
delayedMetergoline=Metergoline(slowIdx',:); % delayed responding astrocytes
clearvars fastIdx slowIdx
respMetergoline=vertcat(fastMetergoline,delayedMetergoline);

AllMetergoline_traces=horzcat(fastMetergoline_traces,slowMetergoline_traces);

% mean traces
fastMetergoline_mean= nanmean(fastMetergoline_traces,2);
slowMetergoline_mean= nanmean(slowMetergoline_traces,2);
AllMetergoline_mean= nanmean(AllMetergoline_traces,2);

% SEM calculations
fastMetergoline_SEM=nanstd(fastMetergoline_traces')/sqrt(size(fastMetergoline_traces,2));
slowMetergoline_SEM=nanstd(slowMetergoline_traces')/sqrt(size(slowMetergoline_traces,2));

AllMetergoline_SEM=nanstd(AllMetergoline_traces')/sqrt(size(AllMetergoline_traces,2));


% Plot
figure('name', 'Metergoline, fast vs delayed')
hold on
axis off
xlim([-3 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

lineProps.col = {red};
mseb(baselineCorrectedTime,smooth(AllMetergoline_mean',9),smooth(AllMetergoline_SEM,9)',lineProps)
lineProps.col = {green};
mseb(baselineCorrectedTime,(smooth(slowMetergoline_mean',9)+1),smooth(slowMetergoline_SEM,9)',lineProps)
lineProps.col = {blue};
mseb(baselineCorrectedTime,(smooth(fastMetergoline_mean',9)+2),smooth(fastMetergoline_SEM,9)',lineProps)
rectangle('Position', [0 -0.3 8 4])
plot([-2.1 -2.1],[0 1], 'k','LineWidth', 1)
legend('All','Delayed','Fast')

%% Trazodone
fastTrazodone_traces=[];
slowTrazodone_traces=[];
for iROI=1:length(Trazodone)
    fastIdx(iROI)=~isempty(find((Trazodone{iROI,15}>0 && Trazodone{iROI,15}<=Fast_OnsetWindow),1));
    slowIdx(iROI)=~isempty(find((Trazodone{iROI,15}>Fast_OnsetWindow && Trazodone{iROI,15}<=AOnsetWindow),1));
    tempY= Trazodone{iROI,11};
    tempY=tempY(1:nLines);
    if fastIdx(iROI)
        fastTrazodone_traces= horzcat(fastTrazodone_traces, tempY);
    end
    if slowIdx(iROI)
        slowTrazodone_traces= horzcat(slowTrazodone_traces, tempY);
    end
end
fastTrazodone=Trazodone(fastIdx',:); % fast responding astrocytes
delayedTrazodone=Trazodone(slowIdx',:); % delayed responding astrocytes
clearvars fastIdx slowIdx
respTrazodone=vertcat(fastTrazodone,delayedTrazodone);

AllTrazodone_traces=horzcat(fastTrazodone_traces,slowTrazodone_traces);

% mean traces
fastTrazodone_mean= nanmean(fastTrazodone_traces,2);
slowTrazodone_mean= nanmean(slowTrazodone_traces,2);
AllTrazodone_mean= nanmean(AllTrazodone_traces,2);

% SEM calculations
fastTrazodone_SEM=nanstd(fastTrazodone_traces')/sqrt(size(fastTrazodone_traces,2));
slowTrazodone_SEM=nanstd(slowTrazodone_traces')/sqrt(size(slowTrazodone_traces,2));

AllTrazodone_SEM=nanstd(AllTrazodone_traces')/sqrt(size(AllTrazodone_traces,2));


% Plot
figure('name', 'Trazodone, fast vs delayed')
hold on
axis off
xlim([-3 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

lineProps.col = {red};
mseb(baselineCorrectedTime,smooth(AllTrazodone_mean',9),smooth(AllTrazodone_SEM,9)',lineProps)
lineProps.col = {green};
mseb(baselineCorrectedTime,(smooth(slowTrazodone_mean',9)+1),smooth(slowTrazodone_SEM,9)',lineProps)
lineProps.col = {blue};
mseb(baselineCorrectedTime,(smooth(fastTrazodone_mean',9)+2),smooth(fastTrazodone_SEM,9)',lineProps)
rectangle('Position', [0 -0.3 8 4])
plot([-2.1 -2.1],[0 1], 'k','LineWidth', 1)
legend('All','Delayed','Fast')
%% Prazosin
fastPrazosin_traces=[];
slowPrazosin_traces=[];
for iROI=1:length(Prazosin)
    fastIdx(iROI)=~isempty(find((Prazosin{iROI,15}>0 && Prazosin{iROI,15}<=Fast_OnsetWindow),1));
    slowIdx(iROI)=~isempty(find((Prazosin{iROI,15}>Fast_OnsetWindow && Prazosin{iROI,15}<=AOnsetWindow),1));
    tempY= Prazosin{iROI,11};
    tempY=tempY(1:nLines);
    if fastIdx(iROI)
        fastPrazosin_traces= horzcat(fastPrazosin_traces, tempY);
    end
    if slowIdx(iROI)
        slowPrazosin_traces= horzcat(slowPrazosin_traces, tempY);
    end
end
fastPrazosin=Prazosin(fastIdx',:); % fast responding astrocytes
delayedPrazosin=Prazosin(slowIdx',:); % delayed responding astrocytes
clearvars fastIdx slowIdx
respPrazosin=vertcat(fastPrazosin,delayedPrazosin);

AllPrazosin_traces=horzcat(fastPrazosin_traces,slowPrazosin_traces);

% mean trace of KO
fastPrazosin_mean= nanmean(fastPrazosin_traces,2);
slowPrazosin_mean= nanmean(slowPrazosin_traces,2);
AllPrazosin_mean= nanmean(AllPrazosin_traces,2);

% SEM calculations
fastPrazosin_SEM=nanstd(fastPrazosin_traces')/sqrt(size(fastPrazosin_traces,2));
slowPrazosin_SEM=nanstd(slowPrazosin_traces')/sqrt(size(slowPrazosin_traces,2));

AllPrazosin_SEM=nanstd(AllPrazosin_traces')/sqrt(size(AllPrazosin_traces,2));


% Plot
figure('name', 'Prazosin, fast vs delayed')
hold on
axis off
xlim([-3 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

lineProps.col = {red};
mseb(baselineCorrectedTime,smooth(AllPrazosin_mean',9),smooth(AllPrazosin_SEM,9)',lineProps)
lineProps.col = {green};
mseb(baselineCorrectedTime,(smooth(slowPrazosin_mean',9)+1),smooth(slowPrazosin_SEM,9)',lineProps)
lineProps.col = {blue};
mseb(baselineCorrectedTime,(smooth(fastPrazosin_mean',9)+2),smooth(fastPrazosin_SEM,9)',lineProps)
rectangle('Position', [0 -0.3 8 4])
plot([-2.1 -2.1],[0 1], 'k','LineWidth', 1)
legend('All','Delayed','Fast')

%% All Drugs FAST
green=[(27/255) (120/255) (55/255)];
purple=[(123/255) (50/255) (148/255)];
blue= [(0/255) (114/255) (178/255)];
red = [(255/255) (0/255) (0/255)];
black = [0 0 0];
yellow = [(210/255) (219/255) (31/255)];

% Plot
figure('name', 'All drugs, fast')
hold on
axis off
xlim([-3 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

lineProps.col = {black};
mseb(baselineCorrectedTime,(smooth(fastControl_mean',9)+10),smooth(fastControl_SEM,9)',lineProps)
lineProps.col = {green};
mseb(baselineCorrectedTime,(smooth(fastPrazosin_mean',9)+8),smooth(fastPrazosin_SEM,9)',lineProps)
%lineProps.col = {blue};
%plot(baselineCorrectedTime,(smooth(fastDSP_mean',9)+6))
lineProps.col = {red};
mseb(baselineCorrectedTime,(smooth(fastTrazodone_mean',9)+4),smooth(fastTrazodone_SEM,9)',lineProps)
lineProps.col = {yellow};
mseb(baselineCorrectedTime,(smooth(fastAtropine_mean',9)+2),smooth(fastAtropine_SEM,9)',lineProps)
lineProps.col = {purple};
mseb(baselineCorrectedTime,smooth(fastMetergoline_mean',9),smooth(fastMetergoline_SEM,9)',lineProps)

rectangle('Position', [0 -0.3 8 12])
plot([-2.1 -2.1],[0 1], 'k','LineWidth', 1)
legend('Control','Prazosin','Trazodone','Atropine','Metergoline')


%% All Drugs Delayed
green=[(27/255) (120/255) (55/255)];
purple=[(123/255) (50/255) (148/255)];
blue= [(0/255) (114/255) (178/255)];
red = [(255/255) (0/255) (0/255)];
black = [0 0 0];
yellow = [(210/255) (219/255) (31/255)];

% Plot
figure('name', 'All drugs, delayed')
hold on
axis off
xlim([-3 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

lineProps.col = {black};
mseb(baselineCorrectedTime,(smooth(slowControl_mean',9)+6),smooth(slowControl_SEM,9)',lineProps)
lineProps.col = {green};
mseb(baselineCorrectedTime,(smooth(slowPrazosin_mean',9)+5),smooth(slowPrazosin_SEM,9)',lineProps)
lineProps.col = {blue};
mseb(baselineCorrectedTime,(smooth(slowDSP_mean',9)+4),smooth(slowDSP_SEM,9)',lineProps)
lineProps.col = {red};
mseb(baselineCorrectedTime,(smooth(slowTrazodone_mean',9)+2.5),smooth(slowTrazodone_SEM,9)',lineProps)
lineProps.col = {yellow};
mseb(baselineCorrectedTime,(smooth(slowAtropine_mean',9)+1),smooth(slowAtropine_SEM,9)',lineProps)
lineProps.col = {purple};
mseb(baselineCorrectedTime,smooth(slowMetergoline_mean',9),smooth(slowMetergoline_SEM,9)',lineProps)

rectangle('Position', [0 -0.3 8 9])
plot([-2.1 -2.1],[0 1], 'k','LineWidth', 1)
legend('Control','Prazosin','DSP4','Trazodone','Atropine','Metergoline')