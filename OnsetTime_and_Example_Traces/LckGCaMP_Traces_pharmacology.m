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

% ONLY CONSIDER STIM DATA
% get rid of astrocyte neuropil traces
for xROI=1:size(Shortstim,1)
    NostimIdx(xROI)= strcmp(Shortstim{xROI,6},'Nostim');
end
Shortstim=Shortstim(~NostimIdx',:);


% Extract info for the required mice
for iROI=1:length(Shortstim)
    str1(iROI)= ~isempty(strfind(Shortstim{iROI,5},'WT_LR1'));
    str2(iROI)= ~isempty(strfind(Shortstim{iROI,5},'IPRG2'));
    str3(iROI)= ~isempty(strfind(Shortstim{iROI,5},'IPRG3'));
    str4(iROI)= ~isempty(strfind(Shortstim{iROI,5},'WT_LR4'));
    str5(iROI)= ~isempty(strfind(Shortstim{iROI,5},'ARG2'));
    str6(iROI)= ~isempty(strfind(Shortstim{iROI,5},'IPRG6'));
end

WT_LR1=Shortstim(str1',:);
IPRG2=Shortstim(str2',:);
IPRG3=Shortstim(str3',:);
WT_LR4=Shortstim(str4',:);
ARG2=Shortstim(str5',:);
IPRG6=Shortstim(str6',:);

Control=vertcat(WT_LR1,WT_LR4, IPRG2, IPRG3, IPRG6, ARG2);


for iROI=1:length(Control)
    %if a process ROI is overlapping with a soma or process, exclude it
    %from data
    nonOverlapIdx2(iROI)=~ischar(Control{iROI,12});
end

% remove overlapping processes
Control = Control(nonOverlapIdx2',:);


% Calculate the first peak onset time and AUC after stim
% peak onsets and AUC in the first second after stim for each ROI
for iROI= 1:length(Control)    
    trace=Control{iROI,9};
    FrameRate=Control{iROI,13};
    nframes2=length(trace);
    TimeX(1:nframes2) = (1:nframes2)/FrameRate;
    BL_time=Control{iROI, 8}/Control{iROI,13};
    baselineCorrectedTime=TimeX-BL_time;
    %first 1 sec after stim onset
    x1=Control{iROI, 8};
    x2=round((BL_time+1)*FrameRate);
    x3= round(FrameRate*(BL_time+10));
    
    % onset time
    Onsets=find_first_onset_time(baselineCorrectedTime(10:end), trace(10:end),2.5,2);
    if isempty(Onsets)
        Onsets=nan(1,1);
    end
    Control{iROI, 14}= Onsets;
    % AUC
    Control{iROI,15}=trapz(trace(x1:x2));
    Control{iROI,16}=trapz(trace(x1:x3));
    
end
%% load  pharmacology trace data

% Load trace data from drug treatments
load('D:\Data\GCaMP_RCaMP\Revision\Lck_GCaMP\FilesforMatlab\Pharmacology\Traces_pharmacology_Lck_nostim_vs_longstim_12_2017.mat');

Pharmacology=All_traces;

% trace info- ROIType etc.
% get info for plotting

% ONLY CONSIDER STIM DATA
% get rid of astrocyte neuropil traces
for xROI=1:size(Pharmacology,1)
    NostimIdx2(xROI)= strcmp(Pharmacology{xROI,6},'Nostim');
end
Pharmacology=Pharmacology(~NostimIdx2',:);


% Combine Control and Pharmacology
Control(:,17)={'Control'};
Pharmacology(:,18)=Pharmacology(:,7);
Pharmacology(:,7)=[];

Lck_traces=vertcat(Control, Pharmacology);
clearvars All_traces Control Pharmacology

for iROI=1:length(Lck_traces)
        %find ROITypes
    N_str= strfind(Lck_traces{iROI, 1},'N');
    D_str= strfind(Lck_traces{iROI, 1},'D'); %hand selected dendrite
    r_str= strfind(Lck_traces{iROI, 1},'r'); %FLIKA ROIs
    EF_str= strfind(Lck_traces{iROI, 1},'E');
    
    if ~isempty(N_str)
        Lck_traces{iROI,20}='Neuron';
    elseif ~isempty(D_str)
        Lck_traces{iROI,20}='Dendrite';
    elseif ~isempty(r_str)
        P_str= strfind(Lck_traces{iROI, 3},'GCaMP');
        if ~isempty(P_str)
            Lck_traces{iROI,20}='Process';
        else
            Lck_traces{iROI,20}='Dendrite';
        end
    elseif ~isempty(EF_str)
        Lck_traces{iROI,20}='Endfeet';
    end
    % make new unique trial names
    Lck_traces{iROI,18}=strcat(Lck_traces{iROI,5},'_',Lck_traces{iROI,4},'_',Lck_traces{iROI,2});
    % unique ROI names
    Lck_traces{iROI,19}=strcat(Lck_traces{iROI,5},'_',Lck_traces{iROI,4},'_',Lck_traces{iROI,2}, '_',Lck_traces{iROI,1});
end

%% Split data for each drug

% make a table for each drug
for iROI=1:length(Lck_traces)
    con_str(iROI)= ~isempty(strfind(Lck_traces{iROI,17},'Control'));
    DSP_str(iROI)= ~isempty(strfind(Lck_traces{iROI,17},'DSP'));
    atr_str(iROI)= ~isempty(strfind(Lck_traces{iROI,17},'Atropine'));
    Met_str(iROI)= ~isempty(strfind(Lck_traces{iROI,17},'Metergoline'));
    Traz_str(iROI)= ~isempty(strfind(Lck_traces{iROI,17},'Trazodone'));
    Praz_str(iROI)= ~isempty(strfind(Lck_traces{iROI,17},'Prazosin'));
end

Control=Lck_traces(con_str',:);
DSP4=Lck_traces(DSP_str',:);
Atropine=Lck_traces(atr_str',:);
Metergoline=Lck_traces(Met_str',:);
Trazodone=Lck_traces(Traz_str',:);
Prazosin=Lck_traces(Praz_str',:);

%% Control

for iROI=1:length(Control)
    rc_str(iROI)= ~isempty(strfind(Control{iROI,3},'RCaMP'));
    gc_str(iROI)= ~isempty(strfind(Control{iROI,3},'GCaMP'));
end

RCaMP=Control(rc_str',:);
GCaMP=Control(gc_str',:);

RCaMPCon_traces=[];
for iROI=1:length(RCaMP)
    respOTIdx1(iROI)=~isempty(find((RCaMP{iROI,14}>0 && RCaMP{iROI,14}<=Fast_OnsetWindow),1));
    tempY= RCaMP{iROI,9};
    if length(tempY)>nframes
        tempY=tempY(1:nframes);
        if respOTIdx1(iROI)
            RCaMPCon_traces= horzcat(RCaMPCon_traces, tempY);
        end
    end
end
RCaMPCon_mean= nanmean(RCaMPCon_traces,2);

fastCon_traces=[];
slowCon_traces=[];
for iROI=1:length(GCaMP)
    fastIdx(iROI)=~isempty(find((GCaMP{iROI,14}>0 && GCaMP{iROI,14}<=Fast_OnsetWindow),1));
    slowIdx(iROI)=~isempty(find((GCaMP{iROI,14}>Fast_OnsetWindow && GCaMP{iROI,14}<=AOnsetWindow),1));
    tempY= GCaMP{iROI,9};
    if length(tempY)>nframes
        tempY=tempY(1:nframes);
        if fastIdx(iROI)
            fastCon_traces= horzcat(fastCon_traces, tempY);
        end
        if slowIdx(iROI)
            slowCon_traces= horzcat(slowCon_traces, tempY);
        end
    end
end
clearvars fastIdx slowIdx

% mean trace of KO
fastControl_mean= nanmean(fastCon_traces,2);
slowControl_mean= nanmean(slowCon_traces,2);
RCaMPControl_mean= nanmean(RCaMPCon_traces,2);

% SEM calculations
fastControl_SEM=nanstd(fastCon_traces')/sqrt(size(fastCon_traces,2));
slowControl_SEM=nanstd(slowCon_traces')/sqrt(size(slowCon_traces,2));

RCaMPControl_SEM=nanstd(RCaMPCon_traces')/sqrt(size(RCaMPCon_traces,2));


% shaded error bar plot
green=[(27/255) (120/255) (55/255)];
purple=[(123/255) (50/255) (148/255)];
blue= [(0/255) (114/255) (178/255)];


% Plot
figure('name', 'Control- fast, delayed, and neurons ')
hold on
axis off
xlim([-3 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

lineProps.col = {purple};
mseb(baselineCorrectedTime(1:nframes),(RCaMPControl_mean'+2),smooth(RCaMPControl_SEM,9)',lineProps)
lineProps.col = {green};
mseb(baselineCorrectedTime(1:nframes),slowControl_mean',smooth(slowControl_SEM,9)',lineProps)
lineProps.col = {blue};
mseb(baselineCorrectedTime(1:nframes),(fastControl_mean'+1),smooth(fastControl_SEM,9)',lineProps)
rectangle('Position', [0 -0.3 8 4])
plot([-2.1 -2.1],[0 1], 'k','LineWidth', 1)
legend('Neuron','Delayed','Fast')

clearvars RCaMP GCaMP rc_str gc_str


%% example traces RCaMP, WT_LR4 spot3

for iROI=1:length(Control)
    wt_str(iROI)= ~isempty(strfind(Control{iROI,18},'WT_LR4'));
end
WT_LR4=Control(wt_str',:);

for iROI=1:length(WT_LR4)
    spt3B_str(iROI)= ~isempty(strfind(WT_LR4{iROI,18},'spot3'));
end
WT_LR4_spot3=WT_LR4(spt3B_str',:);

for iROI=1:length(WT_LR4_spot3)
    rcamp_str(iROI)= ~isempty(strfind(WT_LR4_spot3{iROI,3},'RCaMP'));
end
WT_LR4_spot3_RC=WT_LR4_spot3(rcamp_str',:);

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
plot([0 0],[0 50], 'k','LineWidth', 1)

% PLOT THE ROI MASK

%% DSP4

for iROI=1:length(DSP4)
    rc_str(iROI)= ~isempty(strfind(DSP4{iROI,3},'RCaMP'));
    gc_str(iROI)= ~isempty(strfind(DSP4{iROI,3},'GCaMP'));
end

RCaMP=DSP4(rc_str',:);
GCaMP=DSP4(gc_str',:);

RCaMPDSP_traces=[];
for iROI=1:length(RCaMP)
    respOTIdx1(iROI)=~isempty(find((RCaMP{iROI,14}>0 && RCaMP{iROI,14}<=Fast_OnsetWindow),1));
    tempY= RCaMP{iROI,9};
    if length(tempY)>nframes
        tempY=tempY(1:nframes);
        if respOTIdx1(iROI)
            RCaMPDSP_traces= horzcat(RCaMPDSP_traces, tempY);
        end
    end
end
RCaMPDSP_mean= nanmean(RCaMPDSP_traces,2);

fastDSP_traces=[];
slowDSP_traces=[];
for iROI=1:length(GCaMP)
    fastIdx(iROI)=~isempty(find((GCaMP{iROI,14}>0 && GCaMP{iROI,14}<=Fast_OnsetWindow),1));
    slowIdx(iROI)=~isempty(find((GCaMP{iROI,14}>Fast_OnsetWindow && GCaMP{iROI,14}<=AOnsetWindow),1));
    tempY= GCaMP{iROI,9};
    if length(tempY)>nframes
        tempY=tempY(1:nframes);
        if fastIdx(iROI)
            fastDSP_traces= horzcat(fastDSP_traces, tempY);
        end
        if slowIdx(iROI)
            slowDSP_traces= horzcat(slowDSP_traces, tempY);
        end
    end
end
clearvars fastIdx slowIdx

% mean trace of KO
fastDSP4_mean= nanmean(fastDSP_traces,2);
slowDSP4_mean= nanmean(slowDSP_traces,2);
RCaMPDSP4_mean= nanmean(RCaMPDSP_traces,2);

% SEM calculations
fastDSP4_SEM=nanstd(fastDSP_traces')/sqrt(size(fastDSP_traces,2));
slowDSP4_SEM=nanstd(slowDSP_traces')/sqrt(size(slowDSP_traces,2));

RCaMPDSP4_SEM=nanstd(RCaMPDSP_traces')/sqrt(size(RCaMPDSP_traces,2));


% shaded error bar plot
green=[(27/255) (120/255) (55/255)];
purple=[(123/255) (50/255) (148/255)];
blue= [(0/255) (114/255) (178/255)];


% Plot
figure('name', 'DSP4- fast, delayed, and neurons ')
hold on
axis off
xlim([-3 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

lineProps.col = {purple};
mseb(baselineCorrectedTime(1:nframes),(RCaMPDSP4_mean'+2),smooth(RCaMPDSP4_SEM,9)',lineProps)
lineProps.col = {green};
mseb(baselineCorrectedTime(1:nframes),slowDSP4_mean',smooth(slowDSP4_SEM,9)',lineProps)
lineProps.col = {blue};
mseb(baselineCorrectedTime(1:nframes),(fastDSP4_mean'+1),smooth(fastDSP4_SEM,9)',lineProps)
rectangle('Position', [0 -0.3 8 4])
plot([-2.1 -2.1],[0 1], 'k','LineWidth', 1)
legend('Neuron','Delayed','Fast')

clearvars RCaMP GCaMP rc_str gc_str


%% Atropine
for iROI=1:length(Atropine)
    rc_str(iROI)= ~isempty(strfind(Atropine{iROI,3},'RCaMP'));
    gc_str(iROI)= ~isempty(strfind(Atropine{iROI,3},'GCaMP'));
end

RCaMP=Atropine(rc_str',:);
GCaMP=Atropine(gc_str',:);

RCaMPAtropine_traces=[];
for iROI=1:length(RCaMP)
    respOTIdx1(iROI)=~isempty(find((RCaMP{iROI,14}>0 && RCaMP{iROI,14}<=Fast_OnsetWindow),1));
    tempY= RCaMP{iROI,9};
    if length(tempY)>nframes
        tempY=tempY(1:nframes);
        if respOTIdx1(iROI)
            RCaMPAtropine_traces= horzcat(RCaMPAtropine_traces, tempY);
        end
    end
end
RCaMPAtropine_mean= nanmean(RCaMPAtropine_traces,2);

fastAtropine_traces=[];
slowAtropine_traces=[];
for iROI=1:length(GCaMP)
    fastIdx(iROI)=~isempty(find((GCaMP{iROI,14}>0 && GCaMP{iROI,14}<=Fast_OnsetWindow),1));
    slowIdx(iROI)=~isempty(find((GCaMP{iROI,14}>Fast_OnsetWindow && GCaMP{iROI,14}<=AOnsetWindow),1));
    tempY= GCaMP{iROI,9};
    if length(tempY)>nframes
        tempY=tempY(1:nframes);
        if fastIdx(iROI)
            fastAtropine_traces= horzcat(fastAtropine_traces, tempY);
        end
        if slowIdx(iROI)
            slowAtropine_traces= horzcat(slowAtropine_traces, tempY);
        end
    end
end
clearvars fastIdx slowIdx

% mean trace of KO
fastAtropine_mean= nanmean(fastAtropine_traces,2);
slowAtropine_mean= nanmean(slowAtropine_traces,2);
RCaMPAtropine_mean= nanmean(RCaMPAtropine_traces,2);

% SEM calculations
fastAtropine_SEM=nanstd(fastAtropine_traces')/sqrt(size(fastAtropine_traces,2));
slowAtropine_SEM=nanstd(slowAtropine_traces')/sqrt(size(slowAtropine_traces,2));

RCaMPAtropine_SEM=nanstd(RCaMPAtropine_traces')/sqrt(size(RCaMPAtropine_traces,2));


% shaded error bar plot
green=[(27/255) (120/255) (55/255)];
purple=[(123/255) (50/255) (148/255)];
blue= [(0/255) (114/255) (178/255)];


% Plot
figure('name', 'Atropine- fast, delayed, and neurons ')
hold on
axis off
xlim([-3 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

lineProps.col = {purple};
mseb(baselineCorrectedTime(1:nframes),(RCaMPAtropine_mean'+2),smooth(RCaMPAtropine_SEM,9)',lineProps)
lineProps.col = {green};
mseb(baselineCorrectedTime(1:nframes),slowAtropine_mean',smooth(slowAtropine_SEM,9)',lineProps)
lineProps.col = {blue};
mseb(baselineCorrectedTime(1:nframes),(fastAtropine_mean'+1),smooth(fastAtropine_SEM,9)',lineProps)
rectangle('Position', [0 -0.3 8 4])
plot([-2.1 -2.1],[0 1], 'k','LineWidth', 1)
legend('Neuron','Delayed','Fast')

clearvars RCaMP GCaMP rc_str gc_str


%% Metergoline
for iROI=1:length(Metergoline)
    rc_str(iROI)= ~isempty(strfind(Metergoline{iROI,3},'RCaMP'));
    gc_str(iROI)= ~isempty(strfind(Metergoline{iROI,3},'GCaMP'));
end

RCaMP=Metergoline(rc_str',:);
GCaMP=Metergoline(gc_str',:);

RCaMPMetergoline_traces=[];
for iROI=1:length(RCaMP)
    respOTIdx1(iROI)=~isempty(find((RCaMP{iROI,14}>0 && RCaMP{iROI,14}<=Fast_OnsetWindow),1));
    tempY= RCaMP{iROI,9};
    if length(tempY)>nframes
        tempY=tempY(1:nframes);
        if respOTIdx1(iROI)
            RCaMPMetergoline_traces= horzcat(RCaMPMetergoline_traces, tempY);
        end
    end
end
RCaMPMetergoline_mean= nanmean(RCaMPMetergoline_traces,2);

fastMetergoline_traces=[];
slowMetergoline_traces=[];
for iROI=1:length(GCaMP)
    fastIdx(iROI)=~isempty(find((GCaMP{iROI,14}>0 && GCaMP{iROI,14}<=Fast_OnsetWindow),1));
    slowIdx(iROI)=~isempty(find((GCaMP{iROI,14}>Fast_OnsetWindow && GCaMP{iROI,14}<=AOnsetWindow),1));
    tempY= GCaMP{iROI,9};
    if length(tempY)>nframes
        tempY=tempY(1:nframes);
        if fastIdx(iROI)
            fastMetergoline_traces= horzcat(fastMetergoline_traces, tempY);
        end
        if slowIdx(iROI)
            slowMetergoline_traces= horzcat(slowMetergoline_traces, tempY);
        end
    end
end
clearvars fastIdx slowIdx

% mean trace of KO
fastMetergoline_mean= nanmean(fastMetergoline_traces,2);
slowMetergoline_mean= nanmean(slowMetergoline_traces,2);
RCaMPMetergoline_mean= nanmean(RCaMPMetergoline_traces,2);

% SEM calculations
fastMetergoline_SEM=nanstd(fastMetergoline_traces')/sqrt(size(fastMetergoline_traces,2));
slowMetergoline_SEM=nanstd(slowMetergoline_traces')/sqrt(size(slowMetergoline_traces,2));

RCaMPMetergoline_SEM=nanstd(RCaMPMetergoline_traces')/sqrt(size(RCaMPMetergoline_traces,2));


% shaded error bar plot
green=[(27/255) (120/255) (55/255)];
purple=[(123/255) (50/255) (148/255)];
blue= [(0/255) (114/255) (178/255)];


% Plot
figure('name', 'Metergoline- fast, delayed, and neurons ')
hold on
axis off
xlim([-3 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

lineProps.col = {purple};
mseb(baselineCorrectedTime(1:nframes),(RCaMPMetergoline_mean'+2),smooth(RCaMPMetergoline_SEM,9)',lineProps)
lineProps.col = {green};
mseb(baselineCorrectedTime(1:nframes),slowMetergoline_mean',smooth(slowMetergoline_SEM,9)',lineProps)
lineProps.col = {blue};
mseb(baselineCorrectedTime(1:nframes),(fastMetergoline_mean'+1),smooth(fastMetergoline_SEM,9)',lineProps)
rectangle('Position', [0 -0.3 8 4])
plot([-2.1 -2.1],[0 1], 'k','LineWidth', 1)
legend('Neuron','Delayed','Fast')

clearvars RCaMP GCaMP rc_str gc_str


%% Trazodone
for iROI=1:length(Trazodone)
    rc_str(iROI)= ~isempty(strfind(Trazodone{iROI,3},'RCaMP'));
    gc_str(iROI)= ~isempty(strfind(Trazodone{iROI,3},'GCaMP'));
end

RCaMP=Trazodone(rc_str',:);
GCaMP=Trazodone(gc_str',:);

RCaMPTrazodone_traces=[];
for iROI=1:length(RCaMP)
    respOTIdx1(iROI)=~isempty(find((RCaMP{iROI,14}>0 && RCaMP{iROI,14}<=Fast_OnsetWindow),1));
    tempY= RCaMP{iROI,9};
    if length(tempY)>nframes
        tempY=tempY(1:nframes);
        if respOTIdx1(iROI)
            RCaMPTrazodone_traces= horzcat(RCaMPTrazodone_traces, tempY);
        end
    end
end
RCaMPTrazodone_mean= nanmean(RCaMPTrazodone_traces,2);

fastTrazodone_traces=[];
slowTrazodone_traces=[];
for iROI=1:length(GCaMP)
    fastIdx(iROI)=~isempty(find((GCaMP{iROI,14}>0 && GCaMP{iROI,14}<=Fast_OnsetWindow),1));
    slowIdx(iROI)=~isempty(find((GCaMP{iROI,14}>Fast_OnsetWindow && GCaMP{iROI,14}<=AOnsetWindow),1));
    tempY= GCaMP{iROI,9};
    if length(tempY)>nframes
        tempY=tempY(1:nframes);
        if fastIdx(iROI)
            fastTrazodone_traces= horzcat(fastTrazodone_traces, tempY);
        end
        if slowIdx(iROI)
            slowTrazodone_traces= horzcat(slowTrazodone_traces, tempY);
        end
    end
end
clearvars fastIdx slowIdx

% mean trace of KO
fastTrazodone_mean= nanmean(fastTrazodone_traces,2);
slowTrazodone_mean= nanmean(slowTrazodone_traces,2);
RCaMPTrazodone_mean= nanmean(RCaMPTrazodone_traces,2);

% SEM calculations
fastTrazodone_SEM=nanstd(fastTrazodone_traces')/sqrt(size(fastTrazodone_traces,2));
slowTrazodone_SEM=nanstd(slowTrazodone_traces')/sqrt(size(slowTrazodone_traces,2));

RCaMPTrazodone_SEM=nanstd(RCaMPTrazodone_traces')/sqrt(size(RCaMPTrazodone_traces,2));


% shaded error bar plot
green=[(27/255) (120/255) (55/255)];
purple=[(123/255) (50/255) (148/255)];
blue= [(0/255) (114/255) (178/255)];


% Plot
figure('name', 'Trazodone- fast, delayed, and neurons ')
hold on
axis off
xlim([-3 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

lineProps.col = {purple};
mseb(baselineCorrectedTime(1:nframes),(RCaMPTrazodone_mean'+2),smooth(RCaMPTrazodone_SEM,9)',lineProps)
lineProps.col = {green};
mseb(baselineCorrectedTime(1:nframes),slowTrazodone_mean',smooth(slowTrazodone_SEM,9)',lineProps)
lineProps.col = {blue};
mseb(baselineCorrectedTime(1:nframes),(fastTrazodone_mean'+1),smooth(fastTrazodone_SEM,9)',lineProps)
rectangle('Position', [0 -0.3 8 4])
plot([-2.1 -2.1],[0 1], 'k','LineWidth', 1)
legend('Neuron','Delayed','Fast')

clearvars RCaMP GCaMP rc_str gc_str


%% Prazosin
for iROI=1:length(Prazosin)
    rc_str(iROI)= ~isempty(strfind(Prazosin{iROI,3},'RCaMP'));
    gc_str(iROI)= ~isempty(strfind(Prazosin{iROI,3},'GCaMP'));
end

RCaMP=Prazosin(rc_str',:);
GCaMP=Prazosin(gc_str',:);

RCaMPPrazosin_traces=[];
for iROI=1:length(RCaMP)
    respOTIdx1(iROI)=~isempty(find((RCaMP{iROI,14}>0 && RCaMP{iROI,14}<=Fast_OnsetWindow),1));
    tempY= RCaMP{iROI,9};
    if length(tempY)>nframes
        tempY=tempY(1:nframes);
        if respOTIdx1(iROI)
            RCaMPPrazosin_traces= horzcat(RCaMPPrazosin_traces, tempY);
        end
    end
end
RCaMPPrazosin_mean= nanmean(RCaMPPrazosin_traces,2);

fastPrazosin_traces=[];
slowPrazosin_traces=[];
for iROI=1:length(GCaMP)
    fastIdx(iROI)=~isempty(find((GCaMP{iROI,14}>0 && GCaMP{iROI,14}<=Fast_OnsetWindow),1));
    slowIdx(iROI)=~isempty(find((GCaMP{iROI,14}>Fast_OnsetWindow && GCaMP{iROI,14}<=AOnsetWindow),1));
    tempY= GCaMP{iROI,9};
    if length(tempY)>nframes
        tempY=tempY(1:nframes);
        if fastIdx(iROI)
            fastPrazosin_traces= horzcat(fastPrazosin_traces, tempY);
        end
        if slowIdx(iROI)
            slowPrazosin_traces= horzcat(slowPrazosin_traces, tempY);
        end
    end
end
clearvars fastIdx slowIdx

% mean trace of KO
fastPrazosin_mean= nanmean(fastPrazosin_traces,2);
slowPrazosin_mean= nanmean(slowPrazosin_traces,2);
RCaMPPrazosin_mean= nanmean(RCaMPPrazosin_traces,2);

% SEM calculations
fastPrazosin_SEM=nanstd(fastPrazosin_traces')/sqrt(size(fastPrazosin_traces,2));
slowPrazosin_SEM=nanstd(slowPrazosin_traces')/sqrt(size(slowPrazosin_traces,2));

RCaMPPrazosin_SEM=nanstd(RCaMPPrazosin_traces')/sqrt(size(RCaMPPrazosin_traces,2));


% shaded error bar plot
green=[(27/255) (120/255) (55/255)];
purple=[(123/255) (50/255) (148/255)];
blue= [(0/255) (114/255) (178/255)];


% Plot
figure('name', 'Prazosin- fast, delayed, and neurons ')
hold on
axis off
xlim([-3 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

lineProps.col = {purple};
mseb(baselineCorrectedTime(1:nframes),(RCaMPPrazosin_mean'+2),smooth(RCaMPPrazosin_SEM,9)',lineProps)
lineProps.col = {green};
mseb(baselineCorrectedTime(1:nframes),slowPrazosin_mean',smooth(slowPrazosin_SEM,9)',lineProps)
lineProps.col = {blue};
mseb(baselineCorrectedTime(1:nframes),(fastPrazosin_mean'+1),smooth(fastPrazosin_SEM,9)',lineProps)
rectangle('Position', [0 -0.3 8 4])
plot([-2.1 -2.1],[0 1], 'k','LineWidth', 1)
legend('Neuron','Delayed','Fast')

clearvars RCaMP GCaMP rc_str gc_str


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
mseb(baselineCorrectedTime(1:nframes),(fastControl_mean'+10),fastControl_SEM,lineProps)
lineProps.col = {green};
mseb(baselineCorrectedTime(1:nframes),(fastPrazosin_mean'+8),fastPrazosin_SEM,lineProps)
lineProps.col = {blue};
mseb(baselineCorrectedTime(1:nframes),(fastDSP4_mean'+6),fastDSP4_SEM,lineProps)
lineProps.col = {red};
mseb(baselineCorrectedTime(1:nframes),(fastTrazodone_mean'+4),fastTrazodone_SEM,lineProps)
lineProps.col = {yellow};
mseb(baselineCorrectedTime(1:nframes),(fastAtropine_mean'+2),fastAtropine_SEM,lineProps)
lineProps.col = {purple};
mseb(baselineCorrectedTime(1:nframes),fastMetergoline_mean',fastMetergoline_SEM,lineProps)

rectangle('Position', [0 -0.3 8 12])
plot([-2.1 -2.1],[0 1], 'k','LineWidth', 1)
legend('Control','Prazosin','DSP4','Trazodone','Atropine','Metergoline')


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
mseb(baselineCorrectedTime(1:nframes),(slowControl_mean'+6),slowControl_SEM,lineProps)
lineProps.col = {green};
mseb(baselineCorrectedTime(1:nframes),(slowPrazosin_mean'+5),slowPrazosin_SEM,lineProps)
lineProps.col = {blue};
mseb(baselineCorrectedTime(1:nframes),(slowDSP4_mean'+4),slowDSP4_SEM,lineProps)
lineProps.col = {red};
mseb(baselineCorrectedTime(1:nframes),(slowTrazodone_mean'+2.5),slowTrazodone_SEM,lineProps)
lineProps.col = {yellow};
mseb(baselineCorrectedTime(1:nframes),(slowAtropine_mean'+1),slowAtropine_SEM,lineProps)
lineProps.col = {purple};
mseb(baselineCorrectedTime(1:nframes),slowMetergoline_mean',slowMetergoline_SEM,lineProps)

rectangle('Position', [0 -0.3 8 9])
plot([-2.1 -2.1],[0 1], 'k','LineWidth', 1)
legend('Control','Prazosin','DSP4','Trazodone','Atropine','Metergoline')