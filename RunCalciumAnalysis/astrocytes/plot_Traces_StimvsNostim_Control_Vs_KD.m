
clearvars
close all

% time windows based on stimulation
NOnsetWindow= 2; %1 for short stim % neuronal onset times
AOnsetWindow= 12; % 5 for short stim % astrocyte onset times
Fast_OnsetWindow=1.09;

NPTWindow= 9; % one second longer than stimulation for peak times
APTWindow= 15; % astrocyte longer than stimulation for peak times

stimwindow=20; % 5 s baseline, 15 s imaging


%% load trace data

load('D:\Data\GCaMP_RCaMP\NR1_KD\Results\FilesforMatlab\clean_data\74_traces_longtrials_clean.mat');
mouse74=All_traces;

load('D:\Data\GCaMP_RCaMP\NR1_KD\Results\FilesforMatlab\clean_data\92_traces_longtrials_clean.mat');
mouse92=All_traces;

load('D:\Data\GCaMP_RCaMP\NR1_KD\Results\FilesforMatlab\clean_data\94_traces_longtrials_clean.mat');
mouse94=All_traces;

load('D:\Data\GCaMP_RCaMP\NR1_KD\Results\FilesforMatlab\clean_data\95_traces_longtrials_clean.mat');
mouse95=All_traces;

load('D:\Data\GCaMP_RCaMP\NR1_KD\Results\FilesforMatlab\clean_data\96_traces_longtrials_clean.mat');
mouse96=All_traces;

load('D:\Data\GCaMP_RCaMP\NR1_KD\Results\FilesforMatlab\clean_data\alice_traces_longtrials_clean.mat');
mouseAlice=All_traces;

load('D:\Data\GCaMP_RCaMP\NR1_KD\Results\FilesforMatlab\clean_data\crazy8_traces_longtrials_clean.mat');
mouseCrazy8=All_traces;

%NS=vertcat(mouse74, mouse92, mouseAlice, mouseCrazy8); % non-silencing
NS=vertcat(mouse74, mouse92, mouseAlice, mouseCrazy8,mouse94,mouse95,mouse96);
KD=vertcat(mouse94, mouse95, mouse96); % knock down


% ONLY CONSIDER STIM DATA
for xROI=1:size(NS,1)
    NostimIdx(xROI)= strcmp(NS{xROI,6},'nostim');
end
NS_stim=NS(~NostimIdx',:);

for xROI=1:size(KD,1)
    NostimIdx2(xROI)= strcmp(KD{xROI,6},'nostim');
end
KD_stim=KD(~NostimIdx2',:);


for iROI=1:length(NS_stim)
    
    % make new unique trial names
    NS_stim{iROI,15}=strcat(NS_stim{iROI,5},'_',NS_stim{iROI,4},'_',NS_stim{iROI,2});
    % unique ROI names
    NS_stim{iROI,16}=strcat(NS_stim{iROI,5},'_',NS_stim{iROI,4},'_',NS_stim{iROI,2}, '_',NS_stim{iROI,1});
end

for iROI=1:length(KD_stim)
    
    % make new unique trial names
    KD_stim{iROI,15}=strcat(KD_stim{iROI,5},'_',KD_stim{iROI,4},'_',KD_stim{iROI,2});
    % unique ROI names
    KD_stim{iROI,16}=strcat(KD_stim{iROI,5},'_',KD_stim{iROI,4},'_',KD_stim{iROI,2}, '_',KD_stim{iROI,1});
end

for iROI=1:length(NS_stim)
    %if a process ROI is overlapping with a soma or process, exclude it
    %from data
    nonOverlapIdx(iROI)=~ischar(NS_stim{iROI,14});
end

% remove overlapping processes
NS_stim = NS_stim(nonOverlapIdx',:);

for iROI=1:length(KD_stim)
    %if a process ROI is overlapping with a soma or process, exclude it
    %from data
    nonOverlapIdx2(iROI)=~ischar(KD_stim{iROI,14});
end

% remove overlapping processes
KD_stim = KD_stim(nonOverlapIdx2',:);


% %%
% % Calculate the first peak onset time after stim
% 
% % % peak onsets and AUC in the first second after stim for each ROI
% % for iROI= 1:length(NS_stim)
% %     
    trace=NS_stim{1,8};
    FrameRate=NS_stim{1,11};
    nframes=length(trace);
    TimeX(1:nframes) = (1:nframes)/FrameRate;
    BL_time=4.92;
    baselineCorrectedTime=TimeX-BL_time;
% %     
% %     
% %     % onset time
% %     Onsets=find_first_onset_time(baselineCorrectedTime(10:end), trace(10:end),2.5,2);
% %     if isempty(Onsets)
% %         Onsets=nan(1,1);
% %     end
% %     NS_stim{iROI, 17}= Onsets;
% %     
% % end
% 
% 
% % peak onsets and AUC in the first second after stim for each ROI
% for iROI= 1:length(KD_stim)
%     
%     trace=KD_stim{iROI,8};
%     FrameRate=KD_stim{iROI,11};
%     nframes=length(trace);
%     TimeX(1:nframes) = (1:nframes)/FrameRate;
%     BL_time=4.92;
%     baselineCorrectedTime=TimeX-BL_time;
%     
%     
%     % onset time
%     Onsets=find_first_onset_time(baselineCorrectedTime(10:end), trace(10:end),2.5,2);
%     if isempty(Onsets)
%         Onsets=nan(1,1);
%     end
%     KD_stim{iROI, 17}= Onsets;
%     
% end


%% Mean Traces- Non-silencing Data

% ROI with a response to stimulation from non-silencing

nframes2=round(stimwindow*FrameRate);

for iROI=1:length(NS_stim)
    rc_str(iROI)= ~isempty(strfind(NS_stim{iROI,3},'RCaMP'));
    gc_str(iROI)= ~isempty(strfind(NS_stim{iROI,3},'GCaMP'));
end

RCaMP_NS=NS_stim(rc_str',:);
GCaMP_NS=NS_stim(gc_str',:);

% pull out neuronal traces that start with stimulation
OT_RCaMP_traces=[];
for iROI=1:length(RCaMP_NS)
    respOTIdx1(iROI)=~isempty(find((RCaMP_NS{iROI,12}>0 && RCaMP_NS{iROI,12}<=Fast_OnsetWindow),1));
    tempY= RCaMP_NS{iROI,8};
    BL_time=4.92;
    
    if respOTIdx1(iROI)
        OT_RCaMP_traces= horzcat(OT_RCaMP_traces, tempY);
    end
end
RrespOT=RCaMP_NS(respOTIdx1',:); % responding neurons

OT_RCaMP_mean= nanmean(OT_RCaMP_traces,2); % mean of responding neuronal traces

% pull out responding astrocytes (fast and delayed)
fastAC_traces=[];
slowAC_traces=[];
for iROI=1:length(GCaMP_NS)
    fast_respOTIdx2(iROI)=~isempty(find((GCaMP_NS{iROI,12}>0 && GCaMP_NS{iROI,12}<=Fast_OnsetWindow),1));
    delayed_respOTIdx2(iROI)=~isempty(find((GCaMP_NS{iROI,12}>Fast_OnsetWindow && GCaMP_NS{iROI,12}<=AOnsetWindow),1));
    tempY= GCaMP_NS{iROI,8};
    BL_time=4.92;
    if fast_respOTIdx2(iROI)
        fastAC_traces= horzcat(fastAC_traces, tempY);
    end
    if delayed_respOTIdx2(iROI)
        slowAC_traces= horzcat(slowAC_traces, tempY);
    end
end

fastAC=GCaMP_NS(fast_respOTIdx2',:); % fast responding astrocytes
delayedAC=GCaMP_NS(delayed_respOTIdx2',:); % delayed responding astrocytes

GrespOT =vertcat(fastAC,delayedAC);
% mean trace of astrocytes
fastAC_mean= nanmean(fastAC_traces,2);
slowAC_mean= nanmean(slowAC_traces,2);

allAC_traces=horzcat(slowAC_traces,fastAC_traces);
allAC_mean=nanmean(allAC_traces,2);
allAC_SDTrace = nanstd(allAC_traces');
allAC_SEM=allAC_SDTrace/sqrt(size(allAC_traces,2));

%% shaded error bar with
green=[(27/255) (120/255) (55/255)];
purple=[(123/255) (50/255) (148/255)];
blue= [(0/255) (114/255) (178/255)];

% SEM calculations
fastAC_SDTrace = nanstd(fastAC_traces');
fastAC_SEM=fastAC_SDTrace/sqrt(size(fastAC_traces,2));

slowAC_SDTrace = nanstd(slowAC_traces');
slowAC_SEM=slowAC_SDTrace/sqrt(size(slowAC_traces,2));

RC_SDTrace = nanstd(OT_RCaMP_traces');
RC_SEM=RC_SDTrace/sqrt(size(OT_RCaMP_traces,2));

figure('name', 'non-silencing means- plus SEM')
hold on
axis off
xlim([-1 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

%ylim([-0.2 3]);
lineProps.col = {blue};
mseb(TimeX,slowAC_mean',slowAC_SEM,lineProps)
lineProps.col = {green};
mseb(TimeX,(fastAC_mean'+0.75),fastAC_SEM,lineProps)
lineProps.col = {purple};
mseb(TimeX,(OT_RCaMP_mean'+2),RC_SEM,lineProps)

% H2=shadedErrorBar(TimeX,(fastAC_mean+0.75),fastAC_SEM)%,green,1);
% H3=shadedErrorBar(TimeX,(OT_RCaMP_mean+2),RC_SEM)%,purple,1);
rectangle('Position', [5 -0.3 8 4])
plot([-0.1 -0.1],[0 1], 'k','LineWidth', 1)


figure('name', 'astrocyte ROIs means- plus SEM')
hold on
axis off
xlim([-1 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

%ylim([-0.2 3]);
lineProps.col = {green};
mseb(TimeX(5:end),smooth(allAC_mean(5:end),7),smooth(allAC_SEM(5:end),7)',lineProps)


rectangle('Position', [5 -0.1 8 0.8])
plot([-1 -1],[0 0.5], 'k','LineWidth', 2)


figure('name', 'neuron ROIs means- plus SEM')
hold on
axis off
xlim([-1 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

%ylim([-0.2 3]);
lineProps.col = {purple};
mseb(TimeX(5:end),(smooth(OT_RCaMP_mean(5:end),7)),smooth(RC_SEM(5:end),7)',lineProps)


rectangle('Position', [5 -0.1 8 2])
plot([-1 -1],[0 0.5], 'k','LineWidth', 2)

%% mean traces for knockdown
for iROI=1:length(KD_stim)
    rc_str2(iROI)= ~isempty(strfind(KD_stim{iROI,3},'RCaMP'));
    gc_str2(iROI)= ~isempty(strfind(KD_stim{iROI,3},'GCaMP'));
end

RCaMP_KD=KD_stim(rc_str2',:);
GCaMP_KD=KD_stim(gc_str2',:);

% pull out responding neuron traces based on onset time
OT_RCaMP_traces=[];
for iROI=1:length(RCaMP_KD)
    respOTIdx2(iROI)=~isempty(find((RCaMP_KD{iROI,12}>0 && RCaMP_KD{iROI,12}<=Fast_OnsetWindow),1));
    tempY= RCaMP_KD{iROI,8};
    BL_time=4.92;
    
    if respOTIdx2(iROI)
        OT_RCaMP_traces= horzcat(OT_RCaMP_traces, tempY);
    end
end
RrespOT2=RCaMP_KD(respOTIdx2',:); % responding neurons
OT_RCaMP_mean2= nanmean(OT_RCaMP_traces,2);

% pull out responding neurons
fastAC_traces=[];
slowAC_traces=[];
for iROI=1:length(GCaMP_KD)
    fast_respOTIdx3(iROI)=~isempty(find((GCaMP_KD{iROI,12}>0 && GCaMP_KD{iROI,12}<=Fast_OnsetWindow),1));
    delayed_respOTIdx3(iROI)=~isempty(find((GCaMP_KD{iROI,12}>Fast_OnsetWindow && GCaMP_KD{iROI,12}<=AOnsetWindow),1));
    tempY= GCaMP_KD{iROI,8};
    BL_time=4.92;
    if fast_respOTIdx3(iROI)
        fastAC_traces= horzcat(fastAC_traces, tempY);
    end
    if delayed_respOTIdx3(iROI)
        slowAC_traces= horzcat(slowAC_traces, tempY);
    end
end

fastAC2=GCaMP_KD(fast_respOTIdx3',:); % fast responding astrocytes
delayedAC2=GCaMP_KD(delayed_respOTIdx3',:); % delayed responding astrocytes

GrespOT2 =vertcat(fastAC2,delayedAC2);
% mean trace of astrocytes
fastAC_mean2= nanmean(fastAC_traces,2);
slowAC_mean2= nanmean(slowAC_traces,2);


% SEM calculations
fastAC_SDTrace2 = nanstd(fastAC_traces');
fastAC_SEM2=fastAC_SDTrace2/sqrt(size(fastAC_traces,2));

slowAC_SDTrace2 = nanstd(slowAC_traces');
slowAC_SEM2=slowAC_SDTrace2/sqrt(size(slowAC_traces,2));

RC_SDTrace2 = nanstd(OT_RCaMP_traces');
RC_SEM2=RC_SDTrace2/sqrt(size(OT_RCaMP_traces,2));

figure('name', 'knock down means- plus SEM')
hold on
axis off
xlim([-1 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

%ylim([-0.2 3]);
lineProps.col = {blue};
mseb(TimeX,slowAC_mean2',slowAC_SEM2,lineProps)
lineProps.col = {green};
mseb(TimeX,(fastAC_mean2'+0.75),fastAC_SEM2,lineProps)
lineProps.col = {purple};
mseb(TimeX,(OT_RCaMP_mean2'+2),RC_SEM2,lineProps)

% H2=shadedErrorBar(TimeX,(fastAC_mean+0.75),fastAC_SEM)%,green,1);
% H3=shadedErrorBar(TimeX,(OT_RCaMP_mean+2),RC_SEM)%,purple,1);
rectangle('Position', [5 -0.3 8 4])
plot([-0.1 -0.1],[0 1], 'k','LineWidth', 1)


%% Mean Traces- knockdowns vs controls
% 
% %3
% figure('name', 'Knockout fast vs wildtype fast')
% hold on
% axis off
% xlim([-1 25]);
% lineProps.width = 1;
% lineProps.edgestyle = ':';
% 
% %ylim([-0.2 3]);
% lineProps.col = {blue};
% mseb(TimeX2,fastKO_mean',fastKO_SEM,lineProps)
% lineProps.col = {green};
% mseb(TimeX2,(fastWT_mean'+1),fastWT_SEM,lineProps)
% 
% rectangle('Position', [5 -0.3 8 4])
% plot([-0.1 -0.1],[0 1], 'k','LineWidth', 1)

% %4
% figure('name', 'Knockout delayed vs wildtype delayed')
% hold on
% axis off
% xlim([-1 25]);
% lineProps.width = 1;
% lineProps.edgestyle = ':';
% 
% %ylim([-0.2 3]);
% lineProps.col = {blue};
% mseb(TimeX2,slowKO_mean',slowKO_SEM,lineProps)
% lineProps.col = {green};
% mseb(TimeX2,(slowWT_mean'+1),slowWT_SEM,lineProps)
% 
% rectangle('Position', [5 -0.3 8 4])
% plot([-0.1 -0.1],[0 1], 'k','LineWidth', 1)
% 
% 
% %5
% figure('name', 'Knockout vs wildtype, fast vs delayed')
% hold on
% axis off
% xlim([-1 25]);
% lineProps.width = 1;
% lineProps.edgestyle = ':';
% 
% lineProps.col = {blue};
% mseb(TimeX2,(fastKO_mean'+2),fastKO_SEM,lineProps)
% lineProps.col = {green};
% mseb(TimeX2,(fastWT_mean'+3),fastWT_SEM,lineProps)
% lineProps.col = {blue};
% mseb(TimeX2,slowKO_mean',slowKO_SEM,lineProps)
% lineProps.col = {green};
% mseb(TimeX2,(slowWT_mean'+1),slowWT_SEM,lineProps)
% 
% rectangle('Position', [5 -0.3 8 4])
% plot([-0.1 -0.1],[0 1], 'k','LineWidth', 1)






