
clearvars
close all

%% Load data

% time windows based on stimulation
NOnsetWindow= 2; %1 for short stim % neuronal onset times
AOnsetWindow= 12; % 5 for short stim % astrocyte onset times
Fast_AOnsetWindow=1.09;
Fast_NOnsetWindow=1.09;
NPTWindow= 9; % one second longer than stimulation for peak times
APTWindow= 15; % astrocyte longer than stimulation for peak times

stimwindow=20; % 5 s baseline, 15 s imaging

%AstrocyteExcelFile='E:\Data\Two_Photon_Data\GCaMP_RCaMP\Manuscript\Figures\DataForTraceOverlay_Heatmaps\cytovslck_Nostimvsstim_comparison\cyto-AstrocyteTraces_Stim.xlsx';
%NeuronalExcelFile ='E:\Data\Two_Photon_Data\GCaMP_RCaMP\Manuscript\Figures\DataForTraceOverlay_Heatmaps\cytovslck_Nostimvsstim_comparison\cyto-NeuronalTraces_Stim.xlsx';
%peak data
peaks1=load('J:\Jill_Stobart\P2Y1_Data\Results\lckGC_P2Y1_timepoints_peaks_05_2018.mat');
peaks2=load('J:\Jill_Stobart\P2Y1_Data\Results\RC_P2Y1_timepoints_peaks_05_2018.mat');
% Load trace data
traces1=load('J:\Jill_Stobart\P2Y1_Data\Results\lckGC_P2Y1_timepoints_traces_05_2018.mat');
traces2=load('J:\Jill_Stobart\P2Y1_Data\Results\RC_P2Y1_timepoints_traces_05_2018.mat');

saveFile='J:\Jill_Stobart\P2Y1_Data\Results\lckGC&RC_P2Y1_timepoints_onsetTimes_05_2018.csv';

LckPeaks=peaks1.AllData2(2:end,:);
RCPeaks=peaks2.AllData2(2:end,:);

ShortstimPeaks=vertcat(peaks1.AllData2(2:end,:), peaks2.AllData2(2:end,:));

% fix single digit trial names
for x=1:length(ShortstimPeaks)
    number_pos=regexp(ShortstimPeaks{x,12},'[0-9]');
    if length(number_pos)<2
        trialname=ShortstimPeaks{x,12};
        number=str2double(trialname(number_pos));
        ShortstimPeaks{x,12}=strcat('trial',num2str(number,'%02d'));
    end
end

Shortstim=vertcat(traces1.AllData2(2:end,:), traces2.AllData2(2:end,:));

for x=1:length(Shortstim)
    number_pos=regexp(Shortstim{x,2},'[0-9]');
    if length(number_pos)<2
        trialname=Shortstim{x,2};
        number=str2double(trialname(number_pos));
        Shortstim{x,2}=strcat('trial',num2str(number,'%02d'));
    end
end

%% Peak info- ROiTypes, etc.

for iROI=1:length(ShortstimPeaks)
    %find ROITypes
    N_str= strfind(ShortstimPeaks{iROI, 1},'N');
    D_str= strfind(ShortstimPeaks{iROI, 1},'D'); %hand selected dendrite
    r_str= strfind(ShortstimPeaks{iROI, 1},'r'); %FLIKA ROIs
    EF_str= strfind(ShortstimPeaks{iROI, 1},'E');
    
    if ~isempty(N_str)
        ShortstimPeaks{iROI,23}='Neuron';
    elseif ~isempty(D_str)
        ShortstimPeaks{iROI,23}='Dendrite';
    elseif ~isempty(r_str)
        P_str= strfind(ShortstimPeaks{iROI, 3},'LckGCaMP');
        if ~isempty(P_str)
            ShortstimPeaks{iROI,23}='Process';
        else
            ShortstimPeaks{iROI,23}='Dendrite';
        end
    elseif ~isempty(EF_str)
        ShortstimPeaks{iROI,23}='Endfeet';
    end
    
    % make new unique trial names
    ShortstimPeaks{iROI,24}=strcat(ShortstimPeaks{iROI,13},'_',ShortstimPeaks{iROI,15},'_',ShortstimPeaks{iROI,12});
    % unique ROI names
    ShortstimPeaks{iROI,25}=strcat(ShortstimPeaks{iROI,13},'_',ShortstimPeaks{iROI,15},'_',ShortstimPeaks{iROI,12}, '_',ShortstimPeaks{iROI,10});
end

for iROI=1:length(ShortstimPeaks)
    %if a process ROI is overlapping with a soma or process, exclude it
    %from data
    nonOverlapIdx(iROI)=~ischar(ShortstimPeaks{iROI,19});
end

% remove overlapping processes
ShortstimPeaks = ShortstimPeaks(nonOverlapIdx',:);


%% trace info- ROIType etc.
% get info for plotting


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
        P_str= strfind(Shortstim{iROI, 3},'LckGCaMP');
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



%% Calculate the first peak onset time and AUC after stim


% peak onsets and AUC in the first second after stim for each ROI
for iROI= 1:length(Shortstim)
    
    trace=Shortstim{iROI,8};
    FrameRate=11.84;
    nframes=length(trace);
    TimeX(1:nframes) = (1:nframes)/FrameRate;
    BL_time=round(5*FrameRate)/11.84;
    baselineCorrectedTime=TimeX-BL_time;
    %first 1 sec after stim onset
    x1=round(5*FrameRate);
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

% table for importing into R
if ~exist('saveFiles3','file')
    ShortstimSave=Shortstim;
    ShortstimSave(:,8)=[];
    ShortstimSave(:,8)=[];
    ShortstimSave(:,8)=[];
    ShortstimSave(:,8)=[];
    names2={'ROI','Trial','Channel','Spot','Animal', 'Condition','depth','PixelSize','TimePoint',...
        'ROIType','Spot_trial','ROIs_trial','OnsetTime','TraceAUC1','TraceAUC10'};
    ShortstimSave2=vertcat(names2, ShortstimSave);
    cell2csv(saveFile,ShortstimSave2);
end




%% Responding Neurons and Astrocytes based on onset times

% ROI with a response to stimulation
% must have onset time during stimulation?

for iROI=1:length(Shortstim)
    rc_str(iROI)= ~isempty(strfind(Shortstim{iROI,3},'RCaMP'));
    gc_str(iROI)= ~isempty(strfind(Shortstim{iROI,3},'LckGCaMP'));
end

RCaMP=Shortstim(rc_str',:);
GCaMP=Shortstim(gc_str',:);

OT_RCaMP_traces=[];
for iROI=1:length(RCaMP)
    respOTIdx1(iROI)=~isempty(find(NOnsetWindow>RCaMP{iROI,17}>0,1));
    tempY= RCaMP{iROI,8};
    if length(tempY)>590
        tempY=tempY(1:nframes); %(stimwindow*FrameRate));
        if respOTIdx1(iROI)
            OT_RCaMP_traces= horzcat(OT_RCaMP_traces, tempY);
        end
    end
end
RrespOT=RCaMP(respOTIdx1',:); % responding neurons
OT_RCaMP_mean= mean(OT_RCaMP_traces');

OT_GCaMP_traces=[];
for iROI=1:length(GCaMP)
    respOTIdx2(iROI)=~isempty(find(GCaMP{iROI,17}>0 && GCaMP{iROI,17}<AOnsetWindow));
    tempY= GCaMP{iROI,8};
    if length(tempY)>590
        tempY=tempY(1:nframes); %(stimwindow*FrameRate));
        if respOTIdx2(iROI)
            OT_GCaMP_traces= horzcat(OT_GCaMP_traces, tempY);
        end
    end
end
GrespOT=GCaMP(respOTIdx2',:); % responding astrocytes
% mean trace of astrocytes
OT_GCaMP_mean= mean(OT_GCaMP_traces');

%% RCaMP vs GCaMP plots of responding ROIs based on onset
figure ('name', 'Overlaid traces: All RCaMP responding OT ROIs')
hold on
axis off
for xROI= 1:size(OT_RCaMP_traces,2)
    tempY = OT_RCaMP_traces(:,xROI);
    if length(tempY)>590
        tempY=tempY(1:nframes);
        grey = [0.8,0.8,0.8];
        plot(TimeX,tempY,'Color',grey,'LineWidth',0.01);
    end
end
plot(TimeX, OT_RCaMP_mean, 'Color', 'k','LineWidth',1);

figure ('name', 'Overlaid traces: All GCaMP responding OT ROIs')
hold on
axis off
for xROI= 1:size(OT_GCaMP_traces,2)
    tempY = OT_GCaMP_traces(:,xROI);
    if length(tempY)>590
        tempY=tempY(1:nframes);
        grey = [0.8,0.8,0.8];
        plot(TimeX,tempY,'Color',grey,'LineWidth',0.01);
    end
end
plot(TimeX, OT_GCaMP_mean, 'Color', 'k','LineWidth',1);


%% Histogram of peak time differences between neurons and astrocytes
% compare peak times of responding neurons to responding astrocytes in a field of view
% shift traces of each responding ROI by neuronal peak times

timeDiffs=[];
ShiftedTraces=[];

uniqueGTrials=unique(GrespOT(:,14));

for iTrial=1:size(uniqueGTrials,1)
    
    %find matching trials
    for xTrial= 1:size(GrespOT,1)
        Trial_str(xTrial)= strcmp(GrespOT{xTrial, 14},uniqueGTrials{iTrial});
    end
    GtrialData=GrespOT(Trial_str,:);
    
    
    %find neuron and neuropil data trials
    for rTrial= 1:size(RrespOT,1)
        RTrial_str(rTrial)= strcmp(RrespOT{rTrial, 14},uniqueGTrials{iTrial});
    end
    RtrialData=RrespOT(RTrial_str,:);
    
    %compare neuronal and astrocyte onset times
    for iN= 1:size(RtrialData, 1)
        for iA= 1:size(GtrialData,1)
            timeComparisons{iA,1}=RtrialData{iN,14}; % unique trial name
            timeComparisons{iA,2}=RtrialData{iN,13}; % neuron ROI type
            timeComparisons{iA,3}=RtrialData{iN,15}; % neuron name
            timeComparisons{iA,4}=RtrialData{iN,16}; % neuronal onset time
            timeComparisons{iA,5}=RtrialData{iN,17}; % neuronal auc in 1st s
            timeComparisons{iA,6}=GtrialData{iA,13}; % astrocyte ROI type
            timeComparisons{iA,7}=GtrialData{iA,15}; % astrocyte ROI name
            timeComparisons{iA,8}=GtrialData{iA,16}; % astrocyte onset time
            timeComparisons{iA,9}=GtrialData{iA,17}; % astrocyte auc
            timeComparisons{iA,10}=GtrialData{iA,16}-RtrialData{iN,16}; % onset time (AC) minus neurons onset
            timeComparisons{iA,11}=RtrialData{iN,16}-GtrialData{iA,16}; % neurons time (AC) minus  (N)
        end
        timeDiffs=vertcat(timeDiffs,timeComparisons);
    end
    
    for iN= 1:size(RtrialData, 1)
        % shift traces of all ROIs by each neuronal peak time
        for iTrace=1:size(GtrialData,1)
            ShiftTrace{iTrace,1}=RtrialData{iN,14}; %unique trial name
            ShiftTrace{iTrace,2}=RtrialData{iN,13}; % Neuron Type (soma or neuropil)
            ShiftTrace{iTrace,3}=RtrialData{iN,15}; % Neuron name
            ShiftTrace{iTrace,4}=RtrialData{iN,16}; % Neuron onset time
            ShiftTrace{iTrace,5}=RtrialData{iN,17}; % Neuron auc
            ShiftTrace{iTrace,6}=RtrialData{iN,8}; % Neuron trace
            ShiftTrace{iTrace,7}=GtrialData{iTrace,13}; % ROI type (ROI 2)
            ShiftTrace{iTrace,8}=GtrialData{iTrace,15}; % ROI name (ROI 2)
            ShiftTrace{iTrace,9}=GtrialData{iTrace,8}; % ROI 2 trace
            ShiftTrace{iTrace,10}=TimeX-(5+RtrialData{iN,16}); % minus neuronal onset time
            ShiftTrace{iTrace,11}=TimeX-(5+GtrialData{iTrace,16}); % minus astrocyte onset time
        end
        ShiftedTraces=vertcat(ShiftedTraces,ShiftTrace);
        clear timeComparisons ShiftTrace
    end
    clear GtrialData RtrialData
    
end

for iROI=1:size(timeDiffs,1)
    A_OT_minus_N_OT(iROI)=timeDiffs{iROI,10};
    N_OT_minus_A_OT(iROI)=timeDiffs{iROI,11};
end

% histograms
figure('name', 'histogram: A_OT_minus_N_OT')
histogram(A_OT_minus_N_OT)

figure('name', 'histogram: N_OT_minus_A_OT')
histogram(N_OT_minus_A_OT)






%% responding ROIs based on peak times (max)

% % ROI with a response to stimulation
% % must have peak time around stimulation?
%
% for iROI=1:length(ShortstimPeaks)
%     rc_str2(iROI)= ~isempty(strfind(ShortstimPeaks{iROI,14},'RCaMP'));
%     gc_str2(iROI)= ~isempty(strfind(ShortstimPeaks{iROI,14},'GCaMP'));
% end
%
% RCaMP_Peaks=ShortstimPeaks(rc_str2',:);
% GCaMP_Peaks=ShortstimPeaks(gc_str2',:);
%
% % NOTE: peak times have NOT be corrected for the baseline
%
% for iROI=1:length(RCaMP_Peaks)
%     respPTIdx1(iROI)=~isempty(find(RCaMP_Peaks{iROI,5}>5 && RCaMP_Peaks{iROI,5}<(NPTWindow+5)));%);
% end
% RrespPT=RCaMP_Peaks(respPTIdx1',:); % responding neurons
%
% for iROI=1:length(GCaMP_Peaks)
%     respPTIdx2(iROI)=~isempty(find(GCaMP_Peaks{iROI,5}>5 && GCaMP_Peaks{iROI,5}<(APTWindow+5)));%);
% end
% GrespPT=GCaMP_Peaks(respPTIdx2',:); % responding astrocytes
%
% %% RCaMP vs GCaMP plots of responding ROIs based on peak time
%
% PT_RCaMP_traces=[];
% % exact traces for ROIs responding based on peak time
% for xROI=1:length(RrespPT)
%     CurrentROI=RrespPT(xROI,23);
%     for iROI=1:length(Shortstim)
%         r_roi=strfind(Shortstim{iROI,15},CurrentROI);
%         if ~isempty(r_roi)
%             tempY= Shortstim{iROI,8};
%             if length(tempY)>590
%                 tempY=tempY(1:nframes);
%                 PT_RCaMP_traces= horzcat(PT_RCaMP_traces, tempY);
%             end
%         end
%     end
% end
% PT_RCaMP_mean=mean(PT_RCaMP_traces');
%
% PT_GCaMP_traces=[];
% % exact traces for ROIs responding based on peak time
% for xROI=1:length(GrespPT)
%     CurrentROI=GrespPT(xROI,23);
%     for iROI=1:length(Shortstim)
%         r_roi=strfind(Shortstim{iROI,15},CurrentROI);
%         if ~isempty(r_roi)
%             tempY= Shortstim{iROI,8};
%             if length(tempY)>590
%                 tempY=tempY(1:nframes);
%                 PT_GCaMP_traces= horzcat(PT_GCaMP_traces, tempY);
%             end
%         end
%     end
% end
% PT_GCaMP_mean=mean(PT_GCaMP_traces');
%
% % mean trace of astrocytes
% figure ('name', 'Overlaid traces: All RCaMP responding PT ROIs')
% hold on
% axis off
% for xROI= 1:size(PT_RCaMP_traces,2)
%     tempY = PT_RCaMP_traces(:,xROI);
%     if length(tempY)>590
%         tempY=tempY(1:nframes);
%         grey = [0.8,0.8,0.8];
%         plot(TimeX,tempY,'Color',grey,'LineWidth',0.01);
%     end
% end
% plot(TimeX, PT_RCaMP_mean, 'Color', 'k','LineWidth',1);
%
% figure ('name', 'Overlaid traces: All GCaMP responding PT ROIs')
% hold on
% axis off
% for xROI= 1:size(PT_GCaMP_traces,2)
%     tempY = PT_GCaMP_traces(:,xROI);
%     if length(tempY)>590
%         tempY=tempY(1:nframes);
%         grey = [0.8,0.8,0.8];
%         plot(TimeX,tempY,'Color',grey,'LineWidth',0.01);
%     end
% end
% plot(TimeX, PT_GCaMP_mean, 'Color', 'k','LineWidth',1);

%% how similar are responding groups?

% RrespondingROIs=intersect(RrespOT(:,15), RrespPT(:,23)); %
%
% GrespondingROIs=intersect(GrespOT(:,15), GrespPT(:,23)); %more similar




%% fast ROIs defined by onset time or AUC?

figure('name','1 sec auc vs onset time, astrocytes and neurons with stim onset times');
hold on
scatter(cell2mat(RrespOT(:,16)), cell2mat(RrespOT(:,17)));
scatter(cell2mat(GrespOT(:,17)), cell2mat(GrespOT(:,18)));

figure('name','1 sec auc vs onset time');
hold on
scatter(cell2mat(Shortstim(:,16)), cell2mat(Shortstim(:,17)));

% figure('name','20 sec auc vs onset time');
% hold on
% scatter(cell2mat(Shortstim(:,16)), cell2mat(Shortstim(:,18)));

for iROI=1:length(GrespOT)
    fastIdx(iROI)=~isempty(find(GrespOT{iROI,17}>0 && GrespOT{iROI,17}<=2));
    slowIdx(iROI)=~isempty(find(GrespOT{iROI,17}>2));
end
fastAC=GrespOT(fastIdx',:);
slowAC=GrespOT(slowIdx',:);

fastAC_traces=[];
slowAC_traces=[];
figure ('name', 'fast AC with onset times in first 1 sec')
hold on
axis off
for xROI= 1:size(fastAC,1)
    tempY = fastAC{xROI,8};
    if length(tempY)>590
        tempY=tempY(1:nframes);
        grey = [0.8,0.8,0.8];
        plot(TimeX,tempY,'Color',grey,'LineWidth',0.01);
    end
    fastAC_traces= horzcat(fastAC_traces, tempY);
end
fastAC_mean= mean(fastAC_traces');
plot(TimeX, fastAC_mean, 'Color', 'k','LineWidth',1);
%plot(TimeX, OT_RCaMP_mean, 'Color', 'r','LineWidth',1);
for xROI= 1:size(slowAC,1)
    tempY = slowAC{xROI,8};
    tempY=tempY(1:nframes);
    slowAC_traces= horzcat(slowAC_traces, tempY);
end
slowAC_mean= mean(slowAC_traces');
plot(TimeX, slowAC_mean, 'Color', 'b','LineWidth',1);

% mean onsets
MeanfastACOnset=mean(cell2mat(fastAC(:,16)));

% neurons with onset in first 1 sec of stim
for iROI=1:length(RrespOT)
    fastIdx2(iROI)=~isempty(find(RrespOT{iROI,17}>0 && RrespOT{iROI,17}<=1));
end
fastN=RrespOT(fastIdx2',:);
MeanNeuronOnset=mean(cell2mat(fastN(:,16)));


%% Data for heat maps

% sort astrocytes by onset time
[~, OT_GCidx] = sort([GrespOT{:,16}], 'ascend');
GrespOT_sort=GrespOT(OT_GCidx,:);

% sort neuronss by onset time
[~, OT_RCidx] = sort([RrespOT{:,16}], 'ascend');
RrespOT_sort=RrespOT(OT_RCidx,:);

% astrocyte table
AC_traces=[];
for xROI= 1:length(GrespOT_sort)
    tempOnset = GrespOT_sort{xROI,16}; % individual onset time
    tempY = GrespOT_sort{xROI,8}; % individual trace
    
    tempY=tempY(1:round(25*FrameRate)); % only first 15 seconds of trial (plus 5 sec baseline)
    
    temp=vertcat(tempOnset,tempY);
    AC_traces=vertcat(AC_traces,temp');
end

% neuronal table
N_traces=[];
for xROI= 1:length(RrespOT_sort)
    tempOnset = RrespOT_sort{xROI,16}; % individual onset time
    tempY = RrespOT_sort{xROI,8}; % individual trace
    
    tempY=tempY(1:round(25*FrameRate)); % only first 15 seconds of trial (plus 5 sec baseline)
    
    temp=vertcat(tempOnset,tempY);
    N_traces=vertcat(N_traces,temp');
end


xlswrite(AstrocyteExcelFile, AC_traces)
xlswrite(NeuronalExcelFile, N_traces)

%%
figure ('name', 'fast AC with onset times in first 1 sec-long stim')
hold on
axis off
xlim([-1 25]);
for xROI= 1:size(fastAC,1)
    tempY = fastAC{xROI,8};
    if length(tempY)>590
        tempY=tempY(1:nframes);
        grey = [0.8,0.8,0.8];
        plot(TimeX,tempY,'Color',grey,'LineWidth',0.01);
    end
end
plot(TimeX, fastAC_mean, 'Color', 'k','LineWidth',1);
rectangle('Position', [5 -1 8 10])
plot([-1 -1],[0 1], 'k','LineWidth', 1)

%%
figure('name', 'Lck long stim all means')
hold on
axis off
xlim([0 35]);
ylim([-0.5 2.5]);
plot(TimeX, smooth(fastAC_mean,10), 'Color', 'k','LineWidth',1);
%rectangle('Position', [5 -0.3 8 2.5])
%plot(TimeX, OT_RCaMP_mean, 'Color', 'r','LineWidth',1);
%plot(TimeX, slowAC_mean, 'Color', 'b','LineWidth',1);
%plot([-1 -1],[0 1], 'k','LineWidth', 1)

%% shaded error bar with
green=[(27/255) (120/255) (55/255)];
purple=[(123/255) (50/255) (148/255)];
blue= [(0/255) (114/255) (178/255)];

% SEM calculations
fastAC_SDTrace = std(fastAC_traces');
fastAC_SEM=fastAC_SDTrace/sqrt(size(fastAC_traces,2));

slowAC_SDTrace = std(slowAC_traces');
slowAC_SEM=slowAC_SDTrace/sqrt(size(slowAC_traces,2));

% RC_SDTrace = std(OT_RCaMP_traces');
% RC_SEM=RC_SDTrace/sqrt(size(OT_RCaMP_traces,2));

figure('name', 'Lck short stim all means- plus SEM')
hold on
axis off
xlim([-1 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

%ylim([-0.2 3]);
lineProps.col = {blue};
mseb(TimeX,slowAC_mean,slowAC_SEM,lineProps)
lineProps.col = {green};
mseb(TimeX,(fastAC_mean+0.75),fastAC_SEM,lineProps)
% lineProps.col = {purple};
% mseb(TimeX,(OT_RCaMP_mean+2),RC_SEM,lineProps)

% H2=shadedErrorBar(TimeX,(fastAC_mean+0.75),fastAC_SEM)%,green,1);
% H3=shadedErrorBar(TimeX,(OT_RCaMP_mean+2),RC_SEM)%,purple,1);
rectangle('Position', [5 -0.3 8 4])
plot([-0.1 -0.1],[0 1], 'k','LineWidth', 1)


figure('name', 'Lck short stim all means- plus SD')
hold on
axis off
xlim([-1 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';
%ylim([-0.2 3]);

lineProps.col = {blue};
mseb(TimeX,slowAC_mean,slowAC_SDTrace,lineProps)
lineProps.col = {green};
mseb(TimeX,(fastAC_mean+3),fastAC_SDTrace,lineProps)
lineProps.col = {purple};
mseb(TimeX,(OT_RCaMP_mean+7),RC_SDTrace,lineProps)
rectangle('Position', [5 -0.3 8 10])
plot([-0.1 -0.1],[0 1], 'k','LineWidth', 1)


%% traces for fast endfeet and fast processes

% Group Responders by ROIType
Resp_EF_traces=[];
Resp_processes_traces=[];

%find ROITypes
for xROI= 1:length(fastAC)
    EF_str= strfind(fastAC{xROI, 13},'Endfeet');
    P_str= strfind(fastAC{xROI, 13},'Process');
    
    tempY= fastAC{xROI,8};
    tempY=tempY(3:round(25*FrameRate));
    if ~isempty(EF_str)
        Resp_EF_traces= horzcat(Resp_EF_traces, tempY);
    elseif ~isempty(P_str)
        Resp_processes_traces= horzcat(Resp_processes_traces, tempY);
    end
end

% means and SDs

%endfeet
RespEFmeanTrace = mean(Resp_EF_traces,2)';
RespEFSEMTrace = std(Resp_EF_traces')/sqrt(size(Resp_EF_traces,2));

%processes
RespPmeanTrace = mean(Resp_processes_traces,2)';
RespPSEMTrace = std(Resp_processes_traces')/sqrt(size(Resp_processes_traces,2));


figure('name', 'Lck FAST endfeet vs processes- plus SEM')
hold on
axis off
xlim([-1 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

%ylim([-0.2 3]);
lineProps.col = {blue};
mseb(TimeX(3:round(25*FrameRate)),RespEFmeanTrace,RespEFSEMTrace,lineProps)
lineProps.col = {green};
mseb(TimeX(3:round(25*FrameRate)),(RespPmeanTrace+2.5),RespPSEMTrace,lineProps)
rectangle('Position', [5 -1 8 5.5])
plot([-0.1 -0.1],[0 1], 'k','LineWidth', 1)




%% Responding Neurons and Astrocytes based on onset times

% ROI with a response to stimulation
% must have onset time during stimulation?

% neurons with onset between 7 and 8 s

% stupid subsetting

for iROI=1:length(RCaMP)
    respOTIdx1(iROI)=~isempty(find(RCaMP{iROI,16}<6));
end

RCaMP2=RCaMP(respOTIdx1,:);

OT_RCaMP_traces=[];
for iROI=1:length(RCaMP2)
    respOTIdx3(iROI)=~isempty(find(RCaMP2{iROI,16}>5));
    tempY= RCaMP2{iROI,8};
    if length(tempY)>590
        tempY=tempY(1:nframes); %(stimwindow*FrameRate));
        if respOTIdx3(iROI)
            OT_RCaMP_traces= horzcat(OT_RCaMP_traces, tempY);
        end
    end
end

RrespOT=RCaMP2(respOTIdx3',:); % responding neurons
OT_RCaMP_mean= mean(OT_RCaMP_traces');

x1 = -2;
y1 = -2;
txt1=num2str(length(RrespOT));

figure ('name', 'Overlaid traces: RCaMP ROIs with onset between 5 and 6')
hold on
axis off
for xROI= 1:size(OT_RCaMP_traces,2)
    tempY = OT_RCaMP_traces(:,xROI);
    if length(tempY)>590
        tempY=tempY(1:nframes);
        grey = [0.8,0.8,0.8];
        plot(TimeX,tempY,'Color',grey,'LineWidth',0.01);
    end
end
plot(TimeX, OT_RCaMP_mean, 'Color', 'k','LineWidth',1);
rectangle('Position', [5 -1 8 5.5])
text(x1,y1,txt1)
