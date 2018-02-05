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


% %% compare peak onsets for each RCaMP and GCaMP ROIs
%
% % find astrocytes and neurons that are closest in time
% % find astrocytes and neurons that are closest in space
%
% %if ~exist(saveFiles1, 'file')
% tic
%
% Lck_traces_nonNaNs=~cellfun(@isnan,Lck_traces(:,14));
% Lck_traces=Lck_traces(Lck_traces_nonNaNs,:);
%
% AvsN_SpaceOnsets=[];
% Trials= unique(Lck_traces(:,18));
% Drug=unique(Lck_traces(:,17));
% Condition={'Nostim','Stim'};
%
% % onset times comparisons
% for iDrug=1:length(Drug)
%     CurrentDrug=Drug(iDrug);
%     matchingDrugIdx = find(~cellfun(@isempty, regexp(Lck_traces(:,17), CurrentDrug)));
%     DrugData = Lck_traces(matchingDrugIdx,:);
%
%     for iCond=1:length(Condition)
%         CurrentCondition=Condition(iCond);
%         matchingCondIdx = find(~cellfun(@isempty, regexp(DrugData(:,6), CurrentCondition)));
%         ConditionData = DrugData(matchingCondIdx,:);
%
%         for itrial=1:length(Trials)
%             CurrentTrial=Trials(itrial);
%             % Find the idx of paths matching trial
%             matchingTrialIdx = find(~cellfun(@isempty, regexp(ConditionData(:,18), CurrentTrial)));
%             TrialData = ConditionData(matchingTrialIdx,:);
%
%
%             % find neuronal onset times in this trial
%             NeuronalIdx = find(~cellfun(@isempty, regexp(TrialData(:,3), 'RCaMP')));
%             NeuronalData=TrialData(NeuronalIdx,:);
%             % remove the ones that have NaN onsets (or no detectable signal
%             % onsets in the onset window)
%             Neuro_nonNaNs=~cellfun(@isnan,NeuronalData(:,14));
%             NeuronalData=NeuronalData(Neuro_nonNaNs,:);
%
%
%             % find similar astrocyte onset times
%             AstroIdx = find(~cellfun(@isempty, regexp(TrialData(:,3), 'GCaMP')));
%             AstroData=TrialData(AstroIdx,:);
%             % remove the ones that have NaN onsets (or no detectable signal
%             % onsets in the onset window)
%             Astro_nonNaNs=~cellfun(@isnan,AstroData(:,14));
%             AstroData=AstroData(Astro_nonNaNs,:);
%
%             if ~isempty(AstroData)
%                 if ~isempty(NeuronalData)
%
%                     % considering astrocytes first
%                     for nAstro= 1:size(AstroData,1)
%                         AOnset=AstroData{nAstro,14};
%
%                         for nNeuro=1:size(NeuronalData,1)
%                             NOnsets{nNeuro,1}=NeuronalData{nNeuro,14};
%
%                             % masks for ROI distance calculations
%                             %create a binary image with each pair of ROIs
%                             % for the first ROI
%                             if isnumeric(NeuronalData{nNeuro,10})
%                                 Image1=zeros(128,128);
%                                 Image1(NeuronalData{nNeuro,10})=1;
%                                 %Image1=im2bw(Image1);
%                             elseif islogical(NeuronalData{nNeuro,10})
%                                 Image1= double(NeuronalData{nNeuro,10});
%                             else
%                                 Image1=[];
%                             end
%
%                             % for the second ROI
%                             if isnumeric(AstroData{nAstro,10})
%                                 Image2=zeros(128,128);
%                                 Image2(AstroData{nAstro,10})=1;
%                                 Image2=im2bw(Image2);
%                             elseif islogical(AstroData{nAstro,10})
%                                 Image2= double(AstroData{nAstro,10});
%                             else
%                                 Image2=[];
%                             end
%
%                             Mask=Image1+Image2;
%                             Mask=im2bw(Mask);
%
%                             % find the minimium distance between the edges of the two ROIs
%                             %Pythaogrean theorem method
%
%                             % Define object boundaries
%                             boundaries = bwboundaries(Mask);
%                             numberOfBoundaries = size(boundaries, 1);
%                             if numberOfBoundaries==1
%                                 distance{nNeuro,1} = 0;
%                             elseif numberOfBoundaries>1
%                                 boundary1 = boundaries{1};
%                                 boundary2 = boundaries{2};
%                                 boundary1x = boundary1(:, 2);
%                                 boundary1y = boundary1(:, 1);
%                                 minDistance= zeros(length(boundary2),1);
%                                 for k = 1 : length(boundary2)
%                                     boundary2x = boundary2(k, 2);
%                                     boundary2y = boundary2(k, 1);
%                                     % For this blob, compute distances from boundaries to edge.
%                                     allDistances = sqrt((boundary1x - boundary2x).^2 + (boundary1y - boundary2y).^2);
%                                     % Find closest point, min distance.
%                                     [minDistance(k), indexOfMin] = min(allDistances);
%                                 end
%                                 % Find the overall min distance
%                                 distance{nNeuro,1} = (min(minDistance)*cell2mat(NeuronalData(nNeuro,11)));
%                             else
%                                 distance{nNeuro,1} = [];
%                             end
%                         end
%
%                         %find the neuron that is closest in space
%                         N_spaceMin=min(cell2mat(distance));
%                         N_spaceMinIdx = find([distance{:}] == N_spaceMin);
%                         NSpaceOnset=NOnsets(N_spaceMinIdx,1);
%
%                         % find the astrocyte ROI area
%                         %                 CurrentAstro=AstroData(nAstro,19);
%                         %                 matchingROIIdx= find(~cellfun(@isempty, regexp(Lck_traces(:,23), CurrentAstro)));
%                         %                 A_ROIarea = Lck_traces(matchingROIIdx(1,1),18);
%
%
%                         % preallocation
%                         AvsN_space=cell(1,16);
%
%                         for iSpace=1:length(NSpaceOnset)
%                             % find the spatial astrocyte ROI area
%                             %                     CurrentNeuro2=NeuronalData(N_spaceMinIdx(iSpace),19);
%                             %                     matchingROIIdx3= find(~cellfun(@isempty, regexp(ShortstimPeaks(:,23), CurrentNeuro2)));
%                             %                     N_ROIarea2 = ShortstimPeaks(matchingROIIdx3(1,1),18);
%
%                             AvsN_space(iSpace,1)=AstroData(nAstro,5); % animal
%                             AvsN_space(iSpace,2)=AstroData(nAstro,4); % spot
%                             AvsN_space(iSpace,3)=AstroData(nAstro,2); % trial
%                             AvsN_space(iSpace,4)=AstroData(nAstro,19); % unique astro ROI name
%                             AvsN_space(iSpace,5)=AstroData(nAstro,20); % astro ROI type
%                             AvsN_space(iSpace,6)=CurrentDrug; % astro ROI type
%                             AvsN_space(iSpace,7)=CurrentCondition; % astro ROI type
%                             %                     AvsN_space(iSpace,6)=A_ROIarea;
%
%                             AvsN_space(iSpace,8)=NeuronalData(N_spaceMinIdx(iSpace),19); % unique neuro ROI name
%                             AvsN_space(iSpace,9)=NeuronalData(N_spaceMinIdx(iSpace),20); % neuro ROI type
%                             %                     AvsN_space(iSpace,9)=N_ROIarea2; % astrocyte ROI area
%                             AvsN_space(iSpace,10) =distance(N_spaceMinIdx(iSpace),1);
%                             AvsN_space{iSpace,11}=AOnset; % astro onset time
%                             AvsN_space(iSpace,12)= NSpaceOnset(iSpace);% neuro onset time
%                             AvsN_space{iSpace,13}=AOnset-cell2mat(NSpaceOnset(iSpace)); %astrocytes-neurons
%                             AvsN_space{iSpace,14}=AstroData{nAstro,9}; %astro trace
%                             AvsN_space{iSpace,15}=NeuronalData{N_spaceMinIdx(iSpace),9}; %neuronal trace
%                         end
%
%                         % concatenate all data into a big matrix
%                         AvsN_SpaceOnsets=vertcat(AvsN_SpaceOnsets,AvsN_space);
%
%                         clear boundaries boundary1 boundary2 boundary1x boundary1y boundary2x boundary2y
%                         clear minDistance indexofMin distance AvsN_space NOnsets
%                     end
%
%                 end
%             end
%         end
%     end
% end
% toc
%
% % only consider ROIs that are within -20 and 20 s
%
% AvsN_SpaceOnsetsIdx=find(cell2mat(AvsN_SpaceOnsets(:,13))>=-12 & cell2mat(AvsN_SpaceOnsets(:,13))<=12);
% AvsN_SpaceOnsets=AvsN_SpaceOnsets(AvsN_SpaceOnsetsIdx,:);
%
% Anames={'Animal','Spot','Trial','A_ROI','A_ROIType',...
%     'Drug', 'Condition', 'N_ROI','N_ROIType','Distance','AOnset','NOnset','TimeDiff',...
%     'Atrace','Ntrace'};
%
% AvsN_SpaceOnsets2=vertcat(Anames,AvsN_SpaceOnsets);
% cell2csv('D:\Data\GCaMP_RCaMP\Revision\Lck_GCaMP\FilesforR\Pharmacology\AvsN_SpaceOnsets.csv',AvsN_SpaceOnsets2);

%% Split data for each drug

% make a table for each drug
for iROI=1:length(AvsN_SpaceOnsets)
    con_str(iROI)= ~isempty(strfind(AvsN_SpaceOnsets{iROI,6},'Control'));
    DSP_str(iROI)= ~isempty(strfind(AvsN_SpaceOnsets{iROI,6},'DSP'));
    atr_str(iROI)= ~isempty(strfind(AvsN_SpaceOnsets{iROI,6},'Atropine'));
    Met_str(iROI)= ~isempty(strfind(AvsN_SpaceOnsets{iROI,6},'Metergoline'));
    Traz_str(iROI)= ~isempty(strfind(AvsN_SpaceOnsets{iROI,6},'Trazodone'));
    Praz_str(iROI)= ~isempty(strfind(AvsN_SpaceOnsets{iROI,6},'Prazosin'));
end

Control=AvsN_SpaceOnsets(con_str',:);
DSP4=AvsN_SpaceOnsets(DSP_str',:);
Atropine=AvsN_SpaceOnsets(atr_str',:);
Metergoline=AvsN_SpaceOnsets(Met_str',:);
Trazodone=AvsN_SpaceOnsets(Traz_str',:);
Prazosin=AvsN_SpaceOnsets(Praz_str',:);


%% sort by distance- some graphs with fast/delayed colour code, some graphs without

[~, Dis_idx] = sort([AvsN_SpaceOnsets{:,10}], 'ascend');
AvsN_SpaceOnsets=AvsN_SpaceOnsets(Dis_idx,:);
clearvars Dis_idx

Medians=[];
Drug=unique(Lck_traces(:,17));
Condition={'Stim'};
for iDrug=1:length(Drug)
    CurrentDrug=Drug(iDrug);
    matchingDrugIdx = find(~cellfun(@isempty, regexp(AvsN_SpaceOnsets(:,6), CurrentDrug)));
    DrugData = AvsN_SpaceOnsets(matchingDrugIdx,:);
    
    for iCond=1:length(Condition)
        CurrentCondition=Condition{iCond};
        matchingCondIdx = find(~cellfun(@isempty, regexp(DrugData(:,7), CurrentCondition)));
        ConditionData = DrugData(matchingCondIdx,:);
        %
        figure('name', 'Raster plot AvsN_SpaceOnsets- sorted by distance & coloured by group')
        hold on
        set(gca,'ytick',[])
        set(gca,'YColor',get(gcf,'Color'))
        for iComp=1:length(ConditionData)
            if  ConditionData{iComp,11}<1
                scatter(ConditionData{iComp,13}, iComp, 5, 'filled','b')
            elseif (ConditionData{iComp,11}>=1 &&  ConditionData{iComp,11}<12)
                scatter(ConditionData{iComp,13}, iComp, 5, 'filled','b')
            end
            xlim([-12 12])
            
        end
        plot([0 0],[0 length(ConditionData)], 'r--','LineWidth', 1)
        xlabel('time from neuronal event')
        
        figure('name', 'histogram: A vs N- closest in space')
        histogram(cell2mat(ConditionData(:,13)),158, 'Normalization','pdf')
        xlim([-20 20])
        ylim([0 0.3])
        xlabel('time from neuronal event')
        
       Median_subset(1,1)= ConditionData(1,6);
       Median_subset(1,2)= ConditionData(1,7);
       Median_subset{1,3}=median(cell2mat(ConditionData(:,13)));
       
       Medians=vertcat(Medians,Median_subset);
        
    end
end




