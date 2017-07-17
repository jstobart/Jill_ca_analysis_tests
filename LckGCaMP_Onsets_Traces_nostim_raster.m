
clearvars
close all

%% Load data

% time windows based on stimulation
NOnsetWindow= 8; %1 for short stim % neuronal onset times
AOnsetWindow= 15; % 5 for short stim % astrocyte onset times
Fast_AOnsetWindow=1;
Fast_NOnsetWindow=1;
NPTWindow= 9; % one second longer than stimulation for peak times
APTWindow= 12; % astrocyte longer than stimulation for peak times

% save files names
saveFiles1='E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\Lck_nostim_firstonset_comparisons.mat';
saveFiles2='E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\Lck_nostim_firstonset_comparisons.csv';
saveFiles3= 'E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\Lck_nostim_onset&AUC.csv';

%peak data
load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_2D_nostim_28_04_2017.mat');
%load('D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_2D_longstim_28_04_2017.mat');

% Load trace data
load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_nostim_28_04_2017.mat');
%load('D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_longstim_28_04_2017.mat');


ShortstimPeaks=AllData2(2:end,:);

% fix single digit trial names
for x=1:length(ShortstimPeaks)
    number_pos=regexp(ShortstimPeaks{x,12},'[0-9]');
    if length(number_pos)<2
        trialname=ShortstimPeaks{x,12};
        number=str2double(trialname(number_pos));
        ShortstimPeaks{x,12}=strcat('trial',num2str(number,'%02d'));
    end
end

Shortstim=All_traces;

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
        ShortstimPeaks{iROI,21}='Neuron';
    elseif ~isempty(D_str)
        ShortstimPeaks{iROI,21}='Dendrite';
    elseif ~isempty(r_str)
        P_str= strfind(ShortstimPeaks{iROI, 3},'GCaMP');
        if ~isempty(P_str)
            ShortstimPeaks{iROI,21}='Process';
        else
            ShortstimPeaks{iROI,21}='Dendrite';
        end
    elseif ~isempty(EF_str)
        ShortstimPeaks{iROI,21}='Endfeet';
    end
    
    % make new unique trial names
    ShortstimPeaks{iROI,22}=strcat(ShortstimPeaks{iROI,13},'_',ShortstimPeaks{iROI,15},'_',ShortstimPeaks{iROI,12});
    % unique ROI names
    ShortstimPeaks{iROI,23}=strcat(ShortstimPeaks{iROI,13},'_',ShortstimPeaks{iROI,15},'_',ShortstimPeaks{iROI,12}, '_',ShortstimPeaks{iROI,10});
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
FrameRate=11.84;
nframes=592;
TimeX(1:nframes) = (1:nframes)/FrameRate;

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
        Shortstim{iROI,13}='Neuron';
    elseif ~isempty(D_str)
        Shortstim{iROI,13}='Dendrite';
    elseif ~isempty(r_str)
        P_str= strfind(Shortstim{iROI, 3},'GCaMP');
        if ~isempty(P_str)
            Shortstim{iROI,13}='Process';
        else
            Shortstim{iROI,13}='Dendrite';
        end
    elseif ~isempty(EF_str)
        Shortstim{iROI,13}='Endfeet';
    end
    
    % make new unique trial names
    Shortstim{iROI,14}=strcat(Shortstim{iROI,5},'_',Shortstim{iROI,4},'_',Shortstim{iROI,2});
    % unique ROI names
    Shortstim{iROI,15}=strcat(Shortstim{iROI,5},'_',Shortstim{iROI,4},'_',Shortstim{iROI,2}, '_',Shortstim{iROI,1});
end

for iROI=1:length(Shortstim)
    %if a process ROI is overlapping with a soma or process, exclude it
    %from data
    nonOverlapIdx2(iROI)=~ischar(Shortstim{iROI,12});
end

% remove overlapping processes
Shortstim = Shortstim(nonOverlapIdx2',:);



% %% Calculate the first peak onset time and AUC after stim
% 
% baselineCorrectedTime=TimeX-5;
% 
% % peak onsets and AUC in the first second after stim for each ROI
% for iROI= 1:length(Shortstim)
%     trace=Shortstim{iROI,8};
%     %first 1 sec after stim onset
%     x1=round(FrameRate*5);
%     x2=round(FrameRate*6);
%     x3= round(FrameRate*10);
%     x4= round(FrameRate*15);
%     % onset time
%     if size(trace,1)>590
%         Onsets=find_first_onset_time(baselineCorrectedTime(10:end), trace(10:592,:),2.5,2);
%         if isempty(Onsets)
%             Onsets=nan(1,1);
%         end
%         Shortstim{iROI, 16}= Onsets;
%         % AUC
%         Shortstim{iROI,17}=trapz(trace(x1:x2));
%         Shortstim{iROI,18}=trapz(trace(x3:x4));
%     else
%         Shortstim{iROI,16}=NaN;
%         Shortstim{iROI,17}=NaN;
%         Shortstim{iROI,18}=NaN;
%     end
% end
% 
% % table for importing into R
% if ~exist('saveFiles3','file')
%     ShortstimSave=Shortstim;
%     ShortstimSave(:,8)=[];
%     ShortstimSave(:,8)=[];
%     ShortstimSave(:,8)=[];
%     names2={'ROI','Trial','Channel','Spot','Animal', 'Condition','depth','PixelSize','Overlap',...
%         'ROIType','Spot_trial','ROIs_trial','OnsetTime','TraceAUC1','TraceAUC10'};
%     ShortstimSave2=vertcat(names2, ShortstimSave);
%     cell2csv(saveFiles3,ShortstimSave2);
% end
% 
% 
% 
% %% compare peak onsets for each RCaMP and GCaMP ROIs
% 
% % compare ALL ROIS to ALL ROIS
% 
% if ~exist(saveFiles1, 'file')
%     tic
%     TimeComparisons=[];
%     Trials= unique(Shortstim(:,14));
%     
%     % onset times comparisons
%     for itrial=1:length(Trials)
%         CurrentTrial=Trials(itrial);
%         
%         % Find the idx of paths matching trial
%         matchingTrialIdx = find(~cellfun(@isempty, regexp(Shortstim(:,14), CurrentTrial)));
%         TrialData = Shortstim(matchingTrialIdx,:);
%         
%         
%         % find neuronal onset times in this trial
%         NeuronalIdx = find(~cellfun(@isempty, regexp(TrialData(:,3), 'RCaMP')));
%         NeuronalData=TrialData(NeuronalIdx,:);
%         
%         
%         % find similar astrocyte onset times
%         AstroIdx = find(~cellfun(@isempty, regexp(TrialData(:,3), 'GCaMP')));
%         AstroData=TrialData(AstroIdx,:);
%         
%         
%         parfor nNeuro=1:size(NeuronalData,1)
%             for nAstro= 1:size(AstroData,1)
%                 NOnset=NeuronalData{nNeuro,16};
%                 AOnset=AstroData{nAstro,16};
%                 
%                 if ~isempty(NOnset) || ~isnan(NOnset)
%                     if ~isempty(AOnset) || ~isnan(AOnset)
%                         
%                         % masks for ROI distance calculations
%                         %create a binary image with each pair of ROIs
%                         % for the first ROI
%                         if isnumeric(NeuronalData{nNeuro,10})
%                             Image1=zeros(128,128);
%                             Image1(NeuronalData{nNeuro,10})=1;
%                             %Image1=im2bw(Image1);
%                         elseif islogical(NeuronalData{nNeuro,10})
%                             Image1= double(NeuronalData{nNeuro,10});
%                         else
%                             Image1=[];
%                         end
%                         
%                         % for the second ROI
%                         if isnumeric(AstroData{nAstro,10})
%                             Image2=zeros(128,128);
%                             Image2(AstroData{nAstro,10})=1;
%                             Image2=im2bw(Image2);
%                         elseif islogical(AstroData{nAstro,10})
%                             Image2= double(AstroData{nAstro,10});
%                         else
%                             Image2=[];
%                         end
%                         
%                         Mask=Image1+Image2;
%                         Mask=im2bw(Mask);
%                         
%                         % find the minimium distance between the edges of the two ROIs
%                         %Pythaogrean theorem method
%                         
%                         % Define object boundaries
%                         boundaries = bwboundaries(Mask);
%                         numberOfBoundaries = size(boundaries, 1);
%                         if numberOfBoundaries==1
%                             distance = 0;
%                         elseif numberOfBoundaries>1
%                             boundary1 = boundaries{1};
%                             boundary2 = boundaries{2};
%                             boundary1x = boundary1(:, 2);
%                             boundary1y = boundary1(:, 1);
%                             minDistance= zeros(length(boundary2),1);
%                             for k = 1 : length(boundary2)
%                                 boundary2x = boundary2(k, 2);
%                                 boundary2y = boundary2(k, 1);
%                                 % For this blob, compute distances from boundaries to edge.
%                                 allDistances = sqrt((boundary1x - boundary2x).^2 + (boundary1y - boundary2y).^2);
%                                 % Find closest point, min distance.
%                                 [minDistance(k), indexOfMin] = min(allDistances);
%                             end
%                             % Find the overall min distance
%                             distance = (min(minDistance)*cell2mat(AstroData(nAstro,11)));
%                         else
%                             distance = [];
%                         end
%                         
%                         % find the neuronal ROI area
%                         CurrentNeuron=NeuronalData(nNeuro,15);
%                         matchingROIIdx= find(~cellfun(@isempty, regexp(ShortstimPeaks(:,23), CurrentNeuron)));
%                         N_ROIarea = ShortstimPeaks(matchingROIIdx(1,1),18);
%                         
%                         % find the neuronal ROI area
%                         CurrentAstro=AstroData(nAstro,15);
%                         matchingROIIdx2= find(~cellfun(@isempty, regexp(ShortstimPeaks(:,23), CurrentAstro)));
%                         A_ROIarea = ShortstimPeaks(matchingROIIdx2(1,1),18);
%                         
%                         for iN=1:length(NOnset)
%                             OnsetTimeComparisons=cell(length(AOnset),17);
%                             for iA=1:length(AOnset)
%                                 NeuronalOnset=NOnset(1,iN); % neuronal peak onset
%                                 AstrocyteOnset=AOnset(1,iA); % astrocyte peak onset
%                                 TimeDiff= AstrocyteOnset-NeuronalOnset;
%                                 TimeDiff2 = NeuronalOnset-AstrocyteOnset
%                                 %if TimeDiff>=-10 && TimeDiff<=10
%                                 
%                                 % generate data table with onset comparisons
%                                 OnsetTimeComparisons(iA,1)=NeuronalData(nNeuro,5); % animal
%                                 OnsetTimeComparisons(iA,2)=NeuronalData(nNeuro,4); % spot
%                                 OnsetTimeComparisons(iA,3)=NeuronalData(nNeuro,2); % trial
%                                 OnsetTimeComparisons(iA,4)=NeuronalData(nNeuro,15); % unique neuronal ROI name
%                                 OnsetTimeComparisons(iA,5)=NeuronalData(nNeuro,13); % neuron ROI type
%                                 OnsetTimeComparisons(iA,6)=N_ROIarea; % neuron ROI area
%                                 OnsetTimeComparisons(iA,7)=AstroData(nAstro,15); % unique astrocyte ROI name
%                                 OnsetTimeComparisons(iA,8)=AstroData(nAstro,13); % astrocyte ROI type
%                                 OnsetTimeComparisons(iA,9)=A_ROIarea; % astrocyte ROI area
%                                 OnsetTimeComparisons{iA,10} =distance;
%                                 OnsetTimeComparisons{iA,11} = numberOfBoundaries; % number of ROIs in the mask
%                                 
%                                 
%                                 % Onset times
%                                 OnsetTimeComparisons{iA,12}=strcat('NPeak',num2str(iN,'%02d')); % neuronal peak number
%                                 OnsetTimeComparisons{iA,13}=NeuronalOnset; % neuronal onset time
%                                 OnsetTimeComparisons{iA,14}=strcat('APeak',num2str(iA,'%02d')); % astrocyte peak number
%                                 OnsetTimeComparisons{iA,15}=AstrocyteOnset; % astrocyte onset time
%                                 
%                                 % difference between astrocyte onset and neuronal onset
%                                 OnsetTimeComparisons{iA,16}=TimeDiff; %astrocytes-neurons
%                                 OnsetTimeComparisons{iA,17}=strcat(OnsetTimeComparisons{iA,7},OnsetTimeComparisons{iA,14}); % astrocyte peak name
%                                 %                                 OnsetTimeComparisons{iA,18}=TimeDiff2; %neurons-astrocytes
%                                 %                                 OnsetTimeComparisons{iA,19}=NeuronalData{iN,8}; %neuronal trace
%                                 %                                 OnsetTimeComparisons{iA,20}=AstroData{iA,8};
%                             end
%                             
%                             TimeComparisons=vertcat(TimeComparisons,OnsetTimeComparisons); % concatenate all data into a big matrix
%                             
%                             %clear NeuronalData AstroData boundaries boundary1 boundary2 boundary1x boundary1y boundary2x boundary2y minDistance indexofMin OnsetTimeComparisons
%                         end
%                     end
%                 end
%             end
%         end
%     end
%     
%     save(saveFiles1, 'TimeComparisons','-v7.3');
%     
%     names={'Animal','Spot','Trial','N_ROI','N_ROIType','N_Area','A_ROI','A_ROIType',...
%         'A_Area','distance','ROInum','NPeak','N_Onset','APeak','A_Onset','TimeDiff',...
%         'A_peak_name'};
%     TimeComp2=vertcat(names,TimeComparisons);
%     cell2csv(saveFiles2, TimeComp2);
%     
%     toc
% else
%     load(saveFiles1);
% end
% 
% %% histograms
% figure('name', 'histogram: OnsetTime (A)- OnsetTime (N)')
% histogram(cell2mat(TimeComparisons(:,16)))
% 
% % ACIdx= find(~cellfun(@isempty, regexp(Shortstim(:,3), 'GCaMP')));
% % figure('name', 'Astrocyte onset times')
% % histogram(cell2mat(Shortstim(ACIdx,16)), 'binwidth',0.8450)
% %
% % NIdx= find(~cellfun(@isempty, regexp(Shortstim(:,3), 'RCaMP')));
% % figure('name', 'Neuronal onset times')
% % histogram(cell2mat(Shortstim(NIdx,16)), 'binwidth', 0.845)
% %
% % figure('name', 'Astrocyte AUC')
% % histogram(cell2mat(Shortstim(ACIdx,17)))
% %
% % figure('name', 'Neuronal AUC')
% % histogram(cell2mat(Shortstim(NIdx,17)))
% 
% %% compare peak onsets for each RCaMP and GCaMP ROIs
% 
% % find astrocytes and neurons that are closest in time
% % find astrocytes and neurons that are closest in space
% 
% %if ~exist(saveFiles1, 'file')
% tic
% NvsA_TimeOnsets=[];
% NvsA_SpaceOnsets=[];
% AvsN_TimeOnsets=[];
% AvsN_SpaceOnsets=[];
% Trials= unique(Shortstim(:,14));
% 
% % onset times comparisons
% for itrial=1:length(Trials)
%     CurrentTrial=Trials(itrial);
%     
%     % Find the idx of paths matching trial
%     matchingTrialIdx = find(~cellfun(@isempty, regexp(Shortstim(:,14), CurrentTrial)));
%     TrialData = Shortstim(matchingTrialIdx,:);
%     
%     
%     % find neuronal onset times in this trial
%     NeuronalIdx = find(~cellfun(@isempty, regexp(TrialData(:,3), 'RCaMP')));
%     NeuronalData=TrialData(NeuronalIdx,:);
%     % remove the ones that have NaN onsets (or no detectable signal
%     % onsets in the onset window)
%     Neuro_nonNaNs=~cellfun(@isnan,NeuronalData(:,16));
%     NeuronalData=NeuronalData(Neuro_nonNaNs,:);
%     
%     
%     % find similar astrocyte onset times
%     AstroIdx = find(~cellfun(@isempty, regexp(TrialData(:,3), 'GCaMP')));
%     AstroData=TrialData(AstroIdx,:);
%     % remove the ones that have NaN onsets (or no detectable signal
%     % onsets in the onset window)
%     Astro_nonNaNs=~cellfun(@isnan,AstroData(:,16));
%     AstroData=AstroData(Astro_nonNaNs,:);
%     
%     if ~isempty(AstroData)
%         if ~isempty(NeuronalData)
%             for nNeuro=1:size(NeuronalData,1)
%                 NOnset=NeuronalData{nNeuro,16};
%                 for nAstro= 1:size(AstroData,1)
%                     AOnsets{nAstro,1}=AstroData{nAstro,16};
%                     
%                     % masks for ROI distance calculations
%                     %create a binary image with each pair of ROIs
%                     % for the first ROI
%                     if isnumeric(NeuronalData{nNeuro,10})
%                         Image1=zeros(128,128);
%                         Image1(NeuronalData{nNeuro,10})=1;
%                         %Image1=im2bw(Image1);
%                     elseif islogical(NeuronalData{nNeuro,10})
%                         Image1= double(NeuronalData{nNeuro,10});
%                     else
%                         Image1=[];
%                     end
%                     
%                     % for the second ROI
%                     if isnumeric(AstroData{nAstro,10})
%                         Image2=zeros(128,128);
%                         Image2(AstroData{nAstro,10})=1;
%                         Image2=im2bw(Image2);
%                     elseif islogical(AstroData{nAstro,10})
%                         Image2= double(AstroData{nAstro,10});
%                     else
%                         Image2=[];
%                     end
%                     
%                     Mask=Image1+Image2;
%                     Mask=im2bw(Mask);
%                     
%                     % find the minimium distance between the edges of the two ROIs
%                     %Pythaogrean theorem method
%                     
%                     % Define object boundaries
%                     boundaries = bwboundaries(Mask);
%                     numberOfBoundaries = size(boundaries, 1);
%                     if numberOfBoundaries==1
%                         distance{nAstro,1} = 0;
%                     elseif numberOfBoundaries>1
%                         boundary1 = boundaries{1};
%                         boundary2 = boundaries{2};
%                         boundary1x = boundary1(:, 2);
%                         boundary1y = boundary1(:, 1);
%                         minDistance= zeros(length(boundary2),1);
%                         for k = 1 : length(boundary2)
%                             boundary2x = boundary2(k, 2);
%                             boundary2y = boundary2(k, 1);
%                             % For this blob, compute distances from boundaries to edge.
%                             allDistances = sqrt((boundary1x - boundary2x).^2 + (boundary1y - boundary2y).^2);
%                             % Find closest point, min distance.
%                             [minDistance(k), indexOfMin] = min(allDistances);
%                         end
%                         % Find the overall min distance
%                         distance{nAstro,1} = (min(minDistance)*cell2mat(AstroData(nAstro,11)));
%                     else
%                         distance{nAstro,1} = [];
%                     end
%                 end
%                 
%                 %find the astrocyte that is closest in time
%                 A_time=cell2mat(AOnsets);
%                 [Timediffs, A_timeIdx]=min(abs(A_time-NOnset));
%                 ATimeOnset=AOnsets(A_timeIdx,1);
%                 
%                 %find the astrocyte that is closest in space
%                 A_spaceMin=min(cell2mat(distance));
%                 A_spaceMinIdx = find([distance{:}] == A_spaceMin);
%                 ASpaceOnset=AOnsets(A_spaceMinIdx,1);
%                 
%                 % find the neuronal ROI area
%                 CurrentNeuron=NeuronalData(nNeuro,15);
%                 matchingROIIdx= find(~cellfun(@isempty, regexp(ShortstimPeaks(:,23), CurrentNeuron)));
%                 N_ROIarea = ShortstimPeaks(matchingROIIdx(1,1),18);
%                 
%                 
%                 % preallocation
%                 NvsA_time=cell(1,16);
%                 NvsA_space=cell(1,16);
%                 
%                 for iTime=1:length(ATimeOnset)
%                     % find the time astrocyte ROI area
%                     CurrentAstro1=AstroData(A_timeIdx(iTime),15);
%                     matchingROIIdx2= find(~cellfun(@isempty, regexp(ShortstimPeaks(:,23), CurrentAstro1)));
%                     A_ROIarea1 = ShortstimPeaks(matchingROIIdx2(1,1),18);
%                     
%                     
%                     NvsA_time(iTime,1)=NeuronalData(nNeuro,5); % animal
%                     NvsA_time(iTime,2)=NeuronalData(nNeuro,4); % spot
%                     NvsA_time(iTime,3)=NeuronalData(nNeuro,2); % trial
%                     NvsA_time(iTime,4)=NeuronalData(nNeuro,15); % unique neuronal ROI name
%                     NvsA_time(iTime,5)=NeuronalData(nNeuro,13); % neuron ROI type
%                     NvsA_time(iTime,6)=N_ROIarea;
%                     
%                     NvsA_time(iTime,7)=AstroData(A_timeIdx(iTime),15); % unique astrocyte ROI name
%                     NvsA_time(iTime,8)=AstroData(A_timeIdx(iTime),13); % astrocyte ROI type
%                     NvsA_time(iTime,9)=A_ROIarea1; % astrocyte ROI area
%                     NvsA_time(iTime,10) =distance(A_timeIdx(iTime),1);
%                     NvsA_time{iTime,11}=NOnset; % neuronal onset time
%                     NvsA_time(iTime,12)= ATimeOnset(iTime);% astrocyte onset time
%                     NvsA_time{iTime,13}=cell2mat(ATimeOnset(iTime))-NOnset; %astrocytes-neurons
%                     NvsA_time{iTime,14}=NeuronalData{nNeuro,8}; %neuronal trace
%                     NvsA_time{iTime,15}=AstroData{A_timeIdx(iTime),8}; %neuronal trace
%                 end
%                 
%                 for iSpace=1:length(ASpaceOnset)
%                     % find the spatial astrocyte ROI area
%                     CurrentAstro2=AstroData(A_spaceMinIdx(iSpace),15);
%                     matchingROIIdx3= find(~cellfun(@isempty, regexp(ShortstimPeaks(:,23), CurrentAstro2)));
%                     A_ROIarea2 = ShortstimPeaks(matchingROIIdx3(1,1),18);
%                     
%                     NvsA_space(iSpace,1)=NeuronalData(nNeuro,5); % animal
%                     NvsA_space(iSpace,2)=NeuronalData(nNeuro,4); % spot
%                     NvsA_space(iSpace,3)=NeuronalData(nNeuro,2); % trial
%                     NvsA_space(iSpace,4)=NeuronalData(nNeuro,15); % unique neuronal ROI name
%                     NvsA_space(iSpace,5)=NeuronalData(nNeuro,13); % neuron ROI type
%                     NvsA_space(iSpace,6)=N_ROIarea;
%                     
%                     NvsA_space(iSpace,7)=AstroData(A_spaceMinIdx(iSpace),15); % unique astrocyte ROI name
%                     NvsA_space(iSpace,8)=AstroData(A_spaceMinIdx(iSpace),13); % astrocyte ROI type
%                     NvsA_space(iSpace,9)=A_ROIarea2; % astrocyte ROI area
%                     NvsA_space(iSpace,10) =distance(A_spaceMinIdx(iSpace),1);
%                     NvsA_space{iSpace,11}=NOnset; % neuronal onset time
%                     NvsA_space(iSpace,12)= ASpaceOnset(iSpace);% astrocyte onset time
%                     NvsA_space{iSpace,13}=cell2mat(ASpaceOnset(iSpace))-NOnset; %astrocytes-neurons
%                     NvsA_space{iSpace,14}=NeuronalData{nNeuro,8}; %neuronal trace
%                     NvsA_space{iSpace,15}=AstroData{A_spaceMinIdx(iSpace),8}; %neuronal trace
%                 end
%                 
%                 % concatenate all data into a big matrix
%                 NvsA_TimeOnsets=vertcat(NvsA_TimeOnsets,NvsA_time);
%                 NvsA_SpaceOnsets=vertcat(NvsA_SpaceOnsets,NvsA_space);
%                 
%                 
%                 clear boundaries boundary1 boundary2 boundary1x boundary1y boundary2x boundary2y
%                 clear minDistance indexofMin distance NvsA_time NvsA_space AOnsets
%                 
%             end
%             
%             
%             
%             % considering astrocytes first
%             for nAstro= 1:size(AstroData,1)
%                 AOnset=AstroData{nAstro,16};
%                 
%                 for nNeuro=1:size(NeuronalData,1)
%                     NOnsets{nNeuro,1}=NeuronalData{nNeuro,16};
%                     
%                     % masks for ROI distance calculations
%                     %create a binary image with each pair of ROIs
%                     % for the first ROI
%                     if isnumeric(NeuronalData{nNeuro,10})
%                         Image1=zeros(128,128);
%                         Image1(NeuronalData{nNeuro,10})=1;
%                         %Image1=im2bw(Image1);
%                     elseif islogical(NeuronalData{nNeuro,10})
%                         Image1= double(NeuronalData{nNeuro,10});
%                     else
%                         Image1=[];
%                     end
%                     
%                     % for the second ROI
%                     if isnumeric(AstroData{nAstro,10})
%                         Image2=zeros(128,128);
%                         Image2(AstroData{nAstro,10})=1;
%                         Image2=im2bw(Image2);
%                     elseif islogical(AstroData{nAstro,10})
%                         Image2= double(AstroData{nAstro,10});
%                     else
%                         Image2=[];
%                     end
%                     
%                     Mask=Image1+Image2;
%                     Mask=im2bw(Mask);
%                     
%                     % find the minimium distance between the edges of the two ROIs
%                     %Pythaogrean theorem method
%                     
%                     % Define object boundaries
%                     boundaries = bwboundaries(Mask);
%                     numberOfBoundaries = size(boundaries, 1);
%                     if numberOfBoundaries==1
%                         distance{nNeuro,1} = 0;
%                     elseif numberOfBoundaries>1
%                         boundary1 = boundaries{1};
%                         boundary2 = boundaries{2};
%                         boundary1x = boundary1(:, 2);
%                         boundary1y = boundary1(:, 1);
%                         minDistance= zeros(length(boundary2),1);
%                         for k = 1 : length(boundary2)
%                             boundary2x = boundary2(k, 2);
%                             boundary2y = boundary2(k, 1);
%                             % For this blob, compute distances from boundaries to edge.
%                             allDistances = sqrt((boundary1x - boundary2x).^2 + (boundary1y - boundary2y).^2);
%                             % Find closest point, min distance.
%                             [minDistance(k), indexOfMin] = min(allDistances);
%                         end
%                         % Find the overall min distance
%                         distance{nNeuro,1} = (min(minDistance)*cell2mat(NeuronalData(nNeuro,11)));
%                     else
%                         distance{nNeuro,1} = [];
%                     end
%                 end
%                 
%                 %find the neuron that is closest in time
%                 N_time=cell2mat(NOnsets);
%                 [Timediffs2, N_timeIdx]=min(abs(N_time-AOnset));
%                 NTimeOnset=NOnsets(N_timeIdx,1);
%                 
%                 %find the astrocyte that is closest in space
%                 N_spaceMin=min(cell2mat(distance));
%                 N_spaceMinIdx = find([distance{:}] == N_spaceMin);
%                 NSpaceOnset=NOnsets(N_spaceMinIdx,1);
%                 
%                 % find the neuronal ROI area
%                 CurrentAstro=AstroData(nAstro,15);
%                 matchingROIIdx= find(~cellfun(@isempty, regexp(ShortstimPeaks(:,23), CurrentAstro)));
%                 A_ROIarea = ShortstimPeaks(matchingROIIdx(1,1),18);
%                 
%                 
%                 % preallocation
%                 AvsN_time=cell(1,16);
%                 AvsN_space=cell(1,16);
%                 
%                 for iTime=1:length(NTimeOnset)
%                     % find the time astrocyte ROI area
%                     CurrentNeuro1=NeuronalData(N_timeIdx(iTime),15);
%                     matchingROIIdx2= find(~cellfun(@isempty, regexp(ShortstimPeaks(:,23), CurrentNeuro1)));
%                     N_ROIarea1 = ShortstimPeaks(matchingROIIdx2(1,1),18);
%                     
%                     
%                     AvsN_time(iTime,1)=AstroData(nAstro,5); % animal
%                     AvsN_time(iTime,2)=AstroData(nAstro,4); % spot
%                     AvsN_time(iTime,3)=AstroData(nAstro,2); % trial
%                     AvsN_time(iTime,4)=AstroData(nAstro,15); % unique astrocyte ROI name
%                     AvsN_time(iTime,5)=AstroData(nAstro,13); % astrocyte ROI type
%                     AvsN_time(iTime,6)=A_ROIarea;
%                     
%                     AvsN_time(iTime,7)=NeuronalData(N_timeIdx(iTime),15); % unique neuronal ROI name
%                     AvsN_time(iTime,8)=NeuronalData(N_timeIdx(iTime),13); % astrocyte ROI type
%                     AvsN_time(iTime,9)=N_ROIarea1; % astrocyte ROI area
%                     AvsN_time(iTime,10) =distance(N_timeIdx(iTime),1);
%                     AvsN_time{iTime,11}=AOnset; % astrocyte onset time
%                     AvsN_time(iTime,12)= NTimeOnset(iTime);% neuronal onset time
%                     AvsN_time{iTime,13}=AOnset-cell2mat(NTimeOnset(iTime)); %astrocytes-neurons
%                     AvsN_time{iTime,14}=AstroData{nAstro,8}; %astrocyte trace
%                     AvsN_time{iTime,15}=NeuronalData{N_timeIdx(iTime),8}; %neuronal trace
%                 end
%                 
%                 for iSpace=1:length(NSpaceOnset)
%                     % find the spatial astrocyte ROI area
%                     CurrentNeuro2=NeuronalData(N_spaceMinIdx(iSpace),15);
%                     matchingROIIdx3= find(~cellfun(@isempty, regexp(ShortstimPeaks(:,23), CurrentNeuro2)));
%                     N_ROIarea2 = ShortstimPeaks(matchingROIIdx3(1,1),18);
%                     
%                     AvsN_space(iSpace,1)=AstroData(nAstro,5); % animal
%                     AvsN_space(iSpace,2)=AstroData(nAstro,4); % spot
%                     AvsN_space(iSpace,3)=AstroData(nAstro,2); % trial
%                     AvsN_space(iSpace,4)=AstroData(nAstro,15); % unique astro ROI name
%                     AvsN_space(iSpace,5)=AstroData(nAstro,13); % astro ROI type
%                     AvsN_space(iSpace,6)=A_ROIarea;
%                     
%                     AvsN_space(iSpace,7)=NeuronalData(N_spaceMinIdx(iSpace),15); % unique neuro ROI name
%                     AvsN_space(iSpace,8)=NeuronalData(N_spaceMinIdx(iSpace),13); % neuro ROI type
%                     AvsN_space(iSpace,9)=N_ROIarea2; % astrocyte ROI area
%                     AvsN_space(iSpace,10) =distance(N_spaceMinIdx(iSpace),1);
%                     AvsN_space{iSpace,11}=AOnset; % astro onset time
%                     AvsN_space(iSpace,12)= NSpaceOnset(iSpace);% neuro onset time
%                     AvsN_space{iSpace,13}=AOnset-cell2mat(NSpaceOnset(iSpace)); %astrocytes-neurons
%                     AvsN_space{iSpace,14}=AstroData{nAstro,8}; %astro trace
%                     AvsN_space{iSpace,15}=NeuronalData{N_spaceMinIdx(iSpace),8}; %neuronal trace
%                 end
%                 
%                 % concatenate all data into a big matrix
%                 AvsN_TimeOnsets=vertcat(AvsN_TimeOnsets,AvsN_time);
%                 AvsN_SpaceOnsets=vertcat(AvsN_SpaceOnsets,AvsN_space);
%                 
%                 
%                 clear boundaries boundary1 boundary2 boundary1x boundary1y boundary2x boundary2y
%                 clear minDistance indexofMin distance AvsN_time AvsN_space NOnsets
%             end
%             
%         end
%     end
% end
% 








%%


% save(saveFiles1, 'TimeComparisons','-v7.3');
%
% names={'Animal','Spot','Trial','N_ROI','N_ROIType','N_Area','A_ROI','A_ROIType',...
%     'A_Area','distance','ROInum','NPeak','N_Onset','APeak','A_Onset','TimeDiff',...
%     'A_peak_name'};
% TimeComp2=vertcat(names,TimeComparisons);
% cell2csv(saveFiles2, TimeComp2);
%
% toc
% else
%     load(saveFiles1);
%     end


%% only consider ROIs that are within -20 and 20 s

NvsA_TimeOnsetsIdx=find(cell2mat(NvsA_TimeOnsets(:,13))>=-20 & cell2mat(NvsA_TimeOnsets(:,13))<=20);
NvsA_TimeOnsets=NvsA_TimeOnsets(NvsA_TimeOnsetsIdx,:);

NvsA_SpaceOnsetsIdx=find(cell2mat(NvsA_SpaceOnsets(:,13))>=-20 & cell2mat(NvsA_SpaceOnsets(:,13))<=20);
NvsA_SpaceOnsets=NvsA_SpaceOnsets(NvsA_SpaceOnsetsIdx,:);

AvsN_TimeOnsetsIdx=find(cell2mat(AvsN_TimeOnsets(:,13))>=-20 & cell2mat(AvsN_TimeOnsets(:,13))<=20);
AvsN_TimeOnsets=AvsN_TimeOnsets(AvsN_TimeOnsetsIdx,:);

AvsN_SpaceOnsetsIdx=find(cell2mat(AvsN_SpaceOnsets(:,13))>=-20 & cell2mat(AvsN_SpaceOnsets(:,13))<=20);
AvsN_SpaceOnsets=AvsN_SpaceOnsets(AvsN_SpaceOnsetsIdx,:);

%% sort by astrocyte ROI area

[~, area_idx] = sort([NvsA_SpaceOnsets{:,9}], 'ascend');
NvsA_SpaceOnsetsAr=NvsA_SpaceOnsets(area_idx,:);

figure('name', 'Raster plot NvsA_SpaceOnsets- sorted by astrocyte ROI area')
hold on
set(gca,'ytick',[])
set(gca,'YColor',get(gcf,'Color'))
for iComp=1:length(NvsA_SpaceOnsetsAr)
    
    scatter(cell2mat(NvsA_SpaceOnsetsAr(iComp,13)), iComp, 5, 'filled','k')
    xlim([-20 20])
end
plot([0 0],[0 length(NvsA_SpaceOnsetsAr)], 'r--','LineWidth', 1)
xlabel('time from neuronal event')


%% add astrocyte peak information to onset comparisons cell arrays

[~,~,ib] = intersect(AvsN_SpaceOnsets(:,4),ShortstimPeaks(:,23));% astrocyte ROIs from comparisons vs peak info ROIs

peakInfo=ShortstimPeaks(ib,:);

for iROI=1:length(peakInfo)
    SpaceIdx = find(~cellfun(@isempty, regexp(AvsN_SpaceOnsets(:,4), peakInfo(iROI,23))));
    AvsN_SpaceOnsets(SpaceIdx,16)=peakInfo(iROI,1); % astrocyte amplitude
    AvsN_SpaceOnsets(SpaceIdx,17)=peakInfo(iROI,5); % astrocyte peak time
    AvsN_SpaceOnsets(SpaceIdx,18)=num2cell(cell2mat(peakInfo(iROI,3))*2); % astrocyte duration
end

[~,~,ia] = intersect(AvsN_TimeOnsets(:,4),ShortstimPeaks(:,23));% astrocyte ROIs from comparisons vs peak info ROIs

peakInfo2=ShortstimPeaks(ia,:);

for iROI=1:length(peakInfo2)
    TimeIdx = find(~cellfun(@isempty, regexp(AvsN_TimeOnsets(:,4), peakInfo2(iROI,23))));
    AvsN_TimeOnsets(TimeIdx,16)=peakInfo2(iROI,1); % astrocyte amplitude
    AvsN_TimeOnsets(TimeIdx,17)=peakInfo2(iROI,5); % astrocyte peak time
    AvsN_TimeOnsets(TimeIdx,18)=num2cell(cell2mat(peakInfo2(iROI,3))*2); % astrocyte duration
end

[~,~,ic] = intersect(NvsA_TimeOnsets(:,7),ShortstimPeaks(:,23));% astrocyte ROIs from comparisons vs peak info ROIs

peakInfo3=ShortstimPeaks(ic,:);

for iROI=1:length(peakInfo3)
    TimeIdx2 = find(~cellfun(@isempty, regexp(NvsA_TimeOnsets(:,7), peakInfo3(iROI,23))));
    NvsA_TimeOnsets(TimeIdx2,16)=peakInfo3(iROI,1); % astrocyte amplitude
    NvsA_TimeOnsets(TimeIdx2,17)=peakInfo3(iROI,5); % astrocyte peak time
    NvsA_TimeOnsets(TimeIdx2,18)=num2cell(cell2mat(peakInfo3(iROI,3))*2); % astrocyte duration
end


[~,~,id] = intersect(NvsA_SpaceOnsets(:,7),ShortstimPeaks(:,23));% astrocyte ROIs from comparisons vs peak info ROIs

peakInfo4=ShortstimPeaks(id,:);

for iROI=1:length(peakInfo4)
    SpaceIdx2 = find(~cellfun(@isempty, regexp(NvsA_SpaceOnsets(:,7), peakInfo4(iROI,23))));
    NvsA_SpaceOnsets(SpaceIdx2,16)=peakInfo4(iROI,1); % astrocyte amplitude
    NvsA_SpaceOnsets(SpaceIdx2,17)=peakInfo4(iROI,5); % astrocyte peak time
    NvsA_SpaceOnsets(SpaceIdx2,18)=num2cell(cell2mat(peakInfo4(iROI,3))*2); % astrocyte duration
end


%% unsorted rasters and histograms
figure('name', 'histogram: N vs A- closest in time')
h1=histogram(cell2mat(NvsA_TimeOnsets(:,13)),158, 'Normalization','pdf');
xlim([-20 20])
xlabel('time from neuronal event')

figure('name', 'Raster plot N vs A- closest in time- unsorted')
hold on
set(gca,'ytick',[]) % removes y axis
set(gca,'YColor',get(gcf,'Color'))
%axis off
for iComp=1:length(NvsA_TimeOnsets)
    scatter(cell2mat(NvsA_TimeOnsets(iComp,13)), iComp, 5,'k','o', 'filled')
    xlim([-20 20])

end
    plot([0 0],[0 length(NvsA_TimeOnsets)], 'r--','LineWidth', 1)
    xlabel('time from neuronal event')

%histogram

figure('name', 'histogram: N vs A- closest in space')
h2=histogram(cell2mat(NvsA_SpaceOnsets(:,13)),158, 'Normalization','pdf');
xlim([-20 20])
xlabel('time from neuronal event')

figure('name', 'Raster plot N vs A- closest in space- unsorted')
hold on
set(gca,'ytick',[])
set(gca,'YColor',get(gcf,'Color'))
for iComp=1:length(NvsA_SpaceOnsets)
    
    scatter(cell2mat(NvsA_SpaceOnsets(iComp,13)), iComp, 5, 'filled','k')
    xlim([-20 20])
end
plot([0 0],[0 length(NvsA_SpaceOnsets)], 'r--','LineWidth', 1)
xlabel('time from neuronal event')
%
figure('name', 'histogram: A vs N- closest in time')
histogram(cell2mat(AvsN_TimeOnsets(:,13)),158, 'Normalization','pdf')
xlim([-20 20])
xlabel('time from neuronal event')


figure('name', 'Raster plot A vs N- closest in time- unsorted')
hold on
set(gca,'ytick',[])
set(gca,'YColor',get(gcf,'Color'))
for iComp=1:length(AvsN_TimeOnsets)
    
    scatter(cell2mat(AvsN_TimeOnsets(iComp,13)), iComp, 5, 'filled','k')
    xlim([-20 20])
end
plot([0 0],[0 length(AvsN_TimeOnsets)], 'r--','LineWidth', 1)
xlabel('time from neuronal event')



figure('name', 'histogram: A vs N- closest in space')
histogram(cell2mat(AvsN_SpaceOnsets(:,13)),158, 'Normalization','pdf')
xlim([-20 20])
ylim([0 0.09])
xlabel('time from neuronal event')

figure('name', 'Raster plot A vs N- closest in space- unsorted')
hold on
set(gca,'ytick',[])
set(gca,'YColor',get(gcf,'Color'))
for iComp=1:length(AvsN_SpaceOnsets)
    
    scatter(cell2mat(AvsN_SpaceOnsets(iComp,13)), iComp, 5, 'filled','k')
    xlim([-20 20])
end
plot([0 0],[0 length(AvsN_SpaceOnsets)], 'r--','LineWidth', 1)
xlabel('time from neuronal event')

%% colour code raster by fast or delayed astrocytes

figure('name', 'Raster plot AvsN_SpaceOnsets- coloured by group')
hold on
set(gca,'ytick',[])
set(gca,'YColor',get(gcf,'Color'))
for iComp=1:length(AvsN_SpaceOnsets)
    if AvsN_SpaceOnsets{iComp,11}<1
    scatter(AvsN_SpaceOnsets{iComp,13}, iComp, 5, 'filled','b')
    elseif (AvsN_SpaceOnsets{iComp,11}>=1 && AvsN_SpaceOnsets{iComp,11}<12)
        scatter(AvsN_SpaceOnsets{iComp,13}, iComp, 5, 'filled','g')
    else
        scatter(AvsN_SpaceOnsets{iComp,13}, iComp, 5, 'filled','m')
    end
    xlim([-20 20])
    
end
plot([0 0],[0 length(AvsN_SpaceOnsets)], 'r--','LineWidth', 1)
xlabel('time from neuronal event')

figure('name', 'Raster plot AvsN_TimeOnsets- coloured by group')
hold on
set(gca,'ytick',[])
set(gca,'YColor',get(gcf,'Color'))
for iComp=1:length(AvsN_TimeOnsets)
    if AvsN_TimeOnsets{iComp,11}<1
    scatter(AvsN_TimeOnsets{iComp,13}, iComp, 5, 'filled','b')
    elseif (AvsN_TimeOnsets{iComp,11}>=1 && AvsN_TimeOnsets{iComp,11}<12)
        scatter(AvsN_TimeOnsets{iComp,13}, iComp, 5, 'filled','g')
    else
        scatter(AvsN_TimeOnsets{iComp,13}, iComp, 5, 'filled','m')
    end
    xlim([-20 20])
    
end
plot([0 0],[0 length(AvsN_TimeOnsets)], 'r--','LineWidth', 1)
xlabel('time from neuronal event')



figure('name', 'Raster plot NvsA_TimeOnsets- coloured by group')
hold on
set(gca,'ytick',[])
set(gca,'YColor',get(gcf,'Color'))
for iComp=1:length(NvsA_TimeOnsets)
    if NvsA_TimeOnsets{iComp,12}<1
    scatter(NvsA_TimeOnsets{iComp,13}, iComp, 5, 'filled','b')
    elseif (NvsA_TimeOnsets{iComp,12}>=1 && NvsA_TimeOnsets{iComp,12}<12)
        scatter(NvsA_TimeOnsets{iComp,13}, iComp, 5, 'filled','g')
    else
        scatter(NvsA_TimeOnsets{iComp,13}, iComp, 5, 'filled','m')
    end
    xlim([-20 20])
    
end
plot([0 0],[0 length(NvsA_TimeOnsets)], 'r--','LineWidth', 1)
xlabel('time from neuronal event')


figure('name', 'Raster plot NvsA_SpaceOnsets- coloured by group')
hold on
set(gca,'ytick',[])
set(gca,'YColor',get(gcf,'Color'))
for iComp=1:length(NvsA_SpaceOnsets)
    if NvsA_SpaceOnsets{iComp,12}<1
    scatter(NvsA_SpaceOnsets{iComp,13}, iComp, 5, 'filled','b')
    elseif (NvsA_SpaceOnsets{iComp,12}>=1 && NvsA_SpaceOnsets{iComp,12}<12)
        scatter(NvsA_SpaceOnsets{iComp,13}, iComp, 5, 'filled','g')
    else
        scatter(NvsA_SpaceOnsets{iComp,13}, iComp, 5, 'filled','m')
    end
    xlim([-20 20])
    
end
plot([0 0],[0 length(NvsA_SpaceOnsets)], 'r--','LineWidth', 1)
xlabel('time from neuronal event')


%% sort by distance

[~, Dis_idx] = sort([NvsA_SpaceOnsets{:,10}], 'ascend');
NvsA_SpaceOnsetsD=NvsA_SpaceOnsets(Dis_idx,:);
clearvars Dis_idx

figure('name', 'Raster plot NvsA_SpaceOnsets- sorted by distance & coloured by group')
subplot(1,2,1);
hold on
set(gca,'ytick',[])
set(gca,'YColor',get(gcf,'Color'))
for iComp=1:length(NvsA_SpaceOnsetsD)
    if NvsA_SpaceOnsetsD{iComp,12}<1
    scatter(NvsA_SpaceOnsetsD{iComp,13}, iComp, 5, 'filled','b')
    elseif (NvsA_SpaceOnsetsD{iComp,12}>=1 && NvsA_SpaceOnsetsD{iComp,12}<12)
        scatter(NvsA_SpaceOnsetsD{iComp,13}, iComp, 5, 'filled','g')
    else
        scatter(NvsA_SpaceOnsetsD{iComp,13}, iComp, 5, 'filled','m')
    end
    xlim([-20 20])
    
end
plot([0 0],[0 length(NvsA_SpaceOnsetsD)], 'r--','LineWidth', 1)
xlabel('time from neuronal event')

subplot(1,2,2);
colormap(jet(1000))
clims=[0 30];
hold on
imagesc(cell2mat(NvsA_SpaceOnsetsD(1:length(NvsA_SpaceOnsetsD),10)),clims)
colorbar



[~, Dis_idx] = sort([AvsN_SpaceOnsets{:,10}], 'ascend');
AvsN_SpaceOnsetsD=AvsN_SpaceOnsets(Dis_idx,:);
clearvars Dis_idx

figure('name', 'Raster plot AvsN_SpaceOnsets- sorted by distance & coloured by group')
subplot(1,2,1);
hold on
set(gca,'ytick',[])
set(gca,'YColor',get(gcf,'Color'))
for iComp=1:length(AvsN_SpaceOnsetsD)
    if AvsN_SpaceOnsetsD{iComp,11}<1
    scatter(AvsN_SpaceOnsetsD{iComp,13}, iComp, 5, 'filled','b')
    elseif (AvsN_SpaceOnsetsD{iComp,11}>=1 && AvsN_SpaceOnsetsD{iComp,11}<12)
        scatter(AvsN_SpaceOnsetsD{iComp,13}, iComp, 5, 'filled','g')
    else
        scatter(AvsN_SpaceOnsetsD{iComp,13}, iComp, 5, 'filled','m')
    end
    xlim([-20 20])
    
end
plot([0 0],[0 length(AvsN_SpaceOnsetsD)], 'r--','LineWidth', 1)
xlabel('time from neuronal event')

subplot(1,2,2);
colormap(jet(1000))
clims=[0 30];
hold on
imagesc(cell2mat(AvsN_SpaceOnsetsD(1:length(AvsN_SpaceOnsetsD),10)),clims)
colorbar



[~, Dis_idx] = sort([AvsN_TimeOnsets{:,10}], 'ascend');
AvsN_TimeOnsetsD=AvsN_TimeOnsets(Dis_idx,:);
clearvars Dis_idx

figure('name', 'Raster plot AvsN_TimeOnsets- sorted by distance & coloured by group')
subplot(1,2,1);
hold on
set(gca,'ytick',[])
set(gca,'YColor',get(gcf,'Color'))
for iComp=1:length(AvsN_TimeOnsetsD)
    if AvsN_TimeOnsetsD{iComp,11}<1
    scatter(AvsN_TimeOnsetsD{iComp,13}, iComp, 5, 'filled','b')
    elseif (AvsN_TimeOnsetsD{iComp,11}>=1 && AvsN_TimeOnsetsD{iComp,11}<12)
        scatter(AvsN_TimeOnsetsD{iComp,13}, iComp, 5, 'filled','g')
    else
        scatter(AvsN_TimeOnsetsD{iComp,13}, iComp, 5, 'filled','m')
    end
    xlim([-20 20])
    
end
plot([0 0],[0 length(AvsN_TimeOnsetsD)], 'r--','LineWidth', 1)
xlabel('time from neuronal event')

subplot(1,2,2);
colormap(jet(1000))
clims=[0 60];
hold on
imagesc(cell2mat(AvsN_TimeOnsetsD(1:length(AvsN_TimeOnsetsD),10)),clims)
colorbar




[~, Dis_idx] = sort([NvsA_TimeOnsets{:,10}], 'ascend');
NvsA_TimeOnsetsD = NvsA_TimeOnsets(Dis_idx,:);
clearvars Dis_idx

figure('name', 'Raster plot NvsA_TimeOnsets- sorted by distance & coloured by group')
subplot(1,2,1);
hold on
set(gca,'ytick',[])
set(gca,'YColor',get(gcf,'Color'))
for iComp=1:length(NvsA_TimeOnsetsD)
    if NvsA_TimeOnsetsD{iComp,12}<1
    scatter(NvsA_TimeOnsetsD{iComp,13}, iComp, 5, 'filled','b')
    elseif (NvsA_TimeOnsetsD{iComp,12}>=1 && NvsA_TimeOnsetsD{iComp,12}<12)
        scatter(NvsA_TimeOnsetsD{iComp,13}, iComp, 5, 'filled','g')
    else
        scatter(NvsA_TimeOnsetsD{iComp,13}, iComp, 5, 'filled','m')
    end
    xlim([-20 20])
    
end
plot([0 0],[0 length(NvsA_TimeOnsetsD)], 'r--','LineWidth', 1)
xlabel('time from neuronal event')

subplot(1,2,2);
colormap(jet(1000))
clims=[0 60];
hold on
imagesc(cell2mat(NvsA_TimeOnsetsD(1:length(NvsA_TimeOnsetsD),10)),clims)
colorbar



%% sort by timedifference

% [~, Dis_idx] = sort([NvsA_SpaceOnsets{:,13}], 'ascend');
% NvsA_SpaceOnsets=NvsA_SpaceOnsets(Dis_idx,:);
% 
% figure('name', 'Raster plot NvsA_SpaceOnsets- sorted by timedifference')
% hold on
% set(gca,'ytick',[])
% set(gca,'YColor',get(gcf,'Color'))
% for iComp=1:length(NvsA_SpaceOnsets)
%     
%     scatter(cell2mat(NvsA_SpaceOnsets(iComp,13)), iComp, 5, 'filled','k')
%     xlim([-20 20])
% end
% plot([0 0],[0 length(NvsA_SpaceOnsets)], 'r--','LineWidth', 1)
% xlabel('time from neuronal event')



%% sort by amplitude
[~, amp_idx] = sort([AvsN_SpaceOnsets{:,16}], 'ascend');
AvsN_SpaceOnsetsAmp=AvsN_SpaceOnsets(amp_idx,:);
clearvars amp_idx

figure('name', 'Raster plot AvsN_SpaceOnsets- sorted by AC amplitude & coloured by group')
subplot(1,2,1);
hold on
set(gca,'ytick',[])
set(gca,'YColor',get(gcf,'Color'))
for iComp=1:length(AvsN_SpaceOnsetsAmp)
    if AvsN_SpaceOnsetsAmp{iComp,11}<1
    scatter(AvsN_SpaceOnsetsAmp{iComp,13}, iComp, 5, 'filled','b')
    elseif (AvsN_SpaceOnsetsAmp{iComp,11}>=1 && AvsN_SpaceOnsetsAmp{iComp,11}<12)
        scatter(AvsN_SpaceOnsetsAmp{iComp,13}, iComp, 5, 'filled','g')
    else
        scatter(AvsN_SpaceOnsetsAmp{iComp,13}, iComp, 5, 'filled','m')
    end
    xlim([-20 20])
    
end
plot([0 0],[0 length(AvsN_SpaceOnsetsAmp)], 'r--','LineWidth', 1)
xlabel('time from neuronal event')

subplot(1,2,2);
colormap(jet(1000))
clims=[0 3];
hold on
imagesc(cell2mat(AvsN_SpaceOnsetsAmp(1:length(AvsN_SpaceOnsetsAmp),16)),clims)
colorbar



[~, amp_idx] = sort([NvsA_SpaceOnsets{:,16}], 'ascend');
NvsA_SpaceOnsetsAmp=NvsA_SpaceOnsets(amp_idx,:);
clearvars amp_idx

figure('name', 'Raster plot NvsA_SpaceOnsets- sorted by AC amplitude & coloured by group')
subplot(1,2,1);
hold on
set(gca,'ytick',[])
set(gca,'YColor',get(gcf,'Color'))
for iComp=1:length(NvsA_SpaceOnsetsAmp)
    if NvsA_SpaceOnsetsAmp{iComp,12}<1
    scatter(NvsA_SpaceOnsetsAmp{iComp,13}, iComp, 5, 'filled','b')
    elseif (NvsA_SpaceOnsetsAmp{iComp,12}>=1 && NvsA_SpaceOnsetsAmp{iComp,12}<12)
        scatter(NvsA_SpaceOnsetsAmp{iComp,13}, iComp, 5, 'filled','g')
    else
        scatter(NvsA_SpaceOnsetsAmp{iComp,13}, iComp, 5, 'filled','m')
    end
    xlim([-20 20])
    
end
plot([0 0],[0 length(NvsA_SpaceOnsetsAmp)], 'r--','LineWidth', 1)
xlabel('time from neuronal event')

subplot(1,2,2);
colormap(jet(1000))
clims=[0 3];
hold on
imagesc(cell2mat(NvsA_SpaceOnsetsAmp(1:length(NvsA_SpaceOnsetsAmp),16)),clims)
colorbar




[~, amp_idx] = sort([NvsA_TimeOnsets{:,16}], 'ascend');
NvsA_TimeOnsetsAmp=NvsA_TimeOnsets(amp_idx,:);
clearvars amp_idx

figure('name', 'Raster plot NvsA_TimeOnsets- sorted by AC amplitude & coloured by group')
subplot(1,2,1);
hold on
set(gca,'ytick',[])
set(gca,'YColor',get(gcf,'Color'))
for iComp=1:length(NvsA_TimeOnsetsAmp)
    if NvsA_TimeOnsetsAmp{iComp,12}<1
    scatter(NvsA_TimeOnsetsAmp{iComp,13}, iComp, 5, 'filled','b')
    elseif (NvsA_TimeOnsetsAmp{iComp,12}>=1 && NvsA_TimeOnsetsAmp{iComp,12}<12)
        scatter(NvsA_TimeOnsetsAmp{iComp,13}, iComp, 5, 'filled','g')
    else
        scatter(NvsA_TimeOnsetsAmp{iComp,13}, iComp, 5, 'filled','m')
    end
    xlim([-20 20])
    
end
plot([0 0],[0 length(NvsA_TimeOnsetsAmp)], 'r--','LineWidth', 1)
xlabel('time from neuronal event')

subplot(1,2,2);
colormap(jet(1000))
clims=[0 3];
hold on
imagesc(cell2mat(NvsA_TimeOnsetsAmp(1:length(NvsA_TimeOnsetsAmp),16)),clims)
colorbar



[~, amp_idx] = sort([NvsA_TimeOnsets{:,16}], 'ascend');
NvsA_TimeOnsetsAmp=NvsA_TimeOnsets(amp_idx,:);
clearvars amp_idx

figure('name', 'Raster plot NvsA_TimeOnsets- sorted by AC amplitude & coloured by group')
subplot(1,2,1);
hold on
set(gca,'ytick',[])
set(gca,'YColor',get(gcf,'Color'))
for iComp=1:length(NvsA_TimeOnsetsAmp)
    if NvsA_TimeOnsetsAmp{iComp,12}<1
    scatter(NvsA_TimeOnsetsAmp{iComp,13}, iComp, 5, 'filled','b')
    elseif (NvsA_TimeOnsetsAmp{iComp,12}>=1 && NvsA_TimeOnsetsAmp{iComp,12}<12)
        scatter(NvsA_TimeOnsetsAmp{iComp,13}, iComp, 5, 'filled','g')
    else
        scatter(NvsA_TimeOnsetsAmp{iComp,13}, iComp, 5, 'filled','m')
    end
    xlim([-20 20])
    
end
plot([0 0],[0 length(NvsA_TimeOnsetsAmp)], 'r--','LineWidth', 1)
xlabel('time from neuronal event')

subplot(1,2,2);
colormap(jet(1000))
clims=[0 3];
hold on
imagesc(cell2mat(NvsA_TimeOnsetsAmp(1:length(NvsA_TimeOnsetsAmp),16)),clims)
colorbar

%% Sort by signal duration

[~, dur_idx] = sort([AvsN_SpaceOnsets{:,18}], 'ascend');
AvsN_SpaceOnsetsDur=AvsN_SpaceOnsets(dur_idx,:);
clearvars dur_idx

figure('name', 'Raster plot AvsN_SpaceOnsets- sorted by astrocyte duration')
hold on
set(gca,'ytick',[])
set(gca,'YColor',get(gcf,'Color'))
for iComp=1:length(AvsN_SpaceOnsetsDur)
    if AvsN_SpaceOnsetsDur{iComp,11}<1
    scatter(AvsN_SpaceOnsetsDur{iComp,13}, iComp, 5, 'filled','b')
    elseif (AvsN_SpaceOnsetsDur{iComp,11}>=1 && AvsN_SpaceOnsetsDur{iComp,11}<12)
        scatter(AvsN_SpaceOnsetsDur{iComp,13}, iComp, 5, 'filled','g')
    else
        scatter(AvsN_SpaceOnsetsDur{iComp,13}, iComp, 5, 'filled','m')
    end
    xlim([-20 20])
    
end
plot([0 0],[0 length(AvsN_SpaceOnsetsDur)], 'r--','LineWidth', 1)
xlabel('time from neuronal event')

subplot(1,2,2);
colormap(jet(1000))
clims=[0 10];
hold on
imagesc(cell2mat(AvsN_SpaceOnsetsDur(1:length(AvsN_SpaceOnsetsDur),18)),clims)
colorbar


[~, dur_idx] = sort([NvsA_SpaceOnsets{:,18}], 'ascend');
NvsA_SpaceOnsetsDur=NvsA_SpaceOnsets(dur_idx,:);
clearvars dur_idx

figure('name', 'Raster plot NvsA_SpaceOnsets- sorted by astrocyte duration')
hold on
set(gca,'ytick',[])
set(gca,'YColor',get(gcf,'Color'))
for iComp=1:length(NvsA_SpaceOnsetsDur)
    if NvsA_SpaceOnsetsDur{iComp,12}<1
    scatter(NvsA_SpaceOnsetsDur{iComp,13}, iComp, 5, 'filled','b')
    elseif (NvsA_SpaceOnsetsDur{iComp,12}>=1 && NvsA_SpaceOnsetsDur{iComp,12}<12)
        scatter(NvsA_SpaceOnsetsDur{iComp,13}, iComp, 5, 'filled','g')
    else
        scatter(NvsA_SpaceOnsetsDur{iComp,13}, iComp, 5, 'filled','m')
    end
    xlim([-20 20])
    
end
plot([0 0],[0 length(NvsA_SpaceOnsetsDur)], 'r--','LineWidth', 1)
xlabel('time from neuronal event')

subplot(1,2,2);
colormap(jet(1000))
clims=[0 10];
hold on
imagesc(cell2mat(NvsA_SpaceOnsetsDur(1:length(NvsA_SpaceOnsetsDur),18)),clims)
colorbar


[~, dur_idx] = sort([NvsA_TimeOnsets{:,18}], 'ascend');
NvsA_TimeOnsetsDur=NvsA_TimeOnsets(dur_idx,:);
clearvars dur_idx

figure('name', 'Raster plot NvsA_TimeOnsets- sorted by astrocyte duration')
hold on
set(gca,'ytick',[])
set(gca,'YColor',get(gcf,'Color'))
for iComp=1:length(NvsA_TimeOnsetsDur)
    if NvsA_TimeOnsetsDur{iComp,12}<1
    scatter(NvsA_TimeOnsetsDur{iComp,13}, iComp, 5, 'filled','b')
    elseif (NvsA_TimeOnsetsDur{iComp,12}>=1 && NvsA_TimeOnsetsDur{iComp,12}<12)
        scatter(NvsA_TimeOnsetsDur{iComp,13}, iComp, 5, 'filled','g')
    else
        scatter(NvsA_TimeOnsetsDur{iComp,13}, iComp, 5, 'filled','m')
    end
    xlim([-20 20])
    
end
plot([0 0],[0 length(NvsA_TimeOnsetsDur)], 'r--','LineWidth', 1)
xlabel('time from neuronal event')

subplot(1,2,2);
colormap(jet(1000))
clims=[0 10];
hold on
imagesc(cell2mat(NvsA_TimeOnsetsDur(1:length(NvsA_TimeOnsetsDur),18)),clims)
colorbar



[~, dur_idx] = sort([AvsN_TimeOnsets{:,18}], 'ascend');
AvsN_TimeOnsetsDur=AvsN_TimeOnsets(dur_idx,:);
clearvars dur_idx

figure('name', 'Raster plot AvsN_TimeOnsets- sorted by astrocyte duration')
hold on
set(gca,'ytick',[])
set(gca,'YColor',get(gcf,'Color'))
for iComp=1:length(AvsN_TimeOnsetsDur)
    if AvsN_TimeOnsetsDur{iComp,11}<1
    scatter(AvsN_TimeOnsetsDur{iComp,13}, iComp, 5, 'filled','b')
    elseif (AvsN_TimeOnsetsDur{iComp,11}>=1 && AvsN_TimeOnsetsDur{iComp,11}<12)
        scatter(AvsN_TimeOnsetsDur{iComp,13}, iComp, 5, 'filled','g')
    else
        scatter(AvsN_TimeOnsetsDur{iComp,13}, iComp, 5, 'filled','m')
    end
    xlim([-20 20])
    
end
plot([0 0],[0 length(AvsN_TimeOnsetsDur)], 'r--','LineWidth', 1)
xlabel('time from neuronal event')

subplot(1,2,2);
colormap(jet(1000))
clims=[0 10];
hold on
imagesc(cell2mat(AvsN_TimeOnsetsDur(1:length(AvsN_TimeOnsetsDur),18)),clims)
colorbar


%% Save data
Anames={'Animal','Spot','Trial','A_ROI','A_ROIType',...
    'A_area', 'N_ROI', 'N_ROIType','N_area','Distance','AOnset','NOnset','TimeDiff',...
    'Atrace','Ntrace','amplitude','peakTime','Duration'};
Nnames={'Animal','Spot','Trial','A_ROI','A_ROIType',...
    'A_area', 'N_ROI', 'N_ROIType','N_area','Distance','NOnset','AOnset','TimeDiff',...
    'Atrace','Ntrace','amplitude','peakTime','Duration'};

AvsN_SpaceOnsets=vertcat(Anames,AvsN_SpaceOnsets);
AvsN_TimeOnsets=vertcat(Anames,AvsN_TimeOnsets);

NvsA_SpaceOnsets=vertcat(Nnames,NvsA_SpaceOnsets);
NvsA_TimeOnsets=vertcat(Nnames,NvsA_TimeOnsets);

cell2csv('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\OnsetTimeComps\Nostim\Nostim_AvsN_SpaceOnsets.csv',AvsN_SpaceOnsets);
cell2csv('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\OnsetTimeComps\Nostim\Nostim_AvsN_TimeOnsets.csv',AvsN_TimeOnsets);
cell2csv('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\OnsetTimeComps\Nostim\Nostim_NvsA_SpaceOnsets.csv',NvsA_SpaceOnsets);
cell2csv('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\OnsetTimeComps\Nostim\Nostim_NvsA_TimeOnsets.csv',NvsA_TimeOnsets);
 %% combine nostim and stim
% Stim=NvsA_SpaceOnsets;
% Stim(:,16)={'Stim'};
% Nostim=NvsA_SpaceOnsets;
% Nostim(:,16)={'Nostim'};
% 
% %%
% AllComb=vertcat(Stim,Nostim);
% [~, R_idx] = sort(AllComb(:,2));
% AllComb=AllComb(R_idx,:);
% 
% figure('name', 'histogram: Stim vs Nostim, NvsA_space')
% h1=histogram(cell2mat(Stim(:,13)),158, 'Normalization','pdf');
% hold on
% h3=histogram(cell2mat(Nostim(:,13)),158, 'Normalization','pdf');
% xlim([-20 20])
% xlabel('time from neuronal event')
% 
% figure('name', 'Raster plot Stim vs Nostim- NvsA_space')
% hold on
% set(gca,'ytick',[])
% set(gca,'YColor',get(gcf,'Color'))
% for iComp=1:length(Stim)
%     scatter(cell2mat(Stim(iComp,13)), iComp, 5, 'filled','b')
%     if iComp<=length(Nostim)
%     scatter(cell2mat(Nostim(iComp,13)), iComp, 5, 'filled','r')
%     end
%     xlim([-20 20])
% end
% plot([0 0],[0 length(Stim)], 'r--','LineWidth', 1)
% xlabel('time from neuronal event')
% 
