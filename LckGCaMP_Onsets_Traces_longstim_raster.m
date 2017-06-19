
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
saveFiles1='E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\Lck_longstim_firstonset_comparisons.mat';
saveFiles2='E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\Lck_longstim_firstonset_comparisons.csv';
saveFiles3= 'E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\Lck_longstim_onset&AUC.csv';

%peak data
load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_2D_longstim_28_04_2017.mat');
%load('D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_2D_longstim_28_04_2017.mat');

% Load trace data
load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_longstim_28_04_2017.mat');
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



%% Calculate the first peak onset time and AUC after stim

baselineCorrectedTime=TimeX-5;

% peak onsets and AUC in the first second after stim for each ROI
for iROI= 1:length(Shortstim)
    trace=Shortstim{iROI,8};
    %first 1 sec after stim onset
    x1=round(FrameRate*5);
    x2=round(FrameRate*6);
    x3= round(FrameRate*10);
    x4= round(FrameRate*15);
    % onset time
    if size(trace,1)>590
        Onsets=find_first_onset_time(baselineCorrectedTime(10:end), trace(10:592,:),2.5,2);
        if isempty(Onsets)
            Onsets=nan(1,1);
        end
        Shortstim{iROI, 16}= Onsets;
        % AUC
        Shortstim{iROI,17}=trapz(trace(x1:x2));
        Shortstim{iROI,18}=trapz(trace(x3:x4));
    else
        Shortstim{iROI,16}=NaN;
        Shortstim{iROI,17}=NaN;
        Shortstim{iROI,18}=NaN;
    end
end

% table for importing into R
if ~exist('saveFiles3','file')
    ShortstimSave=Shortstim;
    ShortstimSave(:,8)=[];
    ShortstimSave(:,8)=[];
    ShortstimSave(:,8)=[];
    names2={'ROI','Trial','Channel','Spot','Animal', 'Condition','depth','PixelSize','Overlap',...
        'ROIType','Spot_trial','ROIs_trial','OnsetTime','TraceAUC1','TraceAUC10'};
    ShortstimSave2=vertcat(names2, ShortstimSave);
    cell2csv(saveFiles3,ShortstimSave2);
end



%% compare peak onsets for each RCaMP and GCaMP ROIs

% compare ALL ROIS to ALL ROIS

if ~exist(saveFiles1, 'file')
    tic
    TimeComparisons=[];
    Trials= unique(Shortstim(:,14));
    
    % onset times comparisons
    for itrial=1:length(Trials)
        CurrentTrial=Trials(itrial);
        
        % Find the idx of paths matching trial
        matchingTrialIdx = find(~cellfun(@isempty, regexp(Shortstim(:,14), CurrentTrial)));
        TrialData = Shortstim(matchingTrialIdx,:);
        
        
        % find neuronal onset times in this trial
        NeuronalIdx = find(~cellfun(@isempty, regexp(TrialData(:,3), 'RCaMP')));
        NeuronalData=TrialData(NeuronalIdx,:);
        
        
        % find similar astrocyte onset times
        AstroIdx = find(~cellfun(@isempty, regexp(TrialData(:,3), 'GCaMP')));
        AstroData=TrialData(AstroIdx,:);
        
        
        parfor nNeuro=1:size(NeuronalData,1)
            for nAstro= 1:size(AstroData,1)
                NOnset=NeuronalData{nNeuro,16};
                AOnset=AstroData{nAstro,16};
                
                if ~isempty(NOnset) || ~isnan(NOnset)
                    if ~isempty(AOnset) || ~isnan(AOnset)
                        
                        % masks for ROI distance calculations
                        %create a binary image with each pair of ROIs
                        % for the first ROI
                        if isnumeric(NeuronalData{nNeuro,10})
                            Image1=zeros(128,128);
                            Image1(NeuronalData{nNeuro,10})=1;
                            %Image1=im2bw(Image1);
                        elseif islogical(NeuronalData{nNeuro,10})
                            Image1= double(NeuronalData{nNeuro,10});
                        else
                            Image1=[];
                        end
                        
                        % for the second ROI
                        if isnumeric(AstroData{nAstro,10})
                            Image2=zeros(128,128);
                            Image2(AstroData{nAstro,10})=1;
                            Image2=im2bw(Image2);
                        elseif islogical(AstroData{nAstro,10})
                            Image2= double(AstroData{nAstro,10});
                        else
                            Image2=[];
                        end
                        
                        Mask=Image1+Image2;
                        Mask=im2bw(Mask);
                        
                        % find the minimium distance between the edges of the two ROIs
                        %Pythaogrean theorem method
                        
                        % Define object boundaries
                        boundaries = bwboundaries(Mask);
                        numberOfBoundaries = size(boundaries, 1);
                        if numberOfBoundaries==1
                            distance = 0;
                        elseif numberOfBoundaries>1
                            boundary1 = boundaries{1};
                            boundary2 = boundaries{2};
                            boundary1x = boundary1(:, 2);
                            boundary1y = boundary1(:, 1);
                            minDistance= zeros(length(boundary2),1);
                            for k = 1 : length(boundary2)
                                boundary2x = boundary2(k, 2);
                                boundary2y = boundary2(k, 1);
                                % For this blob, compute distances from boundaries to edge.
                                allDistances = sqrt((boundary1x - boundary2x).^2 + (boundary1y - boundary2y).^2);
                                % Find closest point, min distance.
                                [minDistance(k), indexOfMin] = min(allDistances);
                            end
                            % Find the overall min distance
                            distance = (min(minDistance)*cell2mat(AstroData(nAstro,11)));
                        else
                            distance = [];
                        end
                        
                        % find the neuronal ROI area
                        CurrentNeuron=NeuronalData(nNeuro,15);
                        matchingROIIdx= find(~cellfun(@isempty, regexp(ShortstimPeaks(:,23), CurrentNeuron)));
                        N_ROIarea = ShortstimPeaks(matchingROIIdx(1,1),18);
                        
                        % find the neuronal ROI area
                        CurrentAstro=AstroData(nAstro,15);
                        matchingROIIdx2= find(~cellfun(@isempty, regexp(ShortstimPeaks(:,23), CurrentAstro)));
                        A_ROIarea = ShortstimPeaks(matchingROIIdx2(1,1),18);
                        
                        for iN=1:length(NOnset)
                            OnsetTimeComparisons=cell(length(AOnset),17);
                            for iA=1:length(AOnset)
                                NeuronalOnset=NOnset(1,iN); % neuronal peak onset
                                AstrocyteOnset=AOnset(1,iA); % astrocyte peak onset
                                TimeDiff= AstrocyteOnset-NeuronalOnset;
                                TimeDiff2 = NeuronalOnset-AstrocyteOnset
                                %if TimeDiff>=-10 && TimeDiff<=10
                                
                                % generate data table with onset comparisons
                                OnsetTimeComparisons(iA,1)=NeuronalData(nNeuro,5); % animal
                                OnsetTimeComparisons(iA,2)=NeuronalData(nNeuro,4); % spot
                                OnsetTimeComparisons(iA,3)=NeuronalData(nNeuro,2); % trial
                                OnsetTimeComparisons(iA,4)=NeuronalData(nNeuro,15); % unique neuronal ROI name
                                OnsetTimeComparisons(iA,5)=NeuronalData(nNeuro,13); % neuron ROI type
                                OnsetTimeComparisons(iA,6)=N_ROIarea; % neuron ROI area
                                OnsetTimeComparisons(iA,7)=AstroData(nAstro,15); % unique astrocyte ROI name
                                OnsetTimeComparisons(iA,8)=AstroData(nAstro,13); % astrocyte ROI type
                                OnsetTimeComparisons(iA,9)=A_ROIarea; % astrocyte ROI area
                                OnsetTimeComparisons{iA,10} =distance;
                                OnsetTimeComparisons{iA,11} = numberOfBoundaries; % number of ROIs in the mask
                                
                                
                                % Onset times
                                OnsetTimeComparisons{iA,12}=strcat('NPeak',num2str(iN,'%02d')); % neuronal peak number
                                OnsetTimeComparisons{iA,13}=NeuronalOnset; % neuronal onset time
                                OnsetTimeComparisons{iA,14}=strcat('APeak',num2str(iA,'%02d')); % astrocyte peak number
                                OnsetTimeComparisons{iA,15}=AstrocyteOnset; % astrocyte onset time
                                
                                % difference between astrocyte onset and neuronal onset
                                OnsetTimeComparisons{iA,16}=TimeDiff; %astrocytes-neurons
                                OnsetTimeComparisons{iA,17}=strcat(OnsetTimeComparisons{iA,7},OnsetTimeComparisons{iA,14}); % astrocyte peak name
                                %                                 OnsetTimeComparisons{iA,18}=TimeDiff2; %neurons-astrocytes
                                %                                 OnsetTimeComparisons{iA,19}=NeuronalData{iN,8}; %neuronal trace
                                %                                 OnsetTimeComparisons{iA,20}=AstroData{iA,8};
                            end
                            
                            TimeComparisons=vertcat(TimeComparisons,OnsetTimeComparisons); % concatenate all data into a big matrix
                            
                            %clear NeuronalData AstroData boundaries boundary1 boundary2 boundary1x boundary1y boundary2x boundary2y minDistance indexofMin OnsetTimeComparisons
                        end
                    end
                end
            end
        end
    end
    
    save(saveFiles1, 'TimeComparisons','-v7.3');
    
    names={'Animal','Spot','Trial','N_ROI','N_ROIType','N_Area','A_ROI','A_ROIType',...
        'A_Area','distance','ROInum','NPeak','N_Onset','APeak','A_Onset','TimeDiff',...
        'A_peak_name'};
    TimeComp2=vertcat(names,TimeComparisons);
    cell2csv(saveFiles2, TimeComp2);
    
    toc
else
    load(saveFiles1);
end

%% histograms
figure('name', 'histogram: OnsetTime (A)- OnsetTime (N)')
histogram(cell2mat(TimeComparisons(:,16)))

% ACIdx= find(~cellfun(@isempty, regexp(Shortstim(:,3), 'GCaMP')));
% figure('name', 'Astrocyte onset times')
% histogram(cell2mat(Shortstim(ACIdx,16)), 'binwidth',0.8450)
%
% NIdx= find(~cellfun(@isempty, regexp(Shortstim(:,3), 'RCaMP')));
% figure('name', 'Neuronal onset times')
% histogram(cell2mat(Shortstim(NIdx,16)), 'binwidth', 0.845)
%
% figure('name', 'Astrocyte AUC')
% histogram(cell2mat(Shortstim(ACIdx,17)))
%
% figure('name', 'Neuronal AUC')
% histogram(cell2mat(Shortstim(NIdx,17)))

%% compare peak onsets for each RCaMP and GCaMP ROIs

% find astrocytes and neurons that are closest in time
% find astrocytes and neurons that are closest in space

%if ~exist(saveFiles1, 'file')
tic
NvsA_TimeOnsets=[];
NvsA_SpaceOnsets=[];
AvsN_TimeOnsets=[];
AvsN_SpaceOnsets=[];
Trials= unique(Shortstim(:,14));

% onset times comparisons
for itrial=1:length(Trials)
    CurrentTrial=Trials(itrial);
    
    % Find the idx of paths matching trial
    matchingTrialIdx = find(~cellfun(@isempty, regexp(Shortstim(:,14), CurrentTrial)));
    TrialData = Shortstim(matchingTrialIdx,:);
    
    
    % find neuronal onset times in this trial
    NeuronalIdx = find(~cellfun(@isempty, regexp(TrialData(:,3), 'RCaMP')));
    NeuronalData=TrialData(NeuronalIdx,:);
    % remove the ones that have NaN onsets (or no detectable signal
    % onsets in the onset window)
    Neuro_nonNaNs=~cellfun(@isnan,NeuronalData(:,16));
    NeuronalData=NeuronalData(Neuro_nonNaNs,:);
    
    
    % find similar astrocyte onset times
    AstroIdx = find(~cellfun(@isempty, regexp(TrialData(:,3), 'GCaMP')));
    AstroData=TrialData(AstroIdx,:);
    % remove the ones that have NaN onsets (or no detectable signal
    % onsets in the onset window)
    Astro_nonNaNs=~cellfun(@isnan,AstroData(:,16));
    AstroData=AstroData(Astro_nonNaNs,:);
    
    if ~isempty(AstroData)
        for nNeuro=1:size(NeuronalData,1)
            NOnset=NeuronalData{nNeuro,16};
            for nAstro= 1:size(AstroData,1)
                AOnsets{nAstro,1}=AstroData{nAstro,16};
                
                % masks for ROI distance calculations
                %create a binary image with each pair of ROIs
                % for the first ROI
                if isnumeric(NeuronalData{nNeuro,10})
                    Image1=zeros(128,128);
                    Image1(NeuronalData{nNeuro,10})=1;
                    %Image1=im2bw(Image1);
                elseif islogical(NeuronalData{nNeuro,10})
                    Image1= double(NeuronalData{nNeuro,10});
                else
                    Image1=[];
                end
                
                % for the second ROI
                if isnumeric(AstroData{nAstro,10})
                    Image2=zeros(128,128);
                    Image2(AstroData{nAstro,10})=1;
                    Image2=im2bw(Image2);
                elseif islogical(AstroData{nAstro,10})
                    Image2= double(AstroData{nAstro,10});
                else
                    Image2=[];
                end
                
                Mask=Image1+Image2;
                Mask=im2bw(Mask);
                
                % find the minimium distance between the edges of the two ROIs
                %Pythaogrean theorem method
                
                % Define object boundaries
                boundaries = bwboundaries(Mask);
                numberOfBoundaries = size(boundaries, 1);
                if numberOfBoundaries==1
                    distance{nAstro,1} = 0;
                elseif numberOfBoundaries>1
                    boundary1 = boundaries{1};
                    boundary2 = boundaries{2};
                    boundary1x = boundary1(:, 2);
                    boundary1y = boundary1(:, 1);
                    minDistance= zeros(length(boundary2),1);
                    for k = 1 : length(boundary2)
                        boundary2x = boundary2(k, 2);
                        boundary2y = boundary2(k, 1);
                        % For this blob, compute distances from boundaries to edge.
                        allDistances = sqrt((boundary1x - boundary2x).^2 + (boundary1y - boundary2y).^2);
                        % Find closest point, min distance.
                        [minDistance(k), indexOfMin] = min(allDistances);
                    end
                    % Find the overall min distance
                    distance{nAstro,1} = (min(minDistance)*cell2mat(AstroData(nAstro,11)));
                else
                    distance{nAstro,1} = [];
                end
            end
            
            %find the astrocyte that is closest in time
            A_time=cell2mat(AOnsets);
            [Timediffs, A_timeIdx]=min(abs(A_time-NOnset));
            ATimeOnset=AOnsets(A_timeIdx,1);
            
            %find the astrocyte that is closest in space
            A_spaceMin=min(cell2mat(distance));
            A_spaceMinIdx = find([distance{:}] == A_spaceMin);
            ASpaceOnset=AOnsets(A_spaceMinIdx,1);
            
            % find the neuronal ROI area
            CurrentNeuron=NeuronalData(nNeuro,15);
            matchingROIIdx= find(~cellfun(@isempty, regexp(ShortstimPeaks(:,23), CurrentNeuron)));
            N_ROIarea = ShortstimPeaks(matchingROIIdx(1,1),18);
            
            
            % preallocation
            NvsA_time=cell(1,16);
            NvsA_space=cell(1,16);
            
            for iTime=1:length(ATimeOnset)
                % find the time astrocyte ROI area
                CurrentAstro1=AstroData(A_timeIdx(iTime),15);
                matchingROIIdx2= find(~cellfun(@isempty, regexp(ShortstimPeaks(:,23), CurrentAstro1)));
                A_ROIarea1 = ShortstimPeaks(matchingROIIdx2(1,1),18);
                
                
                NvsA_time(iTime,1)=NeuronalData(nNeuro,5); % animal
                NvsA_time(iTime,2)=NeuronalData(nNeuro,4); % spot
                NvsA_time(iTime,3)=NeuronalData(nNeuro,2); % trial
                NvsA_time(iTime,4)=NeuronalData(nNeuro,15); % unique neuronal ROI name
                NvsA_time(iTime,5)=NeuronalData(nNeuro,13); % neuron ROI type
                NvsA_time(iTime,6)=N_ROIarea;
                
                NvsA_time(iTime,7)=AstroData(A_timeIdx(iTime),15); % unique astrocyte ROI name
                NvsA_time(iTime,8)=AstroData(A_timeIdx(iTime),13); % astrocyte ROI type
                NvsA_time(iTime,9)=A_ROIarea1; % astrocyte ROI area
                NvsA_time(iTime,10) =distance(A_timeIdx(iTime),1);
                NvsA_time{iTime,11}=NOnset; % neuronal onset time
                NvsA_time(iTime,12)= ATimeOnset(iTime);% astrocyte onset time
                NvsA_time{iTime,13}=cell2mat(ATimeOnset(iTime))-NOnset; %astrocytes-neurons
                NvsA_time{iTime,14}=NeuronalData{nNeuro,8}; %neuronal trace
                NvsA_time{iTime,15}=AstroData{A_timeIdx(iTime),8}; %neuronal trace
            end
            
            for iSpace=1:length(ASpaceOnset)
                % find the spatial astrocyte ROI area
                CurrentAstro2=AstroData(A_spaceMinIdx(iSpace),15);
                matchingROIIdx3= find(~cellfun(@isempty, regexp(ShortstimPeaks(:,23), CurrentAstro2)));
                A_ROIarea2 = ShortstimPeaks(matchingROIIdx3(1,1),18);
                
                NvsA_space(iSpace,1)=NeuronalData(nNeuro,5); % animal
                NvsA_space(iSpace,2)=NeuronalData(nNeuro,4); % spot
                NvsA_space(iSpace,3)=NeuronalData(nNeuro,2); % trial
                NvsA_space(iSpace,4)=NeuronalData(nNeuro,15); % unique neuronal ROI name
                NvsA_space(iSpace,5)=NeuronalData(nNeuro,13); % neuron ROI type
                NvsA_space(iSpace,6)=N_ROIarea;
                
                NvsA_space(iSpace,7)=AstroData(A_spaceMinIdx(iSpace),15); % unique astrocyte ROI name
                NvsA_space(iSpace,8)=AstroData(A_spaceMinIdx(iSpace),13); % astrocyte ROI type
                NvsA_space(iSpace,9)=A_ROIarea2; % astrocyte ROI area
                NvsA_space(iSpace,10) =distance(A_spaceMinIdx(iSpace),1);
                NvsA_space{iSpace,11}=NOnset; % neuronal onset time
                NvsA_space(iSpace,12)= ASpaceOnset(iSpace);% astrocyte onset time
                NvsA_space{iSpace,13}=cell2mat(ASpaceOnset(iSpace))-NOnset; %astrocytes-neurons
                NvsA_space{iSpace,14}=NeuronalData{nNeuro,8}; %neuronal trace
                NvsA_space{iSpace,15}=AstroData{A_spaceMinIdx(iSpace),8}; %neuronal trace
            end
            
            % concatenate all data into a big matrix
            NvsA_TimeOnsets=vertcat(NvsA_TimeOnsets,NvsA_time);
            NvsA_SpaceOnsets=vertcat(NvsA_SpaceOnsets,NvsA_space);
            
            
            clear boundaries boundary1 boundary2 boundary1x boundary1y boundary2x boundary2y
            clear minDistance indexofMin distance NvsA_time NvsA_space AOnsets
            
        end
        
        
        
        % considering astrocytes first
        for nAstro= 1:size(AstroData,1)
            AOnset=AstroData{nAstro,16};
            
            for nNeuro=1:size(NeuronalData,1)
                NOnsets{nNeuro,1}=NeuronalData{nNeuro,16};
                
                % masks for ROI distance calculations
                %create a binary image with each pair of ROIs
                % for the first ROI
                if isnumeric(NeuronalData{nNeuro,10})
                    Image1=zeros(128,128);
                    Image1(NeuronalData{nNeuro,10})=1;
                    %Image1=im2bw(Image1);
                elseif islogical(NeuronalData{nNeuro,10})
                    Image1= double(NeuronalData{nNeuro,10});
                else
                    Image1=[];
                end
                
                % for the second ROI
                if isnumeric(AstroData{nAstro,10})
                    Image2=zeros(128,128);
                    Image2(AstroData{nAstro,10})=1;
                    Image2=im2bw(Image2);
                elseif islogical(AstroData{nAstro,10})
                    Image2= double(AstroData{nAstro,10});
                else
                    Image2=[];
                end
                
                Mask=Image1+Image2;
                Mask=im2bw(Mask);
                
                % find the minimium distance between the edges of the two ROIs
                %Pythaogrean theorem method
                
                % Define object boundaries
                boundaries = bwboundaries(Mask);
                numberOfBoundaries = size(boundaries, 1);
                if numberOfBoundaries==1
                    distance{nNeuro,1} = 0;
                elseif numberOfBoundaries>1
                    boundary1 = boundaries{1};
                    boundary2 = boundaries{2};
                    boundary1x = boundary1(:, 2);
                    boundary1y = boundary1(:, 1);
                    minDistance= zeros(length(boundary2),1);
                    for k = 1 : length(boundary2)
                        boundary2x = boundary2(k, 2);
                        boundary2y = boundary2(k, 1);
                        % For this blob, compute distances from boundaries to edge.
                        allDistances = sqrt((boundary1x - boundary2x).^2 + (boundary1y - boundary2y).^2);
                        % Find closest point, min distance.
                        [minDistance(k), indexOfMin] = min(allDistances);
                    end
                    % Find the overall min distance
                    distance{nNeuro,1} = (min(minDistance)*cell2mat(NeuronalData(nNeuro,11)));
                else
                    distance{nNeuro,1} = [];
                end
            end
            
            %find the neuron that is closest in time
            N_time=cell2mat(NOnsets);
            [Timediffs2, N_timeIdx]=min(abs(N_time-AOnset));
            NTimeOnset=NOnsets(N_timeIdx,1);
            
            %find the astrocyte that is closest in space
            N_spaceMin=min(cell2mat(distance));
            N_spaceMinIdx = find([distance{:}] == N_spaceMin);
            NSpaceOnset=NOnsets(N_spaceMinIdx,1);
            
            % find the neuronal ROI area
            CurrentAstro=AstroData(nAstro,15);
            matchingROIIdx= find(~cellfun(@isempty, regexp(ShortstimPeaks(:,23), CurrentAstro)));
            A_ROIarea = ShortstimPeaks(matchingROIIdx(1,1),18);
            
            
            % preallocation
            AvsN_time=cell(1,16);
            AvsN_space=cell(1,16);
            
            for iTime=1:length(NTimeOnset)
                % find the time astrocyte ROI area
                CurrentNeuro1=NeuronalData(N_timeIdx(iTime),15);
                matchingROIIdx2= find(~cellfun(@isempty, regexp(ShortstimPeaks(:,23), CurrentNeuro1)));
                N_ROIarea1 = ShortstimPeaks(matchingROIIdx2(1,1),18);
                
                
                AvsN_time(iTime,1)=AstroData(nAstro,5); % animal
                AvsN_time(iTime,2)=AstroData(nAstro,4); % spot
                AvsN_time(iTime,3)=AstroData(nAstro,2); % trial
                AvsN_time(iTime,4)=AstroData(nAstro,15); % unique astrocyte ROI name
                AvsN_time(iTime,5)=AstroData(nAstro,13); % astrocyte ROI type
                AvsN_time(iTime,6)=A_ROIarea;
                
                AvsN_time(iTime,7)=NeuronalData(N_timeIdx(iTime),15); % unique neuronal ROI name
                AvsN_time(iTime,8)=NeuronalData(N_timeIdx(iTime),13); % astrocyte ROI type
                AvsN_time(iTime,9)=N_ROIarea1; % astrocyte ROI area
                AvsN_time(iTime,10) =distance(N_timeIdx(iTime),1);
                AvsN_time{iTime,11}=AOnset; % astrocyte onset time
                AvsN_time(iTime,12)= NTimeOnset(iTime);% neuronal onset time
                AvsN_time{iTime,13}=AOnset-cell2mat(NTimeOnset(iTime)); %astrocytes-neurons
                AvsN_time{iTime,14}=AstroData{nAstro,8}; %astrocyte trace
                AvsN_time{iTime,15}=NeuronalData{N_timeIdx(iTime),8}; %neuronal trace
            end
            
            for iSpace=1:length(NSpaceOnset)
                % find the spatial astrocyte ROI area
                CurrentNeuro2=NeuronalData(N_spaceMinIdx(iSpace),15);
                matchingROIIdx3= find(~cellfun(@isempty, regexp(ShortstimPeaks(:,23), CurrentNeuro2)));
                N_ROIarea2 = ShortstimPeaks(matchingROIIdx3(1,1),18);
                
                AvsN_space(iSpace,1)=AstroData(nAstro,5); % animal
                AvsN_space(iSpace,2)=AstroData(nAstro,4); % spot
                AvsN_space(iSpace,3)=AstroData(nAstro,2); % trial
                AvsN_space(iSpace,4)=AstroData(nAstro,15); % unique astro ROI name
                AvsN_space(iSpace,5)=AstroData(nAstro,13); % astro ROI type
                AvsN_space(iSpace,6)=A_ROIarea;
                
                AvsN_space(iSpace,7)=NeuronalData(N_spaceMinIdx(iSpace),15); % unique neuro ROI name
                AvsN_space(iSpace,8)=NeuronalData(N_spaceMinIdx(iSpace),13); % neuro ROI type
                AvsN_space(iSpace,9)=N_ROIarea2; % astrocyte ROI area
                AvsN_space(iSpace,10) =distance(N_spaceMinIdx(iSpace),1);
                AvsN_space{iSpace,11}=AOnset; % astro onset time
                AvsN_space(iSpace,12)= NSpaceOnset(iSpace);% neuro onset time
                AvsN_space{iSpace,13}=AOnset-cell2mat(NSpaceOnset(iSpace)); %astrocytes-neurons
                AvsN_space{iSpace,14}=AstroData{nAstro,8}; %astro trace
                AvsN_space{iSpace,15}=NeuronalData{N_spaceMinIdx(iSpace),8}; %neuronal trace
            end
            
            % concatenate all data into a big matrix
            AvsN_TimeOnsets=vertcat(AvsN_TimeOnsets,AvsN_time);
            AvsN_SpaceOnsets=vertcat(AvsN_SpaceOnsets,AvsN_space);
            
            
            clear boundaries boundary1 boundary2 boundary1x boundary1y boundary2x boundary2y
            clear minDistance indexofMin distance AvsN_time AvsN_space NOnsets
        end
        
    end
end









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

%% histograms
figure('name', 'histogram: OnsetTime (A)- OnsetTime (N)- closest in time')
histogram(cell2mat(TimeComparisons(:,16)))



%% Raster Plot with event frequency histogram

% each peak time is a line
% lines are coloured and sorted based on ROI type

figure('name', 'Raster plot all ROIs')% ROIType{iType})
hold on
axis off
for iType=1:5
    ROIType={'Neuron','Neuropil','Endfeet','Dendrite','Process'};
    
    %find ROITypes
    for xROI= 1:size(responders,1)
        ROI_str(xROI)= strcmp(responders{xROI, 24},ROIType{iType});
    end
    
    % only unique ROIs of a given type
    uniqueROIs=unique(responders(ROI_str,24));
    
    colours={[0.6,0,0],[0,0,0.6],[0.2,0.4,0],[0.8,0.2,0],[0.8,0,0.8]};
    % plot raster plot of each ROI and trial
    
    %set(gca,'TickDir','out') % draw the tick marks on the outside
    %set(gca,'YTick', []) % don't draw y-axis ticks
    %set(gca,'PlotBoxAspectRatio',[1 0.05 1]) % short and wide
    %set(gca,'Color',get(gcf,'Color')) % match figure background
    %set(gca,'YColor',get(gcf,'Color')) % hide the y axis
    for iROI= 1:length(uniqueROIs)
        CurrentROI=uniqueROIs(iROI);
        for xROI=1:size(responders,1)
            ROI_idx(xROI)=strcmp(responders{xROI,25},CurrentROI);
        end
        peakTimes=cell2mat(responders(ROI_idx,6));
        if size(peakTimes,1) > size(peakTimes,2)
            peakTimes=peakTimes';
        end
        plot([peakTimes;peakTimes],[ones(size(peakTimes))+(iROI-1);zeros(size(peakTimes))+(iROI-1)],'Color',colours{1,iType},'LineWidth', 2)
    end
    %plot([0 20],[-1 -1], 'k--','LineWidth', 1)
    %axis([0 20 -1 2])
end

%% Responding Neurons and Astrocytes based on onset times

% ROI with a response to stimulation
% must have onset time during stimulation?

for iROI=1:length(Shortstim)
    rc_str(iROI)= ~isempty(strfind(Shortstim{iROI,3},'RCaMP'));
    gc_str(iROI)= ~isempty(strfind(Shortstim{iROI,3},'GCaMP'));
end

RCaMP=Shortstim(rc_str',:);
GCaMP=Shortstim(gc_str',:);

OT_RCaMP_traces=[];
for iROI=1:length(RCaMP)
    respOTIdx1(iROI)=~isempty(find(NOnsetWindow>RCaMP{iROI,16}>0));
    tempY= RCaMP{iROI,8};
    if length(tempY)>590
        tempY=tempY(1:nframes);
        if respOTIdx1(iROI)
            OT_RCaMP_traces= horzcat(OT_RCaMP_traces, tempY);
        end
    end
end
RrespOT=RCaMP(respOTIdx1',:); % responding neurons
OT_RCaMP_mean= mean(OT_RCaMP_traces');

OT_GCaMP_traces=[];
for iROI=1:length(GCaMP)
    respOTIdx2(iROI)=~isempty(find(GCaMP{iROI,16}>0 && GCaMP{iROI,16}<AOnsetWindow));
    tempY= GCaMP{iROI,8};
    if length(tempY)>590
        tempY=tempY(1:nframes);
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



%% ROITypes Shifted Traces- All Astrocyte traces shifted by neuronal onset times

figure ('name', 'Shifted traces- RCaMP vs GCaMP by Neuron onset time')
hold on
axis off
for xROI= 1:size(ShiftedTraces,1)
    tempY1 = ShiftedTraces{xROI,6};
    tempY2 = ShiftedTraces{xROI,9};
    %if length(tempY1)>590 && length(tempY2)>590
    tempY1=tempY1(1:nframes);
    tempY2=tempY2(1:nframes);
    tempX = ShiftedTraces{xROI,10}; % shifted by peak time
    grey = [0.8,0.8,0.8];
    plot(tempX,tempY1','Color',grey,'LineWidth',0.01);
    plot(tempX,(tempY2'+50),'Color',grey,'LineWidth',0.01);
    
    %end
end
plot([0 0],[-1 75], 'k--','LineWidth', 0.5)
%plot([0 8],[-2 -2], 'k','LineWidth', 2)

x1 = -20;
y1 = 5;
y2 = 50;
txt1 = 'RCaMP';
txt2 = 'GCaMP';
text(x1,y1,txt1,'HorizontalAlignment','right')
text(x1,y2,txt2,'HorizontalAlignment','right')

%% ROITypes Shifted Traces- All neuronal traces shifted by astrocyte peak times
stimwindow2=round(FrameRate*30);

figure ('name', 'Shifted traces- RCaMP vs GCaMP by astrocyte onset time')
hold on
axis off
for xROI= 1:size(ShiftedTraces,1)
    tempY1 = ShiftedTraces{xROI,6};
    tempY2 = ShiftedTraces{xROI,9};
    %if length(tempY1)>590 && length(tempY2)>590
    tempY1=tempY1(1:nframes);
    tempY2=tempY2(1:nframes);
    tempX = ShiftedTraces{xROI,11}; % shifted by AC peak time
    grey = [0.8,0.8,0.8];
    plot(tempX,tempY1','Color',grey,'LineWidth',0.01);
    plot(tempX,(tempY2'+50),'Color',grey,'LineWidth',0.01);
    
    %end
end
plot([0 0],[-1 75], 'k--','LineWidth', 0.5)
%plot([0 8],[-2 -2], 'k','LineWidth', 2)

x1 = -20;
y1 = 5;
y2 = 50;
txt1 = 'RCaMP';
txt2 = 'GCaMP';
text(x1,y1,txt1,'HorizontalAlignment','right')
text(x1,y2,txt2,'HorizontalAlignment','right')


%% Astrocyte Onsets Relative to Neuronal Events
% See Poskanzer and Yuste 2016 PNAS Fig. 2

% x axis is time
% y axis is trial number




%% responding ROIs based on peak times (max)

% ROI with a response to stimulation
% must have peak time around stimulation?

for iROI=1:length(ShortstimPeaks)
    rc_str2(iROI)= ~isempty(strfind(ShortstimPeaks{iROI,14},'RCaMP'));
    gc_str2(iROI)= ~isempty(strfind(ShortstimPeaks{iROI,14},'GCaMP'));
end

RCaMP_Peaks=ShortstimPeaks(rc_str2',:);
GCaMP_Peaks=ShortstimPeaks(gc_str2',:);

% NOTE: peak times have NOT be corrected for the baseline

for iROI=1:length(RCaMP_Peaks)
    respPTIdx1(iROI)=~isempty(find(RCaMP_Peaks{iROI,5}>5 && RCaMP_Peaks{iROI,5}<(NPTWindow+5)));%);
end
RrespPT=RCaMP_Peaks(respPTIdx1',:); % responding neurons

for iROI=1:length(GCaMP_Peaks)
    respPTIdx2(iROI)=~isempty(find(GCaMP_Peaks{iROI,5}>5 && GCaMP_Peaks{iROI,5}<(APTWindow+5)));%);
end
GrespPT=GCaMP_Peaks(respPTIdx2',:); % responding astrocytes

%% RCaMP vs GCaMP plots of responding ROIs based on peak time

PT_RCaMP_traces=[];
% exact traces for ROIs responding based on peak time
for xROI=1:length(RrespPT)
    CurrentROI=RrespPT(xROI,23);
    for iROI=1:length(Shortstim)
        r_roi=strfind(Shortstim{iROI,15},CurrentROI);
        if ~isempty(r_roi)
            tempY= Shortstim{iROI,8};
            if length(tempY)>590
                tempY=tempY(1:nframes);
                PT_RCaMP_traces= horzcat(PT_RCaMP_traces, tempY);
            end
        end
    end
end
PT_RCaMP_mean=mean(PT_RCaMP_traces');

PT_GCaMP_traces=[];
% exact traces for ROIs responding based on peak time
for xROI=1:length(GrespPT)
    CurrentROI=GrespPT(xROI,23);
    for iROI=1:length(Shortstim)
        r_roi=strfind(Shortstim{iROI,15},CurrentROI);
        if ~isempty(r_roi)
            tempY= Shortstim{iROI,8};
            if length(tempY)>590
                tempY=tempY(1:nframes);
                PT_GCaMP_traces= horzcat(PT_GCaMP_traces, tempY);
            end
        end
    end
end
PT_GCaMP_mean=mean(PT_GCaMP_traces');

% mean trace of astrocytes
figure ('name', 'Overlaid traces: All RCaMP responding PT ROIs')
hold on
axis off
for xROI= 1:size(PT_RCaMP_traces,2)
    tempY = PT_RCaMP_traces(:,xROI);
    if length(tempY)>590
        tempY=tempY(1:nframes);
        grey = [0.8,0.8,0.8];
        plot(TimeX,tempY,'Color',grey,'LineWidth',0.01);
    end
end
plot(TimeX, PT_RCaMP_mean, 'Color', 'k','LineWidth',1);

figure ('name', 'Overlaid traces: All GCaMP responding PT ROIs')
hold on
axis off
for xROI= 1:size(PT_GCaMP_traces,2)
    tempY = PT_GCaMP_traces(:,xROI);
    if length(tempY)>590
        tempY=tempY(1:nframes);
        grey = [0.8,0.8,0.8];
        plot(TimeX,tempY,'Color',grey,'LineWidth',0.01);
    end
end
plot(TimeX, PT_GCaMP_mean, 'Color', 'k','LineWidth',1);

%% how similar are responding groups?

RrespondingROIs=intersect(RrespOT(:,15), RrespPT(:,23)); %

GrespondingROIs=intersect(GrespOT(:,15), GrespPT(:,23)); %more similar




%% fast ROIs defined by onset time or AUC?

figure('name','1 sec auc vs onset time, astrocytes and neurons with stim onset times');
hold on
scatter(cell2mat(RrespOT(:,16)), cell2mat(RrespOT(:,17)));
scatter(cell2mat(GrespOT(:,16)), cell2mat(GrespOT(:,17)));

figure('name','1 sec auc vs onset time');
hold on
scatter(cell2mat(Shortstim(:,16)), cell2mat(Shortstim(:,17)));

% figure('name','20 sec auc vs onset time');
% hold on
% scatter(cell2mat(Shortstim(:,16)), cell2mat(Shortstim(:,18)));

for iROI=1:length(GrespOT)
    fastIdx(iROI)=~isempty(find(GrespOT{iROI,16}>0 && GrespOT{iROI,16}<=1));
    slowIdx(iROI)=~isempty(find(GrespOT{iROI,16}>1));
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
plot(TimeX, OT_RCaMP_mean, 'Color', 'r','LineWidth',1);
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
    fastIdx2(iROI)=~isempty(find(RrespOT{iROI,16}>0 && RrespOT{iROI,16}<=1));
end
fastN=RrespOT(fastIdx2',:);
MeanNeuronOnset=mean(cell2mat(fastN(:,16)));


%%
figure ('name', 'fast AC with onset times in first 1 sec-long stim')
hold on
axis off
xlim([0 35]);
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

%%
figure('name', 'Lck long stim all means')
hold on
axis off
xlim([0 35]);
ylim([-0.5 2.5]);
plot(TimeX, fastAC_mean, 'Color', 'k','LineWidth',1);
rectangle('Position', [5 -0.3 8 2.5])
plot(TimeX, OT_RCaMP_mean, 'Color', 'r','LineWidth',1);
plot(TimeX, slowAC_mean, 'Color', 'b','LineWidth',1);
plot([-1 -1],[0 1], 'k','LineWidth', 1)

%% shaded error bar with
green=[(27/255) (120/255) (55/255)];
purple=[(123/255) (50/255) (148/255)];
blue= [(0/255) (114/255) (178/255)];

% SEM calculations
fastAC_SDTrace = std(fastAC_traces');
fastAC_SEM=fastAC_SDTrace/sqrt(size(fastAC_traces,2));

slowAC_SDTrace = std(slowAC_traces');
slowAC_SEM=slowAC_SDTrace/sqrt(size(slowAC_traces,2));

RC_SDTrace = std(OT_RCaMP_traces');
RC_SEM=RC_SDTrace/sqrt(size(OT_RCaMP_traces,2));

figure('name', 'Lck short stim all means- plus SEM')
hold on
axis off
xlim([-1 35]);
lineProps.width = 1;
lineProps.edgestyle = ':';

%ylim([-0.2 3]);
lineProps.col = {blue};
mseb(TimeX,slowAC_mean,slowAC_SEM,lineProps)
lineProps.col = {green};
mseb(TimeX,(fastAC_mean+0.75),fastAC_SEM,lineProps)
lineProps.col = {purple};
mseb(TimeX,(OT_RCaMP_mean+2),RC_SEM,lineProps)

% H2=shadedErrorBar(TimeX,(fastAC_mean+0.75),fastAC_SEM)%,green,1);
% H3=shadedErrorBar(TimeX,(OT_RCaMP_mean+2),RC_SEM)%,purple,1);
rectangle('Position', [5 -0.3 8 4])
plot([-0.1 -0.1],[0 1], 'k','LineWidth', 1)


figure('name', 'Lck short stim all means- plus SD')
hold on
axis off
xlim([-1 35]);
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
%%


