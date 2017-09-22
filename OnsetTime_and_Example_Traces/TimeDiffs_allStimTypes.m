
clearvars
%close all

%% Load data

% time windows based on stimulation
NOnsetWindow= 1; %1 for short stim % neuronal onset times
AOnsetWindow= 15; % 5 for short stim % astrocyte onset times
Fast_AOnsetWindow=1;
NPTWindow= 2; % one second longer than stimulation for peak times
APTWindow= 10; % astrocyte longer than stimulation for peak times

% save files names
saveFiles1='D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\shortstim_firstonset_comparisons_fixedDis.mat';
saveFiles2='D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\shortstim_firstonset_comparisons_fixedDis.csv';
saveFiles3= 'D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\shortstim_onset&AUC.csv';

%peak data
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_2D_nostim_05_04_2017.mat');
load('D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_2D_shortstim_28_04_2017.mat');

% Load trace data
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_nostim_05_04_2017.mat');
load('D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_shortstim_28_04_2017.mat');
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_Shortstim_05_04_2017.mat');
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_Shortstim_05_04_2017.mat');


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
        
        
        % find similar astrocyte peak times
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
                            Image=zeros(127,128);
                            Image(NeuronalData{nNeuro,10})=1;
                            ExtraRow=zeros(1,128);
                            Image1=[Image;ExtraRow];
                            Image1=im2bw(Image1);
                        elseif islogical(NeuronalData{nNeuro,10})
                            Image1= double(NeuronalData{nNeuro,10});
                        else
                            Image1=[];
                        end
                        
                        % for the second ROI
                        if isnumeric(AstroData{nAstro,10})
                            Image=zeros(127,128);
                            Image(AstroData{nAstro,10})=1;
                            ExtraRow=zeros(1,128);
                            Image2=[Image;ExtraRow];
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
                                OnsetTimeComparisons{iA,16}=TimeDiff;
                                OnsetTimeComparisons{iA,17}=strcat(OnsetTimeComparisons{iA,7},OnsetTimeComparisons{iA,14}); % astrocyte peak name
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










clearvars
%close all

%% Load data

% time windows based on stimulation
NOnsetWindow= 1; %1 for short stim % neuronal onset times
AOnsetWindow= 15; % 5 for short stim % astrocyte onset times
Fast_AOnsetWindow=1;
NPTWindow= 2; % one second longer than stimulation for peak times
APTWindow= 10; % astrocyte longer than stimulation for peak times

% save files names
saveFiles1='D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\longstim_firstonset_comparisons_fixedDis.mat';
saveFiles2='D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\longstim_firstonset_comparisons_fixedDis.csv';
saveFiles3= 'D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\longstim_onset&AUC.csv';

%peak data
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_2D_nostim_05_04_2017.mat');
load('D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_2D_longstim_28_04_2017.mat');

% Load trace data
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_nostim_05_04_2017.mat');
load('D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_longstim_28_04_2017.mat');
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_Shortstim_05_04_2017.mat');
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_Shortstim_05_04_2017.mat');


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
        
        
        % find similar astrocyte peak times
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
                            Image=zeros(127,128);
                            Image(NeuronalData{nNeuro,10})=1;
                            ExtraRow=zeros(1,128);
                            Image1=[Image;ExtraRow];
                            Image1=im2bw(Image1);
                        elseif islogical(NeuronalData{nNeuro,10})
                            Image1= double(NeuronalData{nNeuro,10});
                        else
                            Image1=[];
                        end
                        
                        % for the second ROI
                        if isnumeric(AstroData{nAstro,10})
                            Image=zeros(127,128);
                            Image(AstroData{nAstro,10})=1;
                            ExtraRow=zeros(1,128);
                            Image2=[Image;ExtraRow];
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
                                OnsetTimeComparisons{iA,16}=TimeDiff;
                                OnsetTimeComparisons{iA,17}=strcat(OnsetTimeComparisons{iA,7},OnsetTimeComparisons{iA,14}); % astrocyte peak name
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






clearvars
%close all

%% Load data

% time windows based on stimulation
NOnsetWindow= 1; %1 for short stim % neuronal onset times
AOnsetWindow= 15; % 5 for short stim % astrocyte onset times
Fast_AOnsetWindow=1;
NPTWindow= 2; % one second longer than stimulation for peak times
APTWindow= 10; % astrocyte longer than stimulation for peak times

% save files names
saveFiles1='D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\nostim_firstonset_comparisons_fixedDis.mat';
saveFiles2='D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\nostim_firstonset_comparisons_fixedDis.csv';
saveFiles3= 'D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\nostim_onset&AUC.csv';

%peak data
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_2D_nostim_05_04_2017.mat');
load('D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_2D_nostim_28_04_2017.mat');

% Load trace data
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_nostim_05_04_2017.mat');
load('D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_nostim_28_04_2017.mat');
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_Shortstim_05_04_2017.mat');
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_Shortstim_05_04_2017.mat');


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
        
        
        % find similar astrocyte peak times
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
                            Image=zeros(127,128);
                            Image(NeuronalData{nNeuro,10})=1;
                            ExtraRow=zeros(1,128);
                            Image1=[Image;ExtraRow];
                            Image1=im2bw(Image1);
                        elseif islogical(NeuronalData{nNeuro,10})
                            Image1= double(NeuronalData{nNeuro,10});
                        else
                            Image1=[];
                        end
                        
                        % for the second ROI
                        if isnumeric(AstroData{nAstro,10})
                            Image=zeros(127,128);
                            Image(AstroData{nAstro,10})=1;
                            ExtraRow=zeros(1,128);
                            Image2=[Image;ExtraRow];
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
                                OnsetTimeComparisons{iA,16}=TimeDiff;
                                OnsetTimeComparisons{iA,17}=strcat(OnsetTimeComparisons{iA,7},OnsetTimeComparisons{iA,14}); % astrocyte peak name
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





clearvars
%close all

%% Load data

% time windows based on stimulation
NOnsetWindow= 1; %1 for short stim % neuronal onset times
AOnsetWindow= 15; % 5 for short stim % astrocyte onset times
Fast_AOnsetWindow=1;
NPTWindow= 2; % one second longer than stimulation for peak times
APTWindow= 10; % astrocyte longer than stimulation for peak times

% save files names
saveFiles1='D:\Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\shortstim_firstonset_comparisons_fixedDis.mat';
saveFiles2='D:\Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\shortstim_firstonset_comparisons_fixedDis.csv';
saveFiles3= 'D:\Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\shortstim_onset&AUC_fixedDis.csv';

%peak data
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_2D_nostim_05_04_2017.mat');
load('D:\Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\cytoGC&RC_2D_shortstim_28_04_2017.mat');

% Load trace data
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_nostim_05_04_2017.mat');
load('D:\Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\cytoGC&RC_traces_2D_shortstim_28_04_2017.mat');
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_Shortstim_05_04_2017.mat');
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_Shortstim_05_04_2017.mat');


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
        
        
        % find similar astrocyte peak times
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
                            Image=zeros(127,128);
                            Image(NeuronalData{nNeuro,10})=1;
                            ExtraRow=zeros(1,128);
                            Image1=[Image;ExtraRow];
                            Image1=im2bw(Image1);
                        elseif islogical(NeuronalData{nNeuro,10})
                            Image1= double(NeuronalData{nNeuro,10});
                        else
                            Image1=[];
                        end
                        
                        % for the second ROI
                        if isnumeric(AstroData{nAstro,10})
                            Image=zeros(127,128);
                            Image(AstroData{nAstro,10})=1;
                            ExtraRow=zeros(1,128);
                            Image2=[Image;ExtraRow];
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
                                OnsetTimeComparisons{iA,16}=TimeDiff;
                                OnsetTimeComparisons{iA,17}=strcat(OnsetTimeComparisons{iA,7},OnsetTimeComparisons{iA,14}); % astrocyte peak name
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






clearvars
%close all

%% Load data

% time windows based on stimulation
NOnsetWindow= 1; %1 for short stim % neuronal onset times
AOnsetWindow= 15; % 5 for short stim % astrocyte onset times
Fast_AOnsetWindow=1;
NPTWindow= 2; % one second longer than stimulation for peak times
APTWindow= 10; % astrocyte longer than stimulation for peak times

% save files names
saveFiles1='D:\Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\longstim_firstonset_comparisons_fixedDis.mat';
saveFiles2='D:\Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\longstim_firstonset_comparisons_fixedDis.csv';
saveFiles3= 'D:\Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\longstim_onset&AUC_fixedDis.csv';

%peak data
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_2D_nostim_05_04_2017.mat');
load('D:\Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\cytoGC&RC_2D_longstim_28_04_2017.mat');

% Load trace data
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_nostim_05_04_2017.mat');
load('D:\Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\cytoGC&RC_traces_2D_longstim_28_04_2017.mat');
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_Shortstim_05_04_2017.mat');
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_Shortstim_05_04_2017.mat');


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
        
        
        % find similar astrocyte peak times
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
                            Image=zeros(127,128);
                            Image(NeuronalData{nNeuro,10})=1;
                            ExtraRow=zeros(1,128);
                            Image1=[Image;ExtraRow];
                            Image1=im2bw(Image1);
                        elseif islogical(NeuronalData{nNeuro,10})
                            Image1= double(NeuronalData{nNeuro,10});
                        else
                            Image1=[];
                        end
                        
                        % for the second ROI
                        if isnumeric(AstroData{nAstro,10})
                            Image=zeros(127,128);
                            Image(AstroData{nAstro,10})=1;
                            ExtraRow=zeros(1,128);
                            Image2=[Image;ExtraRow];
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
                                OnsetTimeComparisons{iA,16}=TimeDiff;
                                OnsetTimeComparisons{iA,17}=strcat(OnsetTimeComparisons{iA,7},OnsetTimeComparisons{iA,14}); % astrocyte peak name
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






clearvars
%close all

%% Load data

% time windows based on stimulation
NOnsetWindow= 1; %1 for short stim % neuronal onset times
AOnsetWindow= 15; % 5 for short stim % astrocyte onset times
Fast_AOnsetWindow=1;
NPTWindow= 2; % one second longer than stimulation for peak times
APTWindow= 10; % astrocyte longer than stimulation for peak times

% save files names
saveFiles1='D:\Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\nostim_firstonset_comparisons_fixedDis.mat';
saveFiles2='D:\Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\nostim_firstonset_comparisons_fixedDis.csv';
saveFiles3= 'D:\Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\nostim_onset&AUC_fixedDis.csv';

%peak data
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_2D_nostim_05_04_2017.mat');
load('D:\Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\cytoGC&RC_2D_nostim_28_04_2017.mat');

% Load trace data
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_nostim_05_04_2017.mat');
load('D:\Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\cytoGC&RC_traces_2D_nostim_28_04_2017.mat');
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_Shortstim_05_04_2017.mat');
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_Shortstim_05_04_2017.mat');


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
        
        
        % find similar astrocyte peak times
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
                            Image=zeros(127,128);
                            Image(NeuronalData{nNeuro,10})=1;
                            ExtraRow=zeros(1,128);
                            Image1=[Image;ExtraRow];
                            Image1=im2bw(Image1);
                        elseif islogical(NeuronalData{nNeuro,10})
                            Image1= double(NeuronalData{nNeuro,10});
                        else
                            Image1=[];
                        end
                        
                        % for the second ROI
                        if isnumeric(AstroData{nAstro,10})
                            Image=zeros(127,128);
                            Image(AstroData{nAstro,10})=1;
                            ExtraRow=zeros(1,128);
                            Image2=[Image;ExtraRow];
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
                                OnsetTimeComparisons{iA,16}=TimeDiff;
                                OnsetTimeComparisons{iA,17}=strcat(OnsetTimeComparisons{iA,7},OnsetTimeComparisons{iA,14}); % astrocyte peak name
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


