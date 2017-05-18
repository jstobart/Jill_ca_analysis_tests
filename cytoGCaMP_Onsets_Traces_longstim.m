
clearvars
close all

%% Load data

% time windows based on stimulation
NOnsetWindow= 30; %1 for short stim % neuronal onset times
AOnsetWindow= 30; % 5 for short stim % astrocyte onset times
Fast_AOnsetWindow=1;
Fast_NOnsetWindow=1;
NPTWindow= 9; % one second longer than stimulation for peak times
APTWindow= 12; % astrocyte longer than stimulation for peak times

% save files names
saveFiles1='D:\Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\cyto_longstim_firstonset_comparisons.mat';
saveFiles2='D:\Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\cyto_longstim_firstonset_comparisons.csv';
saveFiles3= 'D:\Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\cyto_longstim_onset&AUC.csv';

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

ACIdx= find(~cellfun(@isempty, regexp(Shortstim(:,3), 'GCaMP')));
figure('name', 'Astrocyte onset times')
histogram(cell2mat(Shortstim(ACIdx,16)), 'binwidth',0.8450)

NIdx= find(~cellfun(@isempty, regexp(Shortstim(:,3), 'RCaMP')));
figure('name', 'Neuronal onset times')
histogram(cell2mat(Shortstim(NIdx,16)), 'binwidth', 0.845)

figure('name', 'Astrocyte AUC')
histogram(cell2mat(Shortstim(ACIdx,17)))

figure('name', 'Neuronal AUC')
histogram(cell2mat(Shortstim(NIdx,17)))

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
    fastIdx2(iROI)=~isempty(find(RrespOT{iROI,16}>0 && RrespOT{iROI,16}<Fast_NOnsetWindow));
end
fastN=RrespOT(fastIdx2',:);
MeanNeuronOnset=mean(cell2mat(fastN(:,16)));


% y scaled
stimwindow=round(30*11.84);
figure ('name', 'Mean traces only: RCaMP, fast AC, slow AC cyto GCaMP6s')
hold on
axis off
plot(TimeX, fastAC_mean, 'Color', 'k','LineWidth',1);
plot(TimeX, OT_RCaMP_mean, 'Color', 'r','LineWidth',1);
plot(TimeX, slowAC_mean, 'Color', 'b','LineWidth',1);
plot([5 5],[-1 2], 'k--','LineWidth', 1)
plot([5 6],[0 0], 'k','LineWidth', 3)
legend('fast cyto','RCaMP','slow cyto')


% IdxAUC= find(Shortstim(:,17)>5);
% fastAUC=Shortstim(IdxAUC,:);



% label ROIs as responding or not 










%% Group Responders by ROIType
% Resp_RCaMP_traces=[];
% Resp_GCaMP_traces=[];
% Resp_somata_traces=[];
% Resp_EF_traces=[];
% Resp_processes_traces=[];
% Resp_neuron_traces=[];
% Resp_neuropil_traces=[];
% 
% %find GCaMP, RCaMP
% for xROI= 1:length(RespondingROIs)
%     RC_str= strfind(RespondingROIs{xROI, 3},'RCaMP');
%     GC_str= strfind(RespondingROIs{xROI, 3},'GCaMP');
%     
%     tempY= RespondingROIs{xROI,8};
%     if length(tempY)>590
%         tempY=tempY(1:nframes);
%         if ~isempty(RC_str)
%             Resp_RCaMP_traces= horzcat(Resp_RCaMP_traces, tempY);
%         elseif ~isempty(GC_str)
%             Resp_GCaMP_traces= horzcat(Resp_GCaMP_traces, tempY);
%         end
%     end
% end
% 
% 
% %find ROITypes
% for xROI= 1:length(RespondingROIs)
%     N_str= strfind(RespondingROIs{xROI, 13},'Neuron');
%     S_str= strfind(RespondingROIs{xROI, 13},'Dendrite');
%     EF_str= strfind(RespondingROIs{xROI, 13},'Endfeet');
%     P_str= strfind(RespondingROIs{xROI, 13},'Process');
%     NP_str= strfind(RespondingROIs{xROI, 13},'Neuropil');
%     
%     tempY= RespondingROIs{xROI,8};
%     if length(tempY)>590
%         tempY=tempY(1:nframes);
%         if ~isempty(N_str)
%             Resp_neuron_traces= horzcat(Resp_neuron_traces, tempY);
%         elseif ~isempty(S_str)
%             Resp_somata_traces= horzcat(Resp_somata_traces, tempY);
%         elseif ~isempty(EF_str)
%             Resp_EF_traces= horzcat(Resp_EF_traces, tempY);
%         elseif ~isempty(P_str)
%             Resp_processes_traces= horzcat(Resp_processes_traces, tempY);
%         elseif ~isempty(NP_str)
%             Resp_neuropil_traces= horzcat(Resp_neuropil_traces, tempY);
%         end
%     end
% end
% 
% % means and SDs
% %RCaMP
% RespRCmeanTrace = mean(Resp_RCaMP_traces,2);
% RespRCSDTrace = std(Resp_RCaMP_traces');
% 
% %GcaMP
% RespGCmeanTrace = mean(Resp_GCaMP_traces,2);
% RespGCSDTrace = std(Resp_GCaMP_traces');
% 
% %neurons
% RespNmeanTrace = mean(Resp_neuron_traces,2);
% RespNSDTrace = std(Resp_neuron_traces');
% 
% %neuropil
% RespNPmeanTrace = mean(Resp_neuropil_traces,2);
% RespNPSDTrace = std(Resp_neuropil_traces');
% 
% %somata
% RespSmeanTrace = mean(Resp_somata_traces,2);
% RespSSDTrace = std(Resp_somata_traces');
% 
% %endfeet
% RespEFmeanTrace = mean(Resp_EF_traces,2);
% RespEFSDTrace = std(Resp_EF_traces');
% 
% %processes
% RespPmeanTrace = mean(Resp_processes_traces,2);
% RespPSDTrace = std(Resp_processes_traces');
% 
% 
% %% Mean RCaMP vs GCaMP
% figure ('name', 'Mean traces only: GCaMP vs RCaMP responding ROIs')
% hold on
% axis off
% plot(TimeX(1:stimwindow),smooth(RespRCmeanTrace(1:stimwindow)',5),'r','LineWidth',1.5);
% plot(TimeX(1:stimwindow),smooth(RespGCmeanTrace(1:stimwindow)',5),'g','LineWidth', 1.5);
% plot([5 5],[-1 2], 'k--','LineWidth', 1)
% plot([5 13],[0 0], 'k','LineWidth', 3)
% legend('RCaMP','GCaMP')
% 
% 
% 
% %% ROITypes
% figure ('name', 'Overlaid traces + mean:All responding ROItype')
% subplot(1,5,1)
% hold on
% axis off
% ylim([-5 40]);
% for xROI= 1:size(Resp_neuron_traces,2)
%     tempY = Resp_neuron_traces(:,xROI);
%     if length(tempY)>600
%         tempY=tempY(1:nframes);
%     end
%     grey = [0.8,0.8,0.8];
%     plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.01);
%     
% end
% plot([5 13],[-2 -2], 'k','LineWidth', 2)
% plot(TimeX(1:stimwindow), RespNmeanTrace(1:stimwindow), 'k', 'LineWidth',1)
% plot([5 5],[-1 10], 'k--','LineWidth', 0.5)
% x1 = 0;
% y1 = 5;
% txt1 = 'Neurons';
% text(x1,y1,txt1,'HorizontalAlignment','right')
% 
% subplot(1,5,2)
% hold on
% axis off
% ylim([-5 40]);
% for xROI= 1:size(Resp_neuropil_traces,2)
%     tempY = Resp_neuropil_traces(:,xROI);
%     if length(tempY)>600
%         tempY=tempY(1:nframes);
%     end
%     grey = [0.8,0.8,0.8];
%     plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.01);
% end
% plot([5 13],[-2 -2], 'k','LineWidth', 2)
% plot(TimeX(1:stimwindow), RespNPmeanTrace(1:stimwindow), 'k', 'LineWidth',1)
% plot([5 5],[-1 10], 'k--','LineWidth', 0.5)
% x1 = 0;
% y1 = 5;
% txt1 = 'Neuropil';
% text(x1,y1,txt1,'HorizontalAlignment','right')
% 
% 
% subplot(1,5,4)
% hold on
% axis off
% ylim([-5 40]);
% for xROI= 1:size(Resp_EF_traces,2)
%     tempY = Resp_EF_traces(:,xROI);
%     if length(tempY)>600
%         tempY=tempY(1:nframes);
%     end
%     grey = [0.8,0.8,0.8];
%     plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.01);
% end
% plot([5 13],[-2 -2], 'k','LineWidth', 2)
% plot(TimeX(1:stimwindow), RespEFmeanTrace(1:stimwindow), 'k', 'LineWidth',1)
% plot([5 5],[-1 10], 'k--','LineWidth', 0.5)
% x1 = 0;
% y1 = 5;
% txt1 = 'Endfeet';
% text(x1,y1,txt1,'HorizontalAlignment','right')
% 
% 
% subplot(1,5,3)
% hold on
% axis off
% ylim([-5 40]);
% for xROI= 1:size(Resp_somata_traces,2)
%     tempY = Resp_somata_traces(:,xROI);
%     if length(tempY)>600
%         tempY=tempY(1:nframes);
%     end
%     grey = [0.8,0.8,0.8];
%     plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.01);
% end
% plot([5 13],[-2 -2], 'k','LineWidth', 2)
% plot(TimeX(1:stimwindow), RespSmeanTrace(1:stimwindow), 'k', 'LineWidth',1)
% plot([5 5],[-1 10], 'k--','LineWidth', 0.5)
% x1 = 0;
% y1 = 5;
% txt1 = 'Dendrites';
% text(x1,y1,txt1,'HorizontalAlignment','right')
% 
% 
% subplot(1,5,5)
% hold on
% axis off
% ylim([-5 40]);
% for xROI= 1:size(Resp_processes_traces,2)
%     tempY = Resp_processes_traces(:,xROI);
%     if length(tempY)>600
%         tempY=tempY(1:nframes);
%     end
%     grey = [0.8,0.8,0.8];
%     plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.01);
% end
% plot([5 13],[-2 -2], 'k','LineWidth', 2)
% plot(TimeX(1:stimwindow), RespPmeanTrace(1:stimwindow), 'k', 'LineWidth',1)
% plot([5 5],[-1 10], 'k--','LineWidth', 0.5)
% x1 = 0;
% y1 = 5;
% txt1 = 'Processes';
% text(x1,y1,txt1,'HorizontalAlignment','right')
% 
% 


