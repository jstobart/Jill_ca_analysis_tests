
clearvars
close all

%% Load data

%peak data
load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_2D_nostim_05_04_2017.mat');
NostimPeaks=AllData2(2:end,:);

% Load trace data
load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_nostim_05_04_2017.mat');
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_shortstim_05_04_2017.mat');
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_longstim_05_04_2017.mat');

Nostim=All_traces;
%%
% get info for plotting
FrameRate=11.84;
nframes=592;
TimeX(1:nframes) = (1:nframes)/FrameRate;

% get rid of handclicked neuropil ROIs because they are not relevant if
% FLIKA ROIs are included for the neuron channel

% get rid of astrocyte neuropil traces
for xROI=1:size(Nostim,1)
    GCNP2(xROI)= strcmp(Nostim{xROI,1},'np');
end
Nostim=Nostim(~GCNP2',:);

for iROI=1:length(Nostim)
    %find ROITypes
    N_str= strfind(Nostim{iROI, 1},'N');
    D_str= strfind(Nostim{iROI, 1},'D'); %hand selected dendrite
    r_str= strfind(Nostim{iROI, 1},'r'); %FLIKA ROIs
    EF_str= strfind(Nostim{iROI, 1},'E');
    
    if ~isempty(N_str)
        Nostim{iROI,13}='Neuron';
    elseif ~isempty(D_str)
        Nostim{iROI,13}='Dendrite';
    elseif ~isempty(r_str)
        P_str= strfind(Nostim{iROI, 3},'GCaMP');
        if ~isempty(P_str)
            Nostim{iROI,13}='Process';
        else
            Nostim{iROI,13}='Dendrite';
        end
    elseif ~isempty(EF_str)
        Nostim{iROI,13}='Endfeet';
    end
    
    % make new unique trial names
    Nostim{iROI,14}=strcat(Nostim{iROI,5},'_',Nostim{iROI,4},'_',Nostim{iROI,2});
    % unique ROI names
    Nostim{iROI,15}=strcat(Nostim{iROI,5},'_',Nostim{iROI,4},'_',Nostim{iROI,2}, '_',Nostim{iROI,1});
end

for iROI=1:length(Nostim)
    %if a process ROI is overlapping with a soma or process, exclude it
    %from data
    nonOverlapIdx(iROI)=~ischar(Nostim{iROI,12});
end

% remove overlapping processes
Nostim = Nostim(nonOverlapIdx',:);

%%

for iROI=1:length(NostimPeaks)
    %find ROITypes
    N_str= strfind(NostimPeaks{iROI, 1},'N');
    D_str= strfind(NostimPeaks{iROI, 1},'D'); %hand selected dendrite
    r_str= strfind(NostimPeaks{iROI, 1},'r'); %FLIKA ROIs
    EF_str= strfind(NostimPeaks{iROI, 1},'E');
    
    if ~isempty(N_str)
        NostimPeaks{iROI,21}='Neuron';
    elseif ~isempty(D_str)
        NostimPeaks{iROI,21}='Dendrite';
    elseif ~isempty(r_str)
        P_str= strfind(NostimPeaks{iROI, 3},'GCaMP');
        if ~isempty(P_str)
            NostimPeaks{iROI,21}='Process';
        else
            NostimPeaks{iROI,21}='Dendrite';
        end
    elseif ~isempty(EF_str)
        NostimPeaks{iROI,21}='Endfeet';
    end
    
    % make new unique trial names
    NostimPeaks{iROI,22}=strcat(NostimPeaks{iROI,13},'_',NostimPeaks{iROI,15},'_',NostimPeaks{iROI,12});
    % unique ROI names
    NostimPeaks{iROI,23}=strcat(NostimPeaks{iROI,13},'_',NostimPeaks{iROI,15},'_',NostimPeaks{iROI,12}, '_',NostimPeaks{iROI,10});
end

for iROI=1:length(NostimPeaks)
    %if a process ROI is overlapping with a soma or process, exclude it
    %from data
    nonOverlapIdx(iROI)=~ischar(NostimPeaks{iROI,19});
end

% remove overlapping processes
NostimPeaks = NostimPeaks(nonOverlapIdx',:);


%% Find astrocyte peaks that occur around spontaneous neuronal peaks
baselineCorrectedTime=TimeX-5;

% % % find peak onsets for each ROI
% for iROI= 1:length(Nostim)
%     trace=Nostim{iROI,8};
%     Onsets=find_multiple_onset_times(baselineCorrectedTime(10:end), trace(10:592,:),2.5,1,1);
%     if isempty(Onsets)
%         Onsets=nan(1,1);
%     end
%     Nostim{iROI, 16}= Onsets;
% end

%% compare peak onsets for each RCaMP and GCaMP ROI

tic
TimeComparisons=[];
Trials= unique(Nostim(:,14));

% onset times comparisons
for itrial=1:length(Trials)
    CurrentTrial=Trials(itrial);
    
    % Find the idx of paths matching trial
    matchingTrialIdx = find(~cellfun(@isempty, regexp(Nostim(:,14), CurrentTrial)));
    TrialData = Nostim(matchingTrialIdx,:);
    % find neuronal onset times in this trial
    NeuronalIdx = find(~cellfun(@isempty, regexp(TrialData(:,3), 'RCaMP')));
    NeuronalData=TrialData(NeuronalIdx,:);
    
    
    % find similar astrocyte peak times
    AstroIdx = find(~cellfun(@isempty, regexp(TrialData(:,3), 'GCaMP')));
    AstroData=TrialData(AstroIdx,:);
    
    
    parfor nNeuro=1:size(NeuronalData,1)
        for nAstro= 1:size(AstroData,1)
            Ntrace=NeuronalData{nNeuro,8};
            NOnset=find_multiple_onset_times(baselineCorrectedTime(10:end), Ntrace(10:592,:),2.5,1);
            Atrace=AstroData{nAstro,8};
            AOnset=find_multiple_onset_times(baselineCorrectedTime(10:end), Atrace(10:592,:),2.5,1);
            
            if ~isempty(NOnset)
                if ~isempty(AOnset)
                    
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
                    
                    for iN=1:length(NOnset)
                        OnsetTimeComparisons=cell(length(AOnset),14);
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
                            OnsetTimeComparisons(iA,6)=AstroData(nAstro,15); % unique astrocyte ROI name
                            OnsetTimeComparisons(iA,7)=AstroData(nAstro,13); % astrocyte ROI type
                            OnsetTimeComparisons{iA,8} =distance;
                            OnsetTimeComparisons{iA,9} = numberOfBoundaries; % number of ROIs in the mask
                            
                            
                            % Onset times
                            OnsetTimeComparisons{iA,10}=strcat('NPeak',num2str(iN,'%02d')); % neuronal peak number
                            OnsetTimeComparisons{iA,11}=NeuronalOnset; % neuronal onset time
                            OnsetTimeComparisons{iA,12}=strcat('APeak',num2str(iA,'%02d')); % neuronal peak number
                            OnsetTimeComparisons{iA,13}=AstrocyteOnset; % astrocyte onset time
                            
                            % difference between astrocyte onset and neuronal onset
                            OnsetTimeComparisons{iA,14}=TimeDiff;
                            OnsetTimeComparisons{iA,15}=strcat(OnsetTimeComparisons{iA,6},OnsetTimeComparisons{iA,12}); % astrocyte peak name
                            %OnsetTimeComparisons{iA,15}=strcat(OnsetTimeComparisons{iA,6},OnsetTimeComparisons{iA,12}); % astrocyte peak name
                        end
                        
                        TimeComparisons=vertcat(TimeComparisons,OnsetTimeComparisons); % concatenate all data into a big matrix
                        
                        %clear boundaries boundary1 boundary2 boundary1x boundary1y boundary2x boundary2y minDistance indexofMin OnsetTimeComparisons
                    end
                end
            end
        end
    end
end



toc
% % histograms
figure('name', 'histogram: OnsetTime (A)- OnsetTime (N)')
histogram(cell2mat(TimeComparisons(:,14)))

%% Isolate the early traces from the shifted data

for xComp= 1:size(TimeComparisons,1)
    mixedIdx(xComp)= TimeComparisons{xComp,14}>=-5 && TimeComparisons{xComp,14}<=5;
end
mixedOnset=TimeComparisons(mixedIdx,:);

figure('name','astrocyte peaks around neurons');histogram(cell2mat(mixedOnset(:,14)),'BinWidth',0.0845);


%% percent distributions

% before, after, zero
for xComp=1:length(mixedOnset)
   zeroIdx(xComp)= mixedOnset{xComp,14}==0;
    
   beforeNIdx(xComp)= mixedOnset{xComp,14}>=-1 && TimeComparisons{xComp,14}<0;
   
   afterNIdx(xComp)= mixedOnset{xComp,14}>0 && TimeComparisons{xComp,14}<=1;
end
zerodata=mixedOnset(zeroIdx,:);
beforedata=mixedOnset(beforeNIdx,:);
afterdata=mixedOnset(afterNIdx,:);


% ALL Astrocyte Peaks
% percentage of ALL astrocyte peaks that come before, after, or zero
Total_APeak_Num=length(unique(TimeComparisons(:,15)));

NumZeroPeaks=length(unique(zerodata(:,15)));
% number of astrocyte peaks that appear in 1 sec before neurons
NumBeforePeaks=length(unique(beforedata(:,15)));
% number of astrocyte peaks that appear 1 sec after neurons
NumAfterPeaks=length(unique(afterdata(:,15)));

Percent_All_zero=(NumZeroPeaks/Total_APeak_Num)*100;
Percent_All_before=(NumBeforePeaks/Total_APeak_Num)*100;
Percent_All_after=(NumAfterPeaks/Total_APeak_Num)*100;


% ALL Unique neuron-astrocyte peak comparisons
% percentage of unique comparisons that come before, after or zero
Total_Comp_Num=length(TimeComparisons(:,15));

NumZeroComp=length(zerodata(:,15));
% number of astrocyte peaks that appear in 1 sec before neurons
NumBeforeComp=length(beforedata(:,15));
% number of astrocyte peaks that appear 1 sec after neurons
NumAfterComp=length(afterdata(:,15));

Percent_Comp_zero=(NumZeroComp/Total_Comp_Num)*100;
Percent_Comp_before=(NumBeforeComp/Total_Comp_Num)*100;
Percent_Comp_after=(NumAfterComp/Total_Comp_Num)*100;


%% percent astrocyte ROIs per field of view with at least one comparison

% do this in R


%% Distance between ROIs that are compared

figure('name', 'astrocyte onset time diff vs ROI distance');
onsets=cell2mat(mixedOnset(:,14));
distances=cell2mat(mixedOnset(:,8));
scatter(onsets, distances, 10, 'filled');


%% Area under the curve for each ROI with an onset time near the neuron?
% only area in second before or second after?
% based on astrocyte onset time



%% Plot all traces for each Condition and Channel

% individual traces in grey, mean in red

traces=[];
figure ('name', 'Overlaid traces-All ROIs- whole trial')
hold on
axis off
for xROI= 1:nROIs
    tempY = data_traces{xROI,8};
    if length(tempY)>590
        tempY=tempY(1:nframes);
        grey = [0.8,0.8,0.8];
        plot(TimeX,tempY,'Color',grey,'LineWidth',0.1);
        traces= horzcat(traces, tempY);
    end
end

plot([5 13],[-2 -2], 'k','LineWidth', 2)
meanTrace = mean(traces,2);
plot(TimeX, meanTrace, 'k', 'LineWidth',1)
% plot([0 0],[18.5 19.5], 'k','LineWidth', 2)





%% Group Responders by ROIType
Resp_RCaMP_traces=[];
Resp_GCaMP_traces=[];
Resp_somata_traces=[];
Resp_EF_traces=[];
Resp_processes_traces=[];
Resp_neuron_traces=[];
Resp_neuropil_traces=[];

%find GCaMP, RCaMP
for xROI= 1:length(RespondingROIs)
    RC_str= strfind(RespondingROIs{xROI, 3},'RCaMP');
    GC_str= strfind(RespondingROIs{xROI, 3},'GCaMP');
    
    tempY= RespondingROIs{xROI,8};
    if length(tempY)>590
        tempY=tempY(1:nframes);
        if ~isempty(RC_str)
            Resp_RCaMP_traces= horzcat(Resp_RCaMP_traces, tempY);
        elseif ~isempty(GC_str)
            Resp_GCaMP_traces= horzcat(Resp_GCaMP_traces, tempY);
        end
    end
end


%find ROITypes
for xROI= 1:length(RespondingROIs)
    N_str= strfind(RespondingROIs{xROI, 13},'Neuron');
    S_str= strfind(RespondingROIs{xROI, 13},'Dendrite');
    EF_str= strfind(RespondingROIs{xROI, 13},'Endfeet');
    P_str= strfind(RespondingROIs{xROI, 13},'Process');
    NP_str= strfind(RespondingROIs{xROI, 13},'Neuropil');
    
    tempY= RespondingROIs{xROI,8};
    if length(tempY)>590
        tempY=tempY(1:nframes);
        if ~isempty(N_str)
            Resp_neuron_traces= horzcat(Resp_neuron_traces, tempY);
        elseif ~isempty(S_str)
            Resp_somata_traces= horzcat(Resp_somata_traces, tempY);
        elseif ~isempty(EF_str)
            Resp_EF_traces= horzcat(Resp_EF_traces, tempY);
        elseif ~isempty(P_str)
            Resp_processes_traces= horzcat(Resp_processes_traces, tempY);
        elseif ~isempty(NP_str)
            Resp_neuropil_traces= horzcat(Resp_neuropil_traces, tempY);
        end
    end
end

% means and SDs
%RCaMP
RespRCmeanTrace = mean(Resp_RCaMP_traces,2);
RespRCSDTrace = std(Resp_RCaMP_traces');

%GcaMP
RespGCmeanTrace = mean(Resp_GCaMP_traces,2);
RespGCSDTrace = std(Resp_GCaMP_traces');

%neurons
RespNmeanTrace = mean(Resp_neuron_traces,2);
RespNSDTrace = std(Resp_neuron_traces');

%neuropil
RespNPmeanTrace = mean(Resp_neuropil_traces,2);
RespNPSDTrace = std(Resp_neuropil_traces');

%somata
RespSmeanTrace = mean(Resp_somata_traces,2);
RespSSDTrace = std(Resp_somata_traces');

%endfeet
RespEFmeanTrace = mean(Resp_EF_traces,2);
RespEFSDTrace = std(Resp_EF_traces');

%processes
RespPmeanTrace = mean(Resp_processes_traces,2);
RespPSDTrace = std(Resp_processes_traces');

%% Plots
% RCaMP vs GCaMP
figure ('name', 'Overlaid traces: All RCaMP responding ROIs')
hold on
axis off
for xROI= 1:size(Resp_RCaMP_traces,2)
    tempY = Resp_RCaMP_traces(:,xROI);
    if length(tempY)>590
        tempY=tempY(1:nframes);
        grey = [0.8,0.8,0.8];
        plot(TimeX,tempY,'Color',grey,'LineWidth',0.01);
    end
end

plot([5 13],[-2 -2], 'k','LineWidth', 2)
plot(TimeX, RespRCmeanTrace, 'k', 'LineWidth',1)

figure ('name', 'Overlaid traces: All GCaMP responding ROIs')
hold on
axis off
for xROI= 1:size(Resp_GCaMP_traces,2)
    tempY = Resp_GCaMP_traces(:,xROI);
    if length(tempY)>590
        tempY=tempY(1:nframes);
        grey = [0.8,0.8,0.8];
        plot(TimeX,tempY,'Color',grey,'LineWidth',0.01);
    end
end

plot([5 13],[-2 -2], 'k','LineWidth', 2)
plot(TimeX, RespGCmeanTrace, 'k', 'LineWidth',1)


figure ('name', 'Overlaid traces: All RCaMP responding ROIs- stim window')
hold on
axis off
for xROI= 1:size(Resp_RCaMP_traces,2)
    tempY = Resp_RCaMP_traces(:,xROI);
    if length(tempY)>590
        tempY=tempY(1:nframes);
        grey = [0.8,0.8,0.8];
        plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.01);
    end
end

plot([5 13],[-2 -2], 'k','LineWidth', 2)
plot(TimeX(1:stimwindow), RespRCmeanTrace(1:stimwindow), 'k', 'LineWidth',1)
plot([-4 -4],[0 5], 'k','LineWidth', 2)
plot([-4 -2],[0 0], 'k','LineWidth', 2)

figure ('name', 'Overlaid traces: All GCaMP responding ROIs- stim window')
hold on
axis off
for xROI= 1:size(Resp_GCaMP_traces,2)
    tempY = Resp_GCaMP_traces(:,xROI);
    if length(tempY)>590
        tempY=tempY(1:nframes);
        grey = [0.8,0.8,0.8];
        plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.01);
    end
end

plot([5 13],[-2 -2], 'k','LineWidth', 2)
plot(TimeX(1:stimwindow), RespGCmeanTrace(1:stimwindow), 'k', 'LineWidth',1)
plot([-4 -4],[0 5], 'k','LineWidth', 2)
plot([-4 -2],[0 0], 'k','LineWidth', 2)

%% Mean RCaMP vs GCaMP
% GCaMP vs RCaMP
figure ('name', 'Mean traces with error bar: GCaMP vs RCaMP responding ROIs- stim window')
hold on
axis off
plot([5 13],[-1 -1], 'k','LineWidth', 2)
shadedErrorBar(TimeX(1:stimwindow),RespRCmeanTrace(1:stimwindow)',RespRCSDTrace(1:stimwindow),'r');
shadedErrorBar(TimeX(1:stimwindow),(RespGCmeanTrace(1:stimwindow)' + 5),RespGCSDTrace(1:stimwindow),'g');
plot([5 5],[-1 10], 'k--','LineWidth', 0.5)
x1 = -0.5;
y1 = 0;
y2 = 5;

txt1 = 'Neurons';
txt2 = 'Astrocytes';

text(x1,y1,txt1,'HorizontalAlignment','right')
text(x1,y2,txt2,'HorizontalAlignment','right')


figure ('name', 'Mean traces only: GCaMP vs RCaMP responding ROIs')
hold on
axis off
plot(TimeX(1:stimwindow),smooth(RespRCmeanTrace(1:stimwindow)',5),'r','LineWidth',1.5);
plot(TimeX(1:stimwindow),smooth(RespGCmeanTrace(1:stimwindow)',5),'g','LineWidth', 1.5);
plot([5 5],[-1 2], 'k--','LineWidth', 1)
plot([5 13],[0 0], 'k','LineWidth', 3)
legend('RCaMP','GCaMP')



%% ROITypes
figure ('name', 'Overlaid traces + mean:All responding ROItype')
subplot(1,5,1)
hold on
axis off
ylim([-5 40]);
for xROI= 1:size(Resp_neuron_traces,2)
    tempY = Resp_neuron_traces(:,xROI);
    if length(tempY)>600
        tempY=tempY(1:nframes);
    end
    grey = [0.8,0.8,0.8];
    plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.01);
    
end
plot([5 13],[-2 -2], 'k','LineWidth', 2)
plot(TimeX(1:stimwindow), RespNmeanTrace(1:stimwindow), 'k', 'LineWidth',1)
plot([5 5],[-1 10], 'k--','LineWidth', 0.5)
x1 = 0;
y1 = 5;
txt1 = 'Neurons';
text(x1,y1,txt1,'HorizontalAlignment','right')

subplot(1,5,2)
hold on
axis off
ylim([-5 40]);
for xROI= 1:size(Resp_neuropil_traces,2)
    tempY = Resp_neuropil_traces(:,xROI);
    if length(tempY)>600
        tempY=tempY(1:nframes);
    end
    grey = [0.8,0.8,0.8];
    plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.01);
end
plot([5 13],[-2 -2], 'k','LineWidth', 2)
plot(TimeX(1:stimwindow), RespNPmeanTrace(1:stimwindow), 'k', 'LineWidth',1)
plot([5 5],[-1 10], 'k--','LineWidth', 0.5)
x1 = 0;
y1 = 5;
txt1 = 'Neuropil';
text(x1,y1,txt1,'HorizontalAlignment','right')


subplot(1,5,4)
hold on
axis off
ylim([-5 40]);
for xROI= 1:size(Resp_EF_traces,2)
    tempY = Resp_EF_traces(:,xROI);
    if length(tempY)>600
        tempY=tempY(1:nframes);
    end
    grey = [0.8,0.8,0.8];
    plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.01);
end
plot([5 13],[-2 -2], 'k','LineWidth', 2)
plot(TimeX(1:stimwindow), RespEFmeanTrace(1:stimwindow), 'k', 'LineWidth',1)
plot([5 5],[-1 10], 'k--','LineWidth', 0.5)
x1 = 0;
y1 = 5;
txt1 = 'Endfeet';
text(x1,y1,txt1,'HorizontalAlignment','right')


subplot(1,5,3)
hold on
axis off
ylim([-5 40]);
for xROI= 1:size(Resp_somata_traces,2)
    tempY = Resp_somata_traces(:,xROI);
    if length(tempY)>600
        tempY=tempY(1:nframes);
    end
    grey = [0.8,0.8,0.8];
    plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.01);
end
plot([5 13],[-2 -2], 'k','LineWidth', 2)
plot(TimeX(1:stimwindow), RespSmeanTrace(1:stimwindow), 'k', 'LineWidth',1)
plot([5 5],[-1 10], 'k--','LineWidth', 0.5)
x1 = 0;
y1 = 5;
txt1 = 'Dendrites';
text(x1,y1,txt1,'HorizontalAlignment','right')


subplot(1,5,5)
hold on
axis off
ylim([-5 40]);
for xROI= 1:size(Resp_processes_traces,2)
    tempY = Resp_processes_traces(:,xROI);
    if length(tempY)>600
        tempY=tempY(1:nframes);
    end
    grey = [0.8,0.8,0.8];
    plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.01);
end
plot([5 13],[-2 -2], 'k','LineWidth', 2)
plot(TimeX(1:stimwindow), RespPmeanTrace(1:stimwindow), 'k', 'LineWidth',1)
plot([5 5],[-1 10], 'k--','LineWidth', 0.5)
x1 = 0;
y1 = 5;
txt1 = 'Processes';
text(x1,y1,txt1,'HorizontalAlignment','right')



%% Histogram of peak time differences between neurons and astrocytes
% compare peak times of responding neurons to responding astrocytes in a field of view
% shift traces of each responding ROI by neuronal peak times
%
% timeDiffs=[];
% ShiftedTraces=[];
%
% uniqueTrials=unique(responders(:,27));
%
% for iTrial=1:size(uniqueTrials,1)
%
%     %find matching trials
%     for xTrial= 1:size(responders,1)
%         Trial_str(xTrial)= strcmp(responders{xTrial, 27},uniqueTrials{iTrial});
%     end
%     trialData=responders(Trial_str,:);
%
%     for xxTrial= 1:size(RespondingROIs,1)
%         Trial_str2(xxTrial)= strcmp(RespondingROIs{xxTrial, 14},uniqueTrials{iTrial});
%     end
%     trialTraces=RespondingROIs(Trial_str2,:);
%
%     %find neuron and neuropil ROIs
%     for iROI=1:size(trialData,1)
%         Nstr{iROI}= regexp(trialData{iROI,15},'R*');
%     end
%
%     NIndx= ~cellfun(@isempty, Nstr);
%     NeuronData=trialData(NIndx,:);
%     NData=[];
%     %find if there are multiple peaks for a given neuron and choose the
%     %faster peak time
%     if length(unique(NeuronData(:,25)))==length(NeuronData(:,25))
%         NData=NeuronData;
%     else
%         uniqueNeurons=unique(NeuronData(:,25)); % list of unique neurons
%         for xxROI= 1:size(uniqueNeurons,1)
%             for kROI= 1:size(NeuronData,1)
%                 NeuronIdx(kROI)= strcmp(NeuronData{kROI,25},uniqueNeurons{xxROI});
%             end
%             singleNeurons=NeuronData(NeuronIdx,:); % index of each individual neuron
%             clear NeuronIdx
%             if size(singleNeurons,1)>1 % if there is more than one peak
%                 [~,fastPeak]=min(cell2mat(singleNeurons(:,6)));
%                 singleNeurons2=singleNeurons(fastPeak,:); % find the fastest peak
%             else
%                 singleNeurons2=singleNeurons;
%             end
%             NData=vertcat(NData,singleNeurons2);
%             clear singleNeurons singleNeurons2
%         end
%     end
%
%     AstrocyteData=trialData(~NIndx,:);
%     if ~isempty(AstrocyteData)
%         %compare neuronal and astrocyte peaks times
%         for iN= 1:size(NData, 1)
%             for iA= 1:size(AstrocyteData,1)
%                 timeComparisons{iA,1}=NData{iN,27}; % unique trial name
%                 timeComparisons{iA,2}=NData{iN,23}; % neuron ROI type
%                 timeComparisons{iA,3}=NData{iN,11}; % neuron name
%                 timeComparisons{iA,4}=NData{iN,6}; % neuronal max peak time
%                 timeComparisons{iA,5}=NData{iN,8}; % neuronal onset time
%                 timeComparisons{iA,6}=AstrocyteData{iA,23}; % astrocyte ROI type
%                 timeComparisons{iA,7}=AstrocyteData{iA,11}; % astrocyte ROI name
%                 timeComparisons{iA,8}=AstrocyteData{iA,6}; % astrocyte peak time
%                 timeComparisons{iA,9}=AstrocyteData{iA,8}; % astrocyte onset
%                 timeComparisons{iA,10}=AstrocyteData{iA,6}-NData{iN,6}; % max peak time (AC) minus max peak time (N)
%                 timeComparisons{iA,11}=AstrocyteData{iA,6}-NData{iN,8}; % max peak time (AC) minus onset time (N)
%                 timeComparisons{iA,12}=AstrocyteData{iA,8}-NData{iN,8}; % onset time (AC) minus onset time (N)
%                 timeComparisons{iA,13}=AstrocyteData{iA,8}-NData{iN,6}; % onset time (AC) minus max peak time (N)
%             end
%             timeDiffs=vertcat(timeDiffs,timeComparisons);
%         end
%     end
%     for iN= 1:size(NData, 1)
%         % shift traces of all ROIs by each neuronal peak time
%         for iTrace=1:size(trialTraces,1)
%             ShiftTrace{iTrace,1}=NData{iN,27}; %unique trial name
%             ShiftTrace{iTrace,2}=NData{iN,23}; % Neuron Type (soma or neuropil)
%             ShiftTrace{iTrace,3}=NData{iN,11}; % Neuron name
%             ShiftTrace{iTrace,4}=NData{iN,6}; % Neuron peak time
%             ShiftTrace{iTrace,5}=NData{iN,8}; % Neuron onset time (half start)
%             ShiftTrace{iTrace,6}=trialTraces{iTrace,13}; % ROI type (ROI 2)
%             ShiftTrace{iTrace,7}=trialTraces{iTrace,1}; % ROI name (ROI 2)
%             ShiftTrace{iTrace,8}=trialTraces{iTrace,8}; % ROI 2 trace
%             ShiftTrace{iTrace,9}=TimeX-(5+NData{iN,6}); % minus neuronal peak time
%             ShiftTrace{iTrace,10}=TimeX-(5+NData{iN,8}); % minus neuronal onset time
%         end
%         ShiftedTraces=vertcat(ShiftedTraces,ShiftTrace);
%         clear timeComparisons ShiftTrace
%     end
%     clear Nstr NIndx NData trialData trialTraces Trial_str2 Trial_str
% end
%
% for iROI=1:size(timeDiffs,1)
%     peakT_peakT(iROI)=timeDiffs{iROI,10};
%     peakT_onset(iROI)=timeDiffs{iROI,11};
%     onset_onset(iROI)=timeDiffs{iROI,12};
%     onset_peakT(iROI)=timeDiffs{iROI,13};
% end
%
% % histograms
% figure('name', 'histogram: peakTime (A)- peakTime (N)')
% histogram(peakT_peakT)
%
% figure('name', 'histogram: peakTime(A)-peakOnset(N)')
% histogram(peakT_onset)
%
% figure('name', 'histogram: peakOnset(A)-peakOnset(N)')
% histogram(onset_onset)
%
% figure('name', 'histogram: peakOnset(A)-peakTime(N)')
% histogram(onset_peakT)
%
% %% ROITypes Shifted Traces- All Astrocyte traces shifted by neuronal peak times
%
% stimwindow2=round(FrameRate*30);
% figure ('name', 'Shifted traces- All trials by Neuron peak time')
% hold on
% axis off
% for xROI= 1:size(ShiftedTraces,1)
%     tempY = ShiftedTraces{xROI,8};
%     if length(tempY)>600
%         tempY=tempY(1:nframes);
%     end
%     tempX = ShiftedTraces{xROI,9}; % shifted by peak time
%     grey = [0.8,0.8,0.8];
%     plot(tempX(1:stimwindow2),tempY(1:stimwindow2)','Color',grey,'LineWidth',0.01);
% end
% plot([0 0],[-1 100], 'k--','LineWidth', 0.5)
% %plot([0 8],[-2 -2], 'k','LineWidth', 2)
%
%
% figure ('name', 'Shifted traces- RCaMP vs GCaMP by Neuron peak time')
% hold on
% axis off
% for xROI= 1:size(ShiftedTraces,1)
%     tempY = ShiftedTraces{xROI,8};
%     if length(tempY)>600
%         tempY=tempY(1:nframes);
%     end
%     tempX = ShiftedTraces{xROI,9}; % shifted by peak time
%     grey = [0.8,0.8,0.8];
%     if strcmp(ShiftedTraces{xROI,6},'Neuron') || strcmp(ShiftedTraces{xROI,6},'Neuropil')
%         plot(tempX,tempY','Color',grey,'LineWidth',0.01);
%     else
%         plot(tempX,(tempY'+50),'Color',grey,'LineWidth',0.01);
%     end
% end
% plot([0 0],[-1 100], 'k--','LineWidth', 0.5)
% %plot([0 8],[-2 -2], 'k','LineWidth', 2)
%
% x1 = -20;
% y1 = 5;
% y2 = 50;
% txt1 = 'RCaMP';
% txt2 = 'GCaMP';
% text(x1,y1,txt1,'HorizontalAlignment','right')
% text(x1,y2,txt2,'HorizontalAlignment','right')
% %% example traces algining astrocytes and neurons ROITypes
% figure ('name', 'Shifted traces- ROITypes by Neuron peak time')
%
% subplot(1,5,1)
% hold on
% axis off
%
% for xROI= 1:size(ShiftedTraces,1)
%     tempY = ShiftedTraces{xROI,8};
%     if length(tempY)>600
%         tempY=tempY(1:nframes);
%     end
%     tempX = ShiftedTraces{xROI,9}; % shifted by peak time
%     grey = [0.8,0.8,0.8];
%     if strcmp(ShiftedTraces{xROI,6},'Neuron')
%         plot(tempX(1:stimwindow),tempY(1:stimwindow)','Color',grey,...
%             'LineWidth',0.01);
%         ylim([-5 30]);
%     end
% end
% plot([0 0],[-1 100], 'k--','LineWidth', 0.5)
% x1 = -10;
% y1 = 5;
% txt1 = 'Neuron';
% text(x1,y1,txt1,'HorizontalAlignment','right')
%
%
% subplot(1,5,2)
% hold on
% axis off
% for xROI= 1:size(ShiftedTraces,1)
%     tempY = ShiftedTraces{xROI,8};
%     if length(tempY)>600
%         tempY=tempY(1:nframes);
%     end
%     tempX = ShiftedTraces{xROI,9};
%     if strcmp(ShiftedTraces{xROI,6},'Neuropil')
%         plot(tempX(1:stimwindow),tempY(1:stimwindow)','Color',grey,...
%             'LineWidth',0.01);
%         ylim([-5 30]);
%         plot([0 0],[-1 100], 'k--','LineWidth', 0.5)
%     end
% end
% x1 = -10;
% y1 = 5;
% txt1 = 'Neuropil';
% text(x1,y1,txt1,'HorizontalAlignment','right')
%
% subplot(1,5,3)
% hold on
% axis off
% for xROI= 1:size(ShiftedTraces,1)
%     tempY = ShiftedTraces{xROI,8};
%     if length(tempY)>600
%         tempY=tempY(1:nframes);
%     end
%     tempX = ShiftedTraces{xROI,9};
%     if strcmp(ShiftedTraces{xROI,6},'Endfeet')
%         plot(tempX(1:stimwindow),tempY(1:stimwindow)','Color',grey,...
%             'LineWidth',0.01);
%         ylim([-5 30]);
%         plot([0 0],[-1 100], 'k--','LineWidth', 0.5)
%     end
% end
% x1 = -10;
% y1 = 5;
% txt1 = 'Endfeet';
% text(x1,y1,txt1,'HorizontalAlignment','right')
%
% subplot(1,5,4)
% hold on
% axis off
% for xROI= 1:size(ShiftedTraces,1)
%     tempY = ShiftedTraces{xROI,8};
%     if length(tempY)>600
%         tempY=tempY(1:nframes);
%     end
%     tempX = ShiftedTraces{xROI,9};
%     if strcmp(ShiftedTraces{xROI,6},'Soma')
%         plot(tempX(1:stimwindow),tempY(1:stimwindow)','Color',grey,...
%             'LineWidth',0.01);
%         ylim([-5 30]);
%         plot([0 0],[-1 100], 'k--','LineWidth', 0.5)
%     end
% end
% x1 = -10;
% y1 = 5;
% txt1 = 'Somata';
% text(x1,y1,txt1,'HorizontalAlignment','right')
%
% subplot(1,5,5)
% hold on
% axis off
% for xROI= 1:size(ShiftedTraces,1)
%     tempY = ShiftedTraces{xROI,8};
%     if length(tempY)>600
%         tempY=tempY(1:nframes);
%     end
%     tempX = ShiftedTraces{xROI,9};
%     if strcmp(ShiftedTraces{xROI,6},'Process')
%
%         plot(tempX(1:stimwindow),tempY(1:stimwindow)','Color',grey,...
%             'LineWidth',0.01);
%         ylim([-5 30]);
%         plot([0 0],[-1 100], 'k--','LineWidth', 0.5)
%     end
% end
% x1 = -10;
% y1 = 5;
% txt1 = 'Processes';
% text(x1,y1,txt1,'HorizontalAlignment','right')
%
%
%
%
%
%
% %% Isolate the early traces from the shifted data
%
% for xROI= 1:size(timeDiffs,1)
%     earlyIdx(xROI)= timeDiffs{xROI,10}<=0;
% end
% early=timeDiffs(earlyIdx,:);
% late=timeDiffs(~earlyIdx,:);
%
% earlyNames={'TrialName','ROITypeX','ROI_X','ROI_X_peakTime','ROI_X_onset',...
%     'ROITypeY','ROI_Y','ROI_Y_peakTime','ROI_Y_onset',...
%     'peak_peak','peak_onset','onset_onset','onset_peak'};

% AllTimeDiffs=vertcat(earlyNames,timeDiffs);
% cd('D:\Data\GCaMP_RCaMP\cyto_GCaMP6s\Results');
% cell2csv('LongStim_TimeDiffs.csv',AllTimeDiffs);

%% Plot early ROI traces

% all GCaMP ROIs
GC_idx=find(strcmp(RespondingROIs(:,3),'GCaMP'));
GCaMP_ROI= RespondingROIs(GC_idx,:);

%remove duplicate traces
wd=GCaMP_ROI;
[~,idx]=unique(wd(:,15));
GCaMP_ROI=wd(idx,:);

numROI=length(unique(GCaMP_ROI(:,15)));

for iROI=1:size(GCaMP_ROI,1)
    %area under the curve
    tempY=GCaMP_ROI{iROI,8};
    %first 3 sec after stim onset
    x1=round(FrameRate*5);
    x2=round(FrameRate*6);
    GCaMP_ROI{iROI,17}=trapz(tempY(x1:x2));
end

%%
% find the ROIs that have a decent area under the curve for the first 1 sec
% of window
AUC_idx=find(cell2mat(GCaMP_ROI(:,17))>5);
AUC_idx2=find(cell2mat(GCaMP_ROI(:,17))<=5);
earlyGC=GCaMP_ROI(AUC_idx,:);
lateGC=GCaMP_ROI(AUC_idx2,:);

earlyNames=earlyGC(:,15);
%cd('D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results');
%cell2csv('Early_LckGCaMPROIs.csv',earlyNames);


allearlyGC=[];

figure('name', 'AUC greater than 5 in first 1 sec')
hold on
axis off

for xROI= 1:size(earlyGC,1)
    tempY = earlyGC{xROI,8};
    if length(tempY)>600
        tempY=tempY(1:592);
    end
    grey = [0.8,0.8,0.8];
    plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.1);
    allearlyGC= horzcat(allearlyGC, tempY);
end
%plot([5 13],[-2 -2], 'k','LineWidth', 2)
plot([5 5],[-1 20], 'k--','LineWidth', 0.5)
meanEarlyTrace = mean(allearlyGC,2);
plot(TimeX(1:stimwindow), meanEarlyTrace(1:stimwindow), 'k', 'LineWidth',1)
%plot([-4 -4],[0 5], 'k','LineWidth', 2)
%plot([-4 -2],[0 0], 'k','LineWidth', 2)

allLateGC=[];
figure('name', 'late- responding ROIs not in early group')
hold on
axis off
for xROI= 1:size(lateGC,1)
    tempY = lateGC{xROI,8};
    if length(tempY)>600
        tempY=tempY(1:592);
    end
    grey = [0.8,0.8,0.8];
    plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.1);
    allLateGC= horzcat(allLateGC, tempY);
end
plot([5 13],[-2 -2], 'k','LineWidth', 2)
plot([5 5],[-1 20], 'k--','LineWidth', 0.5)
meanLateTrace = mean(allLateGC,2);
plot(TimeX(1:stimwindow), meanLateTrace(1:stimwindow), 'k', 'LineWidth',1)

% lateExport=allLateGC(1:355,:);
% lateExport2=TimeX(1:355);
% lateExport2=[lateExport2,lateExport];

%%
figure('name', 'Mean Late & early traces')
hold on
axis off

plot(TimeX(1:stimwindow), smooth(meanLateTrace(1:stimwindow),3), 'g', 'LineWidth',1)
plot(TimeX(1:stimwindow), smooth(meanEarlyTrace(1:stimwindow),3), 'b', 'LineWidth',1)
% %RCaMP
RespRCmeanTrace = mean(Resp_RCaMP_traces,2);
plot(TimeX(1:stimwindow), smooth(RespRCmeanTrace(1:stimwindow),3), 'r', 'LineWidth',1)
plot([5 13],[-0.5 -0.5], 'k','LineWidth', 2)
plot([5 5],[-0.5 1], 'k--','LineWidth', 0.5)
plot([-4 -4],[0 0.5], 'k','LineWidth', 2)
plot([-4 -2],[0 0], 'k','LineWidth', 2)


figure('name', 'Mean Early traces')
hold on
axis off

plot(TimeX(1:stimwindow), meanEarlyTrace(1:stimwindow), 'g', 'LineWidth',1)
% %RCaMP
RespRCmeanTrace = mean(Resp_RCaMP_traces,2);
plot(TimeX(1:stimwindow), RespRCmeanTrace(1:stimwindow), 'r', 'LineWidth',1)
plot([5 13],[-0.5 -0.5], 'k','LineWidth', 2)
plot([5 5],[-0.5 1], 'k--','LineWidth', 0.5)
plot([-1 -1],[0 0.5], 'k','LineWidth', 1)
%% 3d plot
x1=zeros(1,size(allearlyGC,2));
y1=ones(1,size(allearlyGC,2));
z1=1:size(allearlyGC,2);

earlywindow=allearlyGC(1:stimwindow,:);
figure
h=ribbon(TimeX(1:stimwindow),earlywindow, 0.5);
set(h, {'CData'}, get(h,'ZData'), 'FaceColor','interp','MeshStyle','column')
%hold on
%plot3(x1,y1,z1, 'k','LineWidth',2)
%zlim([0 15])

% figure
% h=ribbon(earlywindow', TimeX(1:stimwindow),0.5);
% set(h, {'CData'}, get(h,'ZData'), 'FaceColor','interp','MeshStyle','column')
% zlim([0 15])

%%
%map=[0 0.5 0.5;0.5 1 1];
h=colormap(gca,'parula');
HeatMap(earlywindow','colormap',h);


%save names
earlyGCnames=earlyGC(:,15);
earlyGCnames2=vertcat('names',earlyGCnames);
cd('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results');
cell2csv('earlyGC_byAUC.csv',earlyGCnames2)

%% subsetting data test plots
allearlyGC=[];

figure('name', 'AUC greater than 5 in first 1 sec')
hold on
axis off

for xROI= 1:size(earlyGC,1)
    tempY = earlyGC{xROI,8};
    if length(tempY)>600
        tempY=tempY(1:592);
    end
    grey = [0.8,0.8,0.8];
    plot(TimeX(55:100),tempY(55:100),'Color',grey,'LineWidth',0.1);
    allearlyGC= horzcat(allearlyGC, tempY);
end
plot([5 5],[-1 20], 'k--','LineWidth', 0.5)
meanEarlyTrace = mean(allearlyGC,2);
plot(TimeX(55:100), meanEarlyTrace(55:100), 'k', 'LineWidth',1)


figure('name', 'responding RCaMP + Early traces GCaMP')
hold on
axis off
ylim([-2 25])
for xROI= 1:size(earlyGC,1)
    tempY = earlyGC{xROI,8};
    grey = [0.8,0.8,0.8];
    plot(TimeX(55:100),(tempY(55:100)')+15,'Color',grey,'LineWidth',0.1);
end
for xROI= 1:size(Resp_RCaMP_traces,2)
    tempY = Resp_RCaMP_traces(:,xROI);
    grey = [0.8,0.8,0.8];
    plot(TimeX(55:100),tempY(55:100)','Color',grey,'LineWidth',0.1);
end
plot([5 5],[-1 100], 'k--','LineWidth', 0.5)
plot([5 13],[-2 -2], 'k','LineWidth', 2)



