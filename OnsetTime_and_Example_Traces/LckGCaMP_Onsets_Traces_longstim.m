
clearvars
close all

%% Load data

% time windows based on stimulation
NOnsetWindow= 2; %1 for short stim % neuronal onset times
AOnsetWindow= 12; % 5 for short stim % astrocyte onset times
Fast_AOnsetWindow=1;
Fast_NOnsetWindow=1;
NPTWindow= 9; % one second longer than stimulation for peak times
APTWindow= 15; % astrocyte longer than stimulation for peak times

stimwindow=20; % 5 s baseline, 15 s imaging 

% save files names
%saveFiles1='E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\Lck_longstim_firstonset_comparisons.mat';
%saveFiles2='E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\Lck_longstim_firstonset_comparisons.csv';
saveFiles3= 'E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\Lck_longstim_onset&AUC.csv';
%saveFiles3= 'E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\Lck_longstim_onset_2sbeforeStim.csv';

%AstrocyteExcelFile='E:\Data\Two_Photon_Data\GCaMP_RCaMP\Manuscript\Figures\DataForTraceOverlay_Heatmaps\cytovslck_Nostimvsstim_comparison\cyto-AstrocyteTraces_Stim.xlsx';
%NeuronalExcelFile ='E:\Data\Two_Photon_Data\GCaMP_RCaMP\Manuscript\Figures\DataForTraceOverlay_Heatmaps\cytovslck_Nostimvsstim_comparison\cyto-NeuronalTraces_Stim.xlsx';
%peak data
load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_2D_longstim_05_04_2017.mat');
%load('D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_2D_longstim_28_04_2017.mat');

% Load trace data
load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_longstim_05_04_2017.mat');


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

%baselineCorrectedTime=TimeX-3;%1.9908;
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
        Onsets=find_first_onset_time(baselineCorrectedTime(10:end), trace(10:592),2.5,2);
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
    


% %% compare peak onsets for each RCaMP and GCaMP ROIs
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
%         % find similar astrocyte peak times
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
% %                                 OnsetTimeComparisons{iA,18}=TimeDiff2; %neurons-astrocytes
% %                                 OnsetTimeComparisons{iA,19}=NeuronalData{iN,8}; %neuronal trace
% %                                 OnsetTimeComparisons{iA,20}=AstroData{iA,8};
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
% ACIdx= find(~cellfun(@isempty, regexp(Shortstim(:,3), 'GCaMP')));
% figure('name', 'Astrocyte onset times')
% histogram(cell2mat(Shortstim(ACIdx,16)), 'binwidth',0.8450)
% 
% NIdx= find(~cellfun(@isempty, regexp(Shortstim(:,3), 'RCaMP')));
% figure('name', 'Neuronal onset times')
% histogram(cell2mat(Shortstim(NIdx,16)), 'binwidth', 0.845)
% 
%  figure('name', 'Astrocyte AUC')
% histogram(cell2mat(Shortstim(ACIdx,17)))
% 
% figure('name', 'Neuronal AUC')
% histogram(cell2mat(Shortstim(NIdx,17)))

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
respOTIdx2(iROI)=~isempty(find(GCaMP{iROI,16}>0 && GCaMP{iROI,16}<AOnsetWindow));
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
scatter(cell2mat(GrespOT(:,16)), cell2mat(GrespOT(:,17)));

figure('name','1 sec auc vs onset time');
hold on
scatter(cell2mat(Shortstim(:,16)), cell2mat(Shortstim(:,17)));

% figure('name','20 sec auc vs onset time');
% hold on
% scatter(cell2mat(Shortstim(:,16)), cell2mat(Shortstim(:,18)));

for iROI=1:length(GrespOT)
    fastIdx(iROI)=~isempty(find(GrespOT{iROI,16}>0 && GrespOT{iROI,16}<=2));
    slowIdx(iROI)=~isempty(find(GrespOT{iROI,16}>2));
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
xlim([-1 25]);
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



