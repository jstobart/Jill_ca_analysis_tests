%% LONG STIM (8sec)
clearvars
close all

%load data traces
 load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\S&LStim_LckGC&RC_traces_17_02_2017.mat');
% %load('D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\S&LStim_LckGC&RC_traces_17_02_2017.mat');
% 
% Long_Short=All_traces;
% 
 load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LStim_LckGC&RC_traces_17_02_2017.mat');
% %load('D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LStim_LckGC&RC_traces_17_02_2017.mat');
% 
% Long=All_traces;
% 
% All_traces=vertcat(Long_Short,Long);
%load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_traces_2D_longstim_05_04_2017.mat')
FrameRate=11.84;
stimwindow=round(FrameRate*30);

searchString='shortstim';

for iROI=1:length(All_traces)
    Stim_str{iROI} = strfind(All_traces{iROI,6},searchString);
    str_idx(iROI)=~isempty(Stim_str{iROI});
    
    %find ROITypes
    N_str= strfind(All_traces{iROI, 1},'N');
    S_str= strfind(All_traces{iROI, 1},'D');
    EF_str= strfind(All_traces{iROI, 1},'E');
    P_str= strfind(All_traces{iROI, 1},'r');
    NP_str= strfind(All_traces{iROI, 1},'np');
    
    if ~isempty(N_str)
        All_traces{iROI,13}='Neuron';
    elseif ~isempty(S_str)
        All_traces{iROI,13}='Dendrite';
    elseif ~isempty(EF_str)
        All_traces{iROI,13}='Endfeet';
    elseif ~isempty(P_str)
        All_traces{iROI,13}='Process';
    elseif ~isempty(NP_str)
        All_traces{iROI,13}='Neuropil';
    end
    
    % make new names
    All_traces{iROI,14}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,2});
    All_traces{iROI,15}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,2}, '_',All_traces{iROI,1});
    All_traces{iROI,16}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,2}, '_',All_traces{iROI,1},'_',All_traces{iROI,6});
    
    
end

data_traces = All_traces(str_idx',:);

for iROI=1:length(data_traces)
    %if a process ROI is overlapping with a soma or process, exclude it
    %from data
    %All_traces{iROI,11}=0;
    nonOverlapIdx(iROI)=~ischar(data_traces{iROI,12});
    
end

% remove overlapping processes
data_traces = data_traces(nonOverlapIdx',:);

% get rid of astrocyte neuropil traces
for xROI=1:size(data_traces,1)
    GCNP2(xROI)=(strcmp(data_traces{xROI,3},'GCaMP')&& strcmp(data_traces{xROI,13},'Neuropil'));
end

data_traces=data_traces(~GCNP2',:);

% get info for plotting
nROIs=length(data_traces);

nframes=592;

TimeX(1:nframes) = (1:nframes)/FrameRate;


%% Plot all traces for the entire trial

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

plot([5 6],[-2 -2], 'k','LineWidth', 2)
meanTrace = mean(traces,2);
plot(TimeX, meanTrace, 'k', 'LineWidth',1)
% plot([0 0],[18.5 19.5], 'k','LineWidth', 2)


%%  Plot only the responding neurons and astrocytes from the same field of view


XLfile = 'E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\respondingROIs_shortstim.xlsx';
%XLfile = 'E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\respondingROIs_longstim.xlsx';

%XLfile = 'D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\respondingROIs_longstim.xlsx';

[~, ~, data] = xlsread(XLfile); %,'Sheet1'); %reads the scoresheet and saves all data in a cell array
responders=data(2:end,:);

for xROI=1:length(responders)
    GCNP(xROI)=(strcmp(responders{xROI,15},'GCaMP')&& strcmp(responders{xROI,23},'Neuropil'));
end

responders=responders(~GCNP',:);

RespondingROIs=[];
for xROI=1:length(responders)
    currentROI=responders{xROI, 25};
    for iROI = 1:size(data_traces,1)
        ROIIndex=strcmp(data_traces{iROI,15},currentROI);
        if ROIIndex
            RespondingROIs=vertcat(RespondingROIs,data_traces(iROI,1:end));
        end
    end
end




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

plot([5 6],[-2 -2], 'k','LineWidth', 2)
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

plot([5 6],[-2 -2], 'k','LineWidth', 2)
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

plot([5 6],[-2 -2], 'k','LineWidth', 2)
plot(TimeX(1:stimwindow), RespRCmeanTrace(1:stimwindow), 'k', 'LineWidth',1)

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

plot([5 6],[-2 -2], 'k','LineWidth', 2)
plot(TimeX(1:stimwindow), RespGCmeanTrace(1:stimwindow), 'k', 'LineWidth',1)

%% Mean RCaMP vs GCaMP
% GCaMP vs RCaMP
figure ('name', 'Mean traces with error bar: GCaMP vs RCaMP responding ROIs- stim window')
hold on
axis off
plot([5 6],[-1 -1], 'k','LineWidth', 2)
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
plot([5 6],[0 0], 'k','LineWidth', 3)
legend('RCaMP','GCaMP')



%% ROITypes
figure ('name', 'Overlaid traces + mean:All responding ROItype')
subplot(1,5,1)
hold on
axis off
ylim([-5 15]);
for xROI= 1:size(Resp_neuron_traces,2)
    tempY = Resp_neuron_traces(:,xROI);
    if length(tempY)>590
        tempY=tempY(1:nframes);
        grey = [0.8,0.8,0.8];
        plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.01);
    end
    
end
plot([5 6],[-2 -2], 'k','LineWidth', 2)
plot(TimeX(1:stimwindow), RespNmeanTrace(1:stimwindow), 'k', 'LineWidth',1)
plot([5 5],[-1 10], 'k--','LineWidth', 0.5)
x1 = 0;
y1 = 5;
txt1 = 'Neurons';
text(x1,y1,txt1,'HorizontalAlignment','right')

subplot(1,5,2)
hold on
axis off
ylim([-5 15]);
for xROI= 1:size(Resp_neuropil_traces,2)
    tempY = Resp_neuropil_traces(:,xROI);
    if length(tempY)>590
        tempY=tempY(1:nframes);
        grey = [0.8,0.8,0.8];
        plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.01);
    end
end
plot([5 6],[-2 -2], 'k','LineWidth', 2)
plot(TimeX(1:stimwindow), RespNPmeanTrace(1:stimwindow), 'k', 'LineWidth',1)
plot([5 5],[-1 10], 'k--','LineWidth', 0.5)
x1 = 0;
y1 = 5;
txt1 = 'Neuropil';
text(x1,y1,txt1,'HorizontalAlignment','right')


subplot(1,5,4)
hold on
axis off
ylim([-5 15]);
for xROI= 1:size(Resp_EF_traces,2)
    tempY = Resp_EF_traces(:,xROI);
    if length(tempY)>590
        tempY=tempY(1:nframes);
        grey = [0.8,0.8,0.8];
        plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.01);
    end
end
plot([5 6],[-2 -2], 'k','LineWidth', 2)
plot(TimeX(1:stimwindow), RespEFmeanTrace(1:stimwindow), 'k', 'LineWidth',1)
plot([5 5],[-1 10], 'k--','LineWidth', 0.5)
x1 = 0;
y1 = 5;
txt1 = 'Endfeet';
text(x1,y1,txt1,'HorizontalAlignment','right')


subplot(1,5,3)
hold on
axis off
ylim([-5 15]);
for xROI= 1:size(Resp_somata_traces,2)
    tempY = Resp_somata_traces(:,xROI);
    if length(tempY)>590
        tempY=tempY(1:nframes);
        grey = [0.8,0.8,0.8];
        plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.01);
    end
end
plot([5 6],[-2 -2], 'k','LineWidth', 2)
plot(TimeX(1:stimwindow), RespSmeanTrace(1:stimwindow), 'k', 'LineWidth',1)
plot([5 5],[-1 10], 'k--','LineWidth', 0.5)
x1 = 0;
y1 = 5;
txt1 = 'Dendrites';
text(x1,y1,txt1,'HorizontalAlignment','right')


subplot(1,5,5)
hold on
axis off
ylim([-5 15]);
for xROI= 1:size(Resp_processes_traces,2)
    tempY = Resp_processes_traces(:,xROI);
    if length(tempY)>590
        tempY=tempY(1:nframes);
        grey = [0.8,0.8,0.8];
        plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.01);
    end
end
plot([5 6],[-2 -2], 'k','LineWidth', 2)
plot(TimeX(1:stimwindow), RespPmeanTrace(1:stimwindow), 'k', 'LineWidth',1)
plot([5 5],[-1 10], 'k--','LineWidth', 0.5)
x1 = 0;
y1 = 5;
txt1 = 'Processes';
text(x1,y1,txt1,'HorizontalAlignment','right')

%% mean ROItype traces

figure ('name', 'Mean traces plus errorbar: ROITypes responding- stim window')
hold on
axis off
plot([5 6],[-1 -1], 'k','LineWidth', 2)
shadedErrorBar(TimeX(1:stimwindow),RespNmeanTrace(1:stimwindow)',RespNSDTrace(1:stimwindow),'r');
shadedErrorBar(TimeX(1:stimwindow),(RespNPmeanTrace(1:stimwindow)' + 3),RespNPSDTrace(1:stimwindow),'b');
shadedErrorBar(TimeX(1:stimwindow),(RespSmeanTrace(1:stimwindow)' + 6),RespSSDTrace(1:stimwindow),'g');
shadedErrorBar(TimeX(1:stimwindow),(RespEFmeanTrace(1:stimwindow)' + 9),RespEFSDTrace(1:stimwindow),'c');
shadedErrorBar(TimeX(1:stimwindow),(RespPmeanTrace(1:stimwindow)' + 15),RespPSDTrace(1:stimwindow),'m');
plot([5 5],[-1 20], 'k--','LineWidth', 0.5)
x1 = -0.5;
y1 = 0;
y2 = 3;
y3 = 6;
y4 = 9;
y5 = 15;

txt1 = 'Neuron';
txt2 = 'Neuropil';
txt3 = 'Dendrites';
txt4 = 'Endfeet';
txt5 = 'Processes';
text(x1,y1,txt1,'HorizontalAlignment','right')
text(x1,y2,txt2,'HorizontalAlignment','right')
text(x1,y3,txt3,'HorizontalAlignment','right')
text(x1,y4,txt4,'HorizontalAlignment','right')
text(x1,y5,txt5,'HorizontalAlignment','right')


%ROITypes
figure ('name', 'Mean traces only: responding ROITypes- stim window')
hold on
axis off
plot(TimeX(1:stimwindow),smooth(RespNmeanTrace(1:stimwindow)',5),'r','LineWidth', 1.5);
plot(TimeX(1:stimwindow),smooth(RespNPmeanTrace(1:stimwindow)',5),'b','LineWidth', 1.5);
plot(TimeX(1:stimwindow),smooth(RespSmeanTrace(1:stimwindow)',5),'g','LineWidth', 1.5);
plot(TimeX(1:stimwindow),smooth(RespEFmeanTrace(1:stimwindow)',5),'k','LineWidth', 1.5);
plot(TimeX(1:stimwindow),smooth(RespPmeanTrace(1:stimwindow)',5),'m','LineWidth', 1.5);
plot([5 5],[-1 2], 'k--','LineWidth', 1)
plot([5 6],[0 0], 'k','LineWidth', 3)
legend('Neuron','Neuropil','Dendrites','Endfeet','Processes')



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

%% Histogram of peak time differences between neurons and astrocytes
% compare peak times of responding neurons to responding astrocytes in a field of view
% shift traces of each responding ROI by neuronal peak times

timeDiffs=[];
ShiftedTraces=[];

uniqueTrials=unique(responders(:,27));

for iTrial=1:size(uniqueTrials,1)
    
    %find matching trials
    for xTrial= 1:size(responders,1)
        Trial_str(xTrial)= strcmp(responders{xTrial, 27},uniqueTrials{iTrial});
    end
    trialData=responders(Trial_str,:);
    
    for xxTrial= 1:size(RespondingROIs,1)
        Trial_str2(xxTrial)= strcmp(RespondingROIs{xxTrial, 14},uniqueTrials{iTrial});
    end
    trialTraces=RespondingROIs(Trial_str2,:);
    
    %find neuron and neuropil ROIs
    for iROI=1:size(trialData,1)
        Nstr{iROI}= regexp(trialData{iROI,15},'R*');
    end
    
    NIndx= ~cellfun(@isempty, Nstr);
    NeuronData=trialData(NIndx,:);
    NData=[];
    %find if there are multiple peaks for a given neuron and choose the
    %faster peak time
    if length(unique(NeuronData(:,25)))==length(NeuronData(:,25))
        NData=NeuronData;
    else
        uniqueNeurons=unique(NeuronData(:,25)); % list of unique neurons
        for xxROI= 1:size(uniqueNeurons,1)
            for kROI= 1:size(NeuronData,1)
                NeuronIdx(kROI)= strcmp(NeuronData{kROI,25},uniqueNeurons{xxROI});
            end
            singleNeurons=NeuronData(NeuronIdx,:); % index of each individual neuron
            clear NeuronIdx
            if size(singleNeurons,1)>1 % if there is more than one peak
                [~,fastPeak]=min(cell2mat(singleNeurons(:,6)));
                singleNeurons2=singleNeurons(fastPeak,:); % find the fastest peak
            else
                singleNeurons2=singleNeurons;
            end
            NData=vertcat(NData,singleNeurons2);
            clear singleNeurons singleNeurons2
        end
    end
    
    AstrocyteData=trialData(~NIndx,:);
    if ~isempty(AstrocyteData)
        %compare neuronal and astrocyte peaks times
        for iN= 1:size(NData, 1)
            for iA= 1:size(AstrocyteData,1)
                timeComparisons{iA,1}=NData{iN,27}; % unique trial name
                timeComparisons{iA,2}=NData{iN,23}; % neuron ROI type
                timeComparisons{iA,3}=NData{iN,11}; % neuron name
                timeComparisons{iA,4}=NData{iN,6}; % neuronal max peak time
                timeComparisons{iA,5}=NData{iN,8}; % neuronal onset time
                timeComparisons{iA,6}=AstrocyteData{iA,23}; % astrocyte ROI type
                timeComparisons{iA,7}=AstrocyteData{iA,11}; % astrocyte ROI name
                timeComparisons{iA,8}=AstrocyteData{iA,6}; % astrocyte peak time
                timeComparisons{iA,9}=AstrocyteData{iA,8}; % astrocyte onset
                timeComparisons{iA,10}=AstrocyteData{iA,6}-NData{iN,6}; % max peak time (AC) minus max peak time (N)
                timeComparisons{iA,11}=AstrocyteData{iA,6}-NData{iN,8}; % max peak time (AC) minus onset time (N)
                timeComparisons{iA,12}=AstrocyteData{iA,8}-NData{iN,8}; % onset time (AC) minus onset time (N)
                timeComparisons{iA,13}=AstrocyteData{iA,8}-NData{iN,6}; % onset time (AC) minus max peak time (N)
            end
            timeDiffs=vertcat(timeDiffs,timeComparisons);
        end
    end
    for iN= 1:size(NData, 1)
        % shift traces of all ROIs by each neuronal peak time
        for iTrace=1:size(trialTraces,1)
            ShiftTrace{iTrace,1}=NData{iN,27}; %unique trial name
            ShiftTrace{iTrace,2}=NData{iN,23}; % Neuron Type (soma or neuropil)
            ShiftTrace{iTrace,3}=NData{iN,11}; % Neuron name
            ShiftTrace{iTrace,4}=NData{iN,6}; % Neuron peak time
            ShiftTrace{iTrace,5}=NData{iN,8}; % Neuron onset time (half start)
            ShiftTrace{iTrace,6}=trialTraces{iTrace,13}; % ROI type (ROI 2)
            ShiftTrace{iTrace,7}=trialTraces{iTrace,1}; % ROI name (ROI 2)
            ShiftTrace{iTrace,8}=trialTraces{iTrace,8}; % ROI 2 trace
            ShiftTrace{iTrace,9}=TimeX-(5+NData{iN,6}); % minus neuronal peak time
            ShiftTrace{iTrace,10}=TimeX-(5+NData{iN,8}); % minus neuronal onset time
        end
        ShiftedTraces=vertcat(ShiftedTraces,ShiftTrace);
        clear timeComparisons ShiftTrace
    end
    clear Nstr NIndx NData trialData trialTraces Trial_str2 Trial_str
end

for iROI=1:size(timeDiffs,1)
    peakT_peakT(iROI)=timeDiffs{iROI,10};
    peakT_onset(iROI)=timeDiffs{iROI,11};
    onset_onset(iROI)=timeDiffs{iROI,12};
    onset_peakT(iROI)=timeDiffs{iROI,13};
end

% histograms
figure('name', 'histogram: peakTime (A)- peakTime (N)')
histogram(peakT_peakT)

figure('name', 'histogram: peakTime(A)-peakOnset(N)')
histogram(peakT_onset)

figure('name', 'histogram: peakOnset(A)-peakOnset(N)')
histogram(onset_onset)

figure('name', 'histogram: peakOnset(A)-peakTime(N)')
histogram(onset_peakT)

%% ROITypes Shifted Traces- All Astrocyte traces shifted by neuronal peak times

stimwindow2=round(FrameRate*30);
figure ('name', 'Shifted traces- All trials by Neuron peak time')
hold on
axis off
for xROI= 1:size(ShiftedTraces,1)
    tempY = ShiftedTraces{xROI,8};
    if length(tempY)>590
        tempY=tempY(1:nframes);
        tempX = ShiftedTraces{xROI,9}; % shifted by peak time
        grey = [0.8,0.8,0.8];
        plot(tempX(1:stimwindow2),tempY(1:stimwindow2)','Color',grey,'LineWidth',0.01);
    end
end
plot([0 0],[-1 100], 'k--','LineWidth', 0.5)
%plot([0 8],[-2 -2], 'k','LineWidth', 2)


figure ('name', 'Shifted traces- RCaMP vs GCaMP by Neuron peak time')
hold on
axis off
for xROI= 1:size(ShiftedTraces,1)
    tempY = ShiftedTraces{xROI,8};
    if length(tempY)>590
        tempY=tempY(1:nframes);
        tempX = ShiftedTraces{xROI,9}; % shifted by peak time
        grey = [0.8,0.8,0.8];
        if strcmp(ShiftedTraces{xROI,6},'Neuron') || strcmp(ShiftedTraces{xROI,6},'Neuropil')
            plot(tempX,tempY','Color',grey,'LineWidth',0.01);
        else
            plot(tempX,(tempY'+50),'Color',grey,'LineWidth',0.01);
        end
    end
end
plot([0 0],[-1 100], 'k--','LineWidth', 0.5)
%plot([0 8],[-2 -2], 'k','LineWidth', 2)

x1 = -20;
y1 = 5;
y2 = 50;
txt1 = 'RCaMP';
txt2 = 'GCaMP';
text(x1,y1,txt1,'HorizontalAlignment','right')
text(x1,y2,txt2,'HorizontalAlignment','right')
%% example traces algining astrocytes and neurons ROITypes
figure ('name', 'Shifted traces- ROITypes by Neuron peak time')

subplot(1,5,1)
hold on
axis off

for xROI= 1:size(ShiftedTraces,1)
    tempY = ShiftedTraces{xROI,8};
    if length(tempY)>590
        tempY=tempY(1:nframes);
        tempX = ShiftedTraces{xROI,9}; % shifted by peak time
        grey = [0.8,0.8,0.8];
        if strcmp(ShiftedTraces{xROI,6},'Neuron')
            plot(tempX(1:stimwindow),tempY(1:stimwindow)','Color',grey,...
                'LineWidth',0.01);
            ylim([-5 30]);
        end
    end
end
plot([0 0],[-1 100], 'k--','LineWidth', 0.5)
x1 = -10;
y1 = 5;
txt1 = 'Neuron';
text(x1,y1,txt1,'HorizontalAlignment','right')


subplot(1,5,2)
hold on
axis off
for xROI= 1:size(ShiftedTraces,1)
    tempY = ShiftedTraces{xROI,8};
    if length(tempY)>590
        tempY=tempY(1:nframes);
        tempX = ShiftedTraces{xROI,9};
        if strcmp(ShiftedTraces{xROI,6},'Neuropil')
            plot(tempX(1:stimwindow),tempY(1:stimwindow)','Color',grey,...
                'LineWidth',0.01);
            ylim([-5 30]);
            plot([0 0],[-1 100], 'k--','LineWidth', 0.5)
        end
    end
end
x1 = -10;
y1 = 5;
txt1 = 'Neuropil';
text(x1,y1,txt1,'HorizontalAlignment','right')

subplot(1,5,3)
hold on
axis off
for xROI= 1:size(ShiftedTraces,1)
    tempY = ShiftedTraces{xROI,8};
    if length(tempY)>590
        tempY=tempY(1:nframes);
        tempX = ShiftedTraces{xROI,9};
        if strcmp(ShiftedTraces{xROI,6},'Endfeet')
            plot(tempX(1:stimwindow),tempY(1:stimwindow)','Color',grey,...
                'LineWidth',0.01);
            ylim([-5 30]);
            plot([0 0],[-1 100], 'k--','LineWidth', 0.5)
        end
    end
end
x1 = -10;
y1 = 5;
txt1 = 'Endfeet';
text(x1,y1,txt1,'HorizontalAlignment','right')

subplot(1,5,4)
hold on
axis off
for xROI= 1:size(ShiftedTraces,1)
    tempY = ShiftedTraces{xROI,8};
    if length(tempY)>590
        tempY=tempY(1:nframes);
        tempX = ShiftedTraces{xROI,9};
        if strcmp(ShiftedTraces{xROI,6},'Soma')
            plot(tempX(1:stimwindow),tempY(1:stimwindow)','Color',grey,...
                'LineWidth',0.01);
            ylim([-5 30]);
            plot([0 0],[-1 100], 'k--','LineWidth', 0.5)
        end
    end
end
x1 = -10;
y1 = 5;
txt1 = 'Somata';
text(x1,y1,txt1,'HorizontalAlignment','right')

subplot(1,5,5)
hold on
axis off
for xROI= 1:size(ShiftedTraces,1)
    tempY = ShiftedTraces{xROI,8};
    if length(tempY)>590
        tempY=tempY(1:nframes);
        tempX = ShiftedTraces{xROI,9};
        if strcmp(ShiftedTraces{xROI,6},'Process')
            
            plot(tempX(1:stimwindow),tempY(1:stimwindow)','Color',grey,...
                'LineWidth',0.01);
            ylim([-5 30]);
            plot([0 0],[-1 100], 'k--','LineWidth', 0.5)
        end
    end
end
x1 = -10;
y1 = 5;
txt1 = 'Processes';
text(x1,y1,txt1,'HorizontalAlignment','right')






%% Isolate the early traces from the shifted data

for xROI= 1:size(timeDiffs,1)
    earlyIdx(xROI)= timeDiffs{xROI,10}<=0;
end
early=timeDiffs(earlyIdx,:);
late=timeDiffs(~earlyIdx,:);

earlyNames={'TrialName','ROITypeX','ROI_X','ROI_X_peakTime','ROI_X_onset',...
    'ROITypeY','ROI_Y','ROI_Y_peakTime','ROI_Y_onset',...
    'peak_peak','peak_onset','onset_onset','onset_peak'};

% AllTimeDiffs=vertcat(earlyNames,timeDiffs);
% cd('D:\Data\GCaMP_RCaMP\cyto_GCaMP6s\Results');
% cell2csv('LongStim_TimeDiffs.csv',AllTimeDiffs);

%% calculate area under the curve for the first 1 sec
GC_idx=find(strcmp(RespondingROIs(:,3),'GCaMP'));
GCaMP_ROI= RespondingROIs(GC_idx,:);

%remove duplicate traces
wd=GCaMP_ROI;
[~,idx]=unique(wd(:,15));
GCaMP_ROI=wd(idx,:);

for iROI=1:size(GCaMP_ROI,1)
    %area under the curve
    tempY=GCaMP_ROI{iROI,8};
    %first 2 sec after stim onset
    x1=round(FrameRate*5);
    x2=round(FrameRate*6);
    GCaMP_ROI{iROI,17}=trapz(tempY(x1:x2));
end

% find the ROIs that have a decent area under the curve for the first 1 sec
% of window
AUC_idx=find(cell2mat(GCaMP_ROI(:,17))>5);
AUC_idx2=find(cell2mat(GCaMP_ROI(:,17))<=5);
earlyGC=GCaMP_ROI(AUC_idx,:);
lateGC=GCaMP_ROI(AUC_idx2,:);

earlyNames=earlyGC(:,15);
%cd('D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results');
%cell2csv('Early_LckGCaMPROIs.csv',earlyNames);


% control
% check area outside the stim window
for iROI=1:size(RespondingROIs,1)
    %area under the curve
    tempY=RespondingROIs{iROI,8};
    %first 2 sec after stim onset
    if length(tempY)>590
    x1=round(FrameRate*30);
    x2=round(FrameRate*31);
    RespondingROIs{iROI,17}=trapz(tempY(x1:x2));
    end
end

% find the ROIs that have a decent area under the curve for the first 1 sec
% of window
gcamp_idx2=find(strcmp(RespondingROIs(:,3),'GCaMP'));
controlGC=RespondingROIs(gcamp_idx2,:);

controlGC_idx=find(cell2mat(controlGC(:,17))>5);
controlGC=controlGC(controlGC_idx,:);

earlyNames=earlyGC(:,15);

%% plots
allearlyGC=[];

figure('name', 'AUC greater than 5 in first 1 sec')
hold on
axis off

for xROI= 1:size(earlyGC,1)
    tempY = earlyGC{xROI,8};
    if length(tempY)>590
        tempY=tempY(1:592);
        grey = [0.8,0.8,0.8];
        plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.1);
        allearlyGC= horzcat(allearlyGC, tempY);
    end
end
plot([5 6],[-2 -2], 'k','LineWidth', 2)
plot([5 5],[-1 20], 'k--','LineWidth', 0.5)
meanEarlyTrace = mean(allearlyGC,2);
plot(TimeX(1:stimwindow), meanEarlyTrace(1:stimwindow), 'k', 'LineWidth',1)

allLateGC=[];
figure('name', 'late- responding ROIs not in early group')
hold on
axis off
for xROI= 1:size(lateGC,1)
    tempY = lateGC{xROI,8};
    if length(tempY)>590
        tempY=tempY(1:592);
    grey = [0.8,0.8,0.8];
    plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.1);
    allLateGC= horzcat(allLateGC, tempY);
    end
end
plot([5 6],[-2 -2], 'k','LineWidth', 2)
plot([5 5],[-1 20], 'k--','LineWidth', 0.5)
meanLateTrace = mean(allLateGC,2);
plot(TimeX(1:stimwindow), meanLateTrace(1:stimwindow), 'k', 'LineWidth',1)

% lateExport=allLateGC(1:355,:);
% lateExport2=TimeX(1:355);
% lateExport2=[lateExport2,lateExport];

% control plots
allcontrolGC=[];
figure('name', 'control with AUC greater than 5 in sec 30-31')
hold on
axis off

for xROI= 1:size(controlGC,1)
    tempY = controlGC{xROI,8};
    if length(tempY)>590
        tempY=tempY(1:592);
        grey = [0.8,0.8,0.8];
        plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.1);
        allcontrolGC= horzcat(allcontrolGC, tempY);
    end
end
plot([5 6],[-2 -2], 'k','LineWidth', 2)
plot([5 5],[-1 20], 'k--','LineWidth', 0.5)
meanEarlyTrace = mean(allcontrolGC,2);
plot(TimeX(1:stimwindow), meanEarlyTrace(1:stimwindow), 'k', 'LineWidth',1)


%%
figure('name', 'Mean Late and early traces')
hold on
axis off

plot(TimeX(1:stimwindow), smooth(meanLateTrace(1:stimwindow),3), 'g', 'LineWidth',1)
plot(TimeX(1:stimwindow), smooth(meanEarlyTrace(1:stimwindow),3), 'b', 'LineWidth',1)
% %RCaMP
RespRCmeanTrace = mean(Resp_RCaMP_traces,2);
plot(TimeX(1:stimwindow), smooth(RespRCmeanTrace(1:stimwindow),3), 'r', 'LineWidth',1)
plot([5 6],[-0.5 -0.5], 'k','LineWidth', 2)
plot([5 5],[-0.5 1], 'k--','LineWidth', 0.5)
plot([-1 -1],[0 0.5], 'k','LineWidth', 1)


figure('name', 'Mean Early traces')
hold on
axis off


% %RCaMP
RespRCmeanTrace = mean(Resp_RCaMP_traces,2);
plot(TimeX(1:stimwindow), RespRCmeanTrace(1:stimwindow), 'r', 'LineWidth',1)
plot([5 6],[-0.5 -0.5], 'k','LineWidth', 2)
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
cell2csv('earlyGC_byAUCshort.csv',earlyGCnames2)

% %% early based on time shifts
% % make new names
% for iROI=1:size(early,1)
%     early{iROI,14}=strcat(early{iROI,1},'_',early{iROI,7});
% end
% earlyROIs=unique(early(:,14));
% 
% % find the index of ROIs at are in the early ROI group
% ind1=find(ismember(RespondingROIs(:,15),earlyROIs(:,1)));
% earlyROIs_traces=RespondingROIs(ind1,:);
% 
% ind2=find(ismember(responders(:,26),earlyROIs(:,1)));
% earlyROIs_peaks=responders(ind2,:);
% 
% %remove duplicate traces
% wd=earlyROIs_traces;
% [~,idx]=unique(wd(:,15));
% earlyROIs_traces=wd(idx,:);
% 
% %separate out traces of neurons and astrocytes that are unique
% % [b, m, n]=unique(earlyROIs_traces(:,15));
% % dupindx=find(diff(sort(m))>1)+1 ;
% % dupvals=earlyROIs_traces(dupindx,15);
% 
% %save names
% earlyGCnames=earlyROIs_traces(:,15);
% earlyGCnames2=vertcat('names',earlyGCnames);
% cd('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results');
% cell2csv('earlyGC_byTimeDiff.csv',earlyGCnames2)
% 
% %%
% 
% figure('name', 'Early Extra traces- unshifted 20')
% hold on
% axis off
% ylim([-2 10])
% 
% for xROI= 1:size(earlyROIs_traces,1)
%     tempY = earlyROIs_traces{xROI,8};
%     if length(tempY)>590
%         tempY=tempY(1:592);
%         if earlyROIs_traces{xROI,17}>5
%             grey = [0.8,0.8,0.8];
%             plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.1);
%             earlyExtratraces= horzcat(earlyExtratraces, tempY);
%         end
%     end
% end
% plot([5 6],[-2 -2], 'k','LineWidth', 2)
% plot([5 5],[-1 20], 'k--','LineWidth', 0.5)
% meanEarlyTrace = mean(earlyExtratraces,2);
% plot(TimeX(1:stimwindow), meanEarlyTrace(1:stimwindow), 'k', 'LineWidth',1)
% 
% 
% 
% 
% 
% %%
% 
% figure('name', 'responding RCaMP + Early traces GCaMP')
% hold on
% axis off
% ylim([-2 25])
% for xROI= 1:size(earlyROIs_traces,1)
%     tempY = earlyROIs_traces{xROI,8};
%     grey = [0.8,0.8,0.8];
%     plot(TimeX(1:stimwindow),(tempY(1:stimwindow)')+15,'Color',grey,'LineWidth',0.1);
% end
% for xROI= 1:size(Resp_RCaMP_traces,2)
%     tempY = Resp_RCaMP_traces(:,xROI);
%     grey = [0.8,0.8,0.8];
%     plot(TimeX(1:stimwindow),tempY(1:stimwindow)','Color',grey,'LineWidth',0.1);
% end
% plot([5 5],[-1 100], 'k--','LineWidth', 0.5)
% plot([5 6],[-2 -2], 'k','LineWidth', 2)
% %%
% % early_early=earlyROIs_peaks{:,6}<5;
% for iR=1:size(earlyROIs_peaks,1)
%     test=find(earlyROIs_peaks{iR,6}<5);
% end
