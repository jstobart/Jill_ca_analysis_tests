%% LONG STIM (8sec)
clearvars

load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\S&LStim_cGC&RC_traces_01_30_2017.mat');

FrameRate=11.84;

searchString='Stim';

for iROI=1:length(All_traces)
    Stim_str{iROI} = strfind(All_traces{iROI,6},searchString);
    str_idx(iROI)=~isempty(Stim_str{iROI});
    
    %find ROITypes
    N_str= strfind(All_traces{iROI, 1},'N');
    S_str= strfind(All_traces{iROI, 1},'S');
    EF_str= strfind(All_traces{iROI, 1},'E');
    P_str= strfind(All_traces{iROI, 1},'r');
    NP_str= strfind(All_traces{iROI, 1},'np');
    
    if ~isempty(N_str)
        All_traces{iROI,12}='Neuron';
    elseif ~isempty(S_str)
        All_traces{iROI,12}='Soma';
    elseif ~isempty(EF_str)
        All_traces{iROI,12}='Endfeet';
    elseif ~isempty(P_str)
        All_traces{iROI,12}='Process';
    elseif ~isempty(NP_str)
        All_traces{iROI,12}='Neuropil';
    end
    
    % make new names
    All_traces{iROI,13}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,1});
    All_traces{iROI,14}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,2}, '_',All_traces{iROI,1});
    All_traces{iROI,15}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,2}, '_',All_traces{iROI,1},'_',All_traces{iROI,6});
    All_traces{iROI,16}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,2}, '_',All_traces{iROI,6});
    All_traces{iROI,17}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,2});
    
    %if a process ROI is overlapping with a soma or process, exclude it
    %from data
    %All_traces{iROI,11}=0;
    nonOverlapIdx(iROI)=~ischar(All_traces{iROI,11});
    
end


% remove overlapping processes
All_traces = All_traces(nonOverlapIdx',:);

% only trials for this condition
data_traces = All_traces(str_idx',:);

%% Plot all traces for the entire trial

% individual traces in grey, mean in red

nROIs=length(data_traces);

nframes=length(All_traces{1,8});

TimeX(1:nframes) = (1:nframes)/FrameRate;

traces=[];
figure ('name', 'All ROIs- whole trial')
hold on
axis off
for xROI= 1:nROIs
    tempY = data_traces{xROI,8};
    grey = [0.8,0.8,0.8];
    plot(TimeX,tempY,'Color',grey,'LineWidth',0.1);
    traces= horzcat(traces, tempY);
end

plot([5 13],[-2 -2], 'k','LineWidth', 2)
meanTrace = mean(traces,2);
plot(TimeX, meanTrace, 'k', 'LineWidth',1)
% plot([0 0],[18.5 19.5], 'k','LineWidth', 2)



%% Plot only the responding neurons and astrocytes from the same field of view


XLfile = 'E:\Data\Two_Photon_Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\respondingROIs_longstim.xlsx';

[~, ~, data] = xlsread(XLfile); %,'Sheet1'); %reads the scoresheet and saves all data in a cell array
data=data(2:end,:);

%separate out responding trials from responding neurons and astrocytes
RespondingROIs=[];
for xROI=1:length(data)
    currentROI=data{xROI, 24};
    for iROI = 1:length(data_traces)
        ROIIndex=strcmp(data_traces{iROI,14},currentROI);
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
    if ~isempty(RC_str)
        Resp_RCaMP_traces= horzcat(Resp_RCaMP_traces, tempY);
    elseif ~isempty(GC_str)
        Resp_GCaMP_traces= horzcat(Resp_GCaMP_traces, tempY);
    end
end


%find ROITypes
for xROI= 1:length(RespondingROIs)
    N_str= strfind(RespondingROIs{xROI, 12},'Neuron');
    S_str= strfind(RespondingROIs{xROI, 12},'Soma');
    EF_str= strfind(RespondingROIs{xROI, 12},'Endfeet');
    P_str= strfind(RespondingROIs{xROI, 12},'Process');
    NP_str= strfind(RespondingROIs{xROI, 12},'Neuropil');
    
    tempY= RespondingROIs{xROI,8};
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
figure ('name', 'All RCaMP responding ROIs')
hold on
axis off
for xROI= 1:size(Resp_RCaMP_traces,2)
    tempY = Resp_RCaMP_traces(:,xROI);
    grey = [0.8,0.8,0.8];
    plot(TimeX,tempY,'Color',grey,'LineWidth',0.01);
end

plot([5 13],[-2 -2], 'k','LineWidth', 2)
plot(TimeX, RespRCmeanTrace, 'k', 'LineWidth',1)

figure ('name', 'All GCaMP responding ROIs')
hold on
axis off
for xROI= 1:size(Resp_GCaMP_traces,2)
    tempY = Resp_GCaMP_traces(:,xROI);
    grey = [0.8,0.8,0.8];
    plot(TimeX,tempY,'Color',grey,'LineWidth',0.01);
end

plot([5 13],[-2 -2], 'k','LineWidth', 2)
plot(TimeX, RespGCmeanTrace, 'k', 'LineWidth',1)

stimwindow=round(FrameRate*30);
figure ('name', 'All RCaMP responding ROIs- stim window')
hold on
axis off
for xROI= 1:size(Resp_RCaMP_traces,2)
    tempY = Resp_RCaMP_traces(:,xROI);
    grey = [0.8,0.8,0.8];
    plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.01);
end

plot([5 13],[-2 -2], 'k','LineWidth', 2)
plot(TimeX(1:stimwindow), RespRCmeanTrace(1:stimwindow), 'k', 'LineWidth',1)

figure ('name', 'All GCaMP responding ROIs- stim window')
hold on
axis off
for xROI= 1:size(Resp_GCaMP_traces,2)
    tempY = Resp_GCaMP_traces(:,xROI);
    grey = [0.8,0.8,0.8];
    plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.01);
end

plot([5 13],[-2 -2], 'k','LineWidth', 2)
plot(TimeX(1:stimwindow), RespGCmeanTrace(1:stimwindow), 'k', 'LineWidth',1)

%% Mean RCaMP vs GCaMP
% GCaMP vs RCaMP
figure ('name', 'GCaMP vs RCaMP responding ROIs- stim window')
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


figure ('name', 'GCaMP vs RCaMP responding ROIs- mean trace only')
hold on
axis off
plot(TimeX(1:stimwindow),smooth(RespRCmeanTrace(1:stimwindow)',5),'r','LineWidth',1.5);
plot(TimeX(1:stimwindow),smooth(RespGCmeanTrace(1:stimwindow)',5),'g','LineWidth', 1.5);
plot([5 5],[-1 2], 'k--','LineWidth', 1)
plot([5 13],[0 0], 'k','LineWidth', 3)
legend('RCaMP','GCaMP')



%% ROITypes
figure ('name', 'All responding ROItype')
subplot(1,5,1)
hold on
axis off
for xROI= 1:size(Resp_neuron_traces,2)
    tempY = Resp_neuron_traces(:,xROI);
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
for xROI= 1:size(Resp_neuropil_traces,2)
    tempY = Resp_neuropil_traces(:,xROI);
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


subplot(1,5,3)
hold on
axis off
for xROI= 1:size(Resp_EF_traces,2)
    tempY = Resp_EF_traces(:,xROI);
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


subplot(1,5,4)
hold on
axis off
for xROI= 1:size(Resp_somata_traces,2)
    tempY = Resp_somata_traces(:,xROI);
    grey = [0.8,0.8,0.8];
    plot(TimeX(1:stimwindow),tempY(1:stimwindow),'Color',grey,'LineWidth',0.01);
end
plot([5 13],[-2 -2], 'k','LineWidth', 2)
plot(TimeX(1:stimwindow), RespSmeanTrace(1:stimwindow), 'k', 'LineWidth',1)
plot([5 5],[-1 10], 'k--','LineWidth', 0.5)
x1 = 0;
y1 = 5;
txt1 = 'Somata';
text(x1,y1,txt1,'HorizontalAlignment','right')


subplot(1,5,5)
hold on
axis off
for xROI= 1:size(Resp_processes_traces,2)
    tempY = Resp_processes_traces(:,xROI);
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

%% mean ROItype traces

figure ('name', 'ROITypes responding-means- stim window')
hold on
axis off
plot([5 13],[-1 -1], 'k','LineWidth', 2)
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
txt3 = 'Somata';
txt4 = 'Endfeet';
txt5 = 'Processes';
text(x1,y1,txt1,'HorizontalAlignment','right')
text(x1,y2,txt2,'HorizontalAlignment','right')
text(x1,y3,txt3,'HorizontalAlignment','right')
text(x1,y4,txt4,'HorizontalAlignment','right')
text(x1,y5,txt5,'HorizontalAlignment','right')


%ROITypes
figure ('name', 'responding ROITypes- mean traces only- stim window')
hold on
axis off
plot(TimeX(1:stimwindow),smooth(RespNmeanTrace(1:stimwindow)',5),'r','LineWidth', 1.5);
plot(TimeX(1:stimwindow),smooth(RespNPmeanTrace(1:stimwindow)',5),'b','LineWidth', 1.5);
plot(TimeX(1:stimwindow),smooth(RespSmeanTrace(1:stimwindow)',5),'g','LineWidth', 1.5);
plot(TimeX(1:stimwindow),smooth(RespEFmeanTrace(1:stimwindow)',5),'k','LineWidth', 1.5);
plot(TimeX(1:stimwindow),smooth(RespPmeanTrace(1:stimwindow)',5),'m','LineWidth', 1.5);
plot([5 5],[-1 2], 'k--','LineWidth', 1)
plot([5 13],[0 0], 'k','LineWidth', 3)
legend('Neuron','Neuropil','Somata','Endfeet','Processes')



%% Raster Plot with event frequency histogram

% each peak time is a line
% lines are coloured and sorted based on ROI type

    figure('name', 'all ROIs')% ROIType{iType})
    hold on
    axis off
for iType=1:5
    ROIType={'Neuron','Neuropil','Endfoot','Soma','Process'};
    
    %find ROITypes
    for xROI= 1:size(data,1)
        ROI_str(xROI)= strcmp(data{xROI, 23},ROIType{iType});
    end
    
    % only unique ROIs of a given type
    uniqueROIs=unique(data(ROI_str,24));
    
    colours={[0.6,0,0],[0,0,0.6],[0.2,0.4,0],[0.8,0.2,0],[0.8,0,0.8]};
    % plot raster plot of each ROI and trial

    %set(gca,'TickDir','out') % draw the tick marks on the outside
    %set(gca,'YTick', []) % don't draw y-axis ticks
    %set(gca,'PlotBoxAspectRatio',[1 0.05 1]) % short and wide
    %set(gca,'Color',get(gcf,'Color')) % match figure background
    %set(gca,'YColor',get(gcf,'Color')) % hide the y axis
    for iROI= 1:length(uniqueROIs)
        CurrentROI=uniqueROIs(iROI);
        for xROI=1:size(data,1)
            ROI_idx(xROI)=strcmp(data{xROI,24},CurrentROI);
        end
        peakTimes=cell2mat(data(ROI_idx,6));
        if size(peakTimes,1) > size(peakTimes,2)
            peakTimes=peakTimes';
        end
        plot([peakTimes;peakTimes],[ones(size(peakTimes))+(iROI-1);zeros(size(peakTimes))+(iROI-1)],'Color',colours{1,iType},'LineWidth', 2)
    end
    %plot([0 20],[-1 -1], 'k--','LineWidth', 1)
    %axis([0 20 -1 2])
end

%% Histogram of peak time differences between neurons and astrocytes

%% example traces aligning astrocyte and neuronal peaks


