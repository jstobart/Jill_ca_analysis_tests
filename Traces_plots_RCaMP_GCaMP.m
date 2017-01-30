%% Sort out conditions- no stim, long stim, short stim

numROIs=length(All_traces);

searchString='Stim';

for xROI=1:numROIs
Stim_str{xROI} = strfind(All_traces{xROI,6},searchString);
str_idx(xROI)=~isempty(Stim_str{xROI});
end

data_traces = All_traces(str_idx',:);

%% Plot all traces for the entire trial

% individual traces in grey, mean in red

% plot each ROItype separately (neurons, neuropil, AC processes, somata,
% endfeet)

nROIs=length(data_traces);

nframes=length(All_traces{1,8});

FrameRate = 11.84;
TimeX(1:nframes) = (1:nframes)/FrameRate;

traces=[];
figure ('name', 'All ROIs- whole trial')
hold on
axis off
for xROI= 1:nROIs
    tempY = data_traces{xROI,8};
    grey = [0.7,0.7,0.7];
    plot(TimeX,tempY,'Color',grey,'LineWidth',0.5);
    traces= horzcat(traces, tempY);
end

plot([5 13],[-2 -2], 'k','LineWidth', 2)
meanTrace = mean(traces,2);
plot(TimeX, meanTrace, 'k', 'LineWidth',1)
% plot([0 0],[18.5 19.5], 'k','LineWidth', 2)

%% plot neurons and astrocytes separately
gcamp_traces=[];
rcamp_traces=[];

figure ('name', 'GCaMP vs RCaMP- whole trial')
hold on
axis off
%find GCaMP and RCaMP
for xROI= 1:nROIs
    GC_str(xROI,1)= strcmp(data_traces{xROI, 3},'GCaMP');
    tempY= data_traces{xROI,8};
   if GC_str(xROI,1)
       plot(TimeX, tempY, 'g', 'LineWidth', 0.5)
       gcamp_traces= horzcat(gcamp_traces, tempY);
   else
       plot(TimeX, tempY, 'r', 'LineWidth', 0.5)
       rcamp_traces= horzcat(rcamp_traces, tempY);
   end
       
end
plot([5 13],[-2 -2], 'k','LineWidth', 2)

%gcamp mean and sd
GCmeanTrace = mean(gcamp_traces,2);
GCSDTrace = std(gcamp_traces');
plot(TimeX, GCmeanTrace, 'k', 'LineWidth',1)

%rcamp mean and sd
RCmeanTrace = mean(rcamp_traces,2);
RCSDTrace = std(rcamp_traces');
plot(TimeX, RCmeanTrace, 'k', 'LineWidth',1)

figure ('name', 'GCaMP vs RCaMP 2- whole trial')
hold on
axis off
plot([5 13],[-2 -2], 'k','LineWidth', 2)
shadedErrorBar(TimeX,RCmeanTrace',RCSDTrace,'r');    
shadedErrorBar(TimeX,(GCmeanTrace' + 2.5),GCSDTrace,'g');
plot([5 5],[-2 5], 'k--','LineWidth', 0.5)
x1 = -2;
y1 = 0;
y2 = 2.5;

txt1 = 'Neurons';
txt2 = 'Astrocytes';

text(x1,y1,txt1,'HorizontalAlignment','right')
text(x1,y2,txt2,'HorizontalAlignment','right')

%% Plot individual ROI types

neuron_traces=[];
neuropil_traces=[];
somata_traces=[];
EF_traces=[];
processes_traces=[];

%find ROITypes
for xROI= 1:nROIs
    N_str= strfind(data_traces{xROI, 1},'N');
    S_str= strfind(data_traces{xROI, 1},'S');
    EF_str= strfind(data_traces{xROI, 1},'E');
    P_str= strfind(data_traces{xROI, 1},'r');
    NP_str= strfind(data_traces{xROI, 1},'np');
    
    tempY= data_traces{xROI,8};
    if ~isempty(N_str)
        neuron_traces= horzcat( neuron_traces, tempY);
        data_traces{xROI,11}='Neuron';
    elseif ~isempty(S_str)
        somata_traces= horzcat(somata_traces, tempY);
        data_traces{xROI,11}='Soma';
    elseif ~isempty(EF_str)
        EF_traces= horzcat(EF_traces, tempY);
        data_traces{xROI,11}='Endfeet';
    elseif ~isempty(P_str)
        processes_traces= horzcat(processes_traces, tempY);
        data_traces{xROI,11}='Process';
    elseif ~isempty(NP_str)
        neuropil_traces= horzcat(neuropil_traces, tempY);
        data_traces{xROI,11}='Neuropil';
    end
end

%neuron mean and sd
NmeanTrace = mean(neuron_traces,2);
NSDTrace = std(neuron_traces');

NPmeanTrace = mean(neuropil_traces,2);
NPSDTrace = std(neuropil_traces');

%astrocytes mean and sd
SmeanTrace = mean(somata_traces,2);
SSDTrace = std(somata_traces');

EmeanTrace = mean(EF_traces,2);
ESDTrace = std(EF_traces');

PmeanTrace = mean(processes_traces,2);
PSDTrace = std(processes_traces');


figure ('name', 'ROITypes whole trial')
hold on
axis off
plot([5 13],[-2 -2], 'k','LineWidth', 2)
shadedErrorBar(TimeX,NmeanTrace',NSDTrace,'r');    
shadedErrorBar(TimeX,(NPmeanTrace' + 3),NPSDTrace,'b');
shadedErrorBar(TimeX,(SmeanTrace' + 6),SSDTrace,'g');
shadedErrorBar(TimeX,(EmeanTrace' + 9),ESDTrace,'c');
shadedErrorBar(TimeX,(PmeanTrace' + 12),PSDTrace,'m');
plot([5 5],[-2 15], 'k--','LineWidth', 0.5)
x1 = -2;
y1 = 0;
y2 = 3;
y3 = 6;
y4 = 9;
y5 = 12;

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


%% Plot only the window around the stimulus- 20 sec after onset

stimwindow= FrameRate*25;

% GCaMP vs RCaMP
figure ('name', 'GCaMP vs RCaMP 2- stim window')
hold on
axis off
plot([5 13],[-1 -1], 'k','LineWidth', 2)
shadedErrorBar(TimeX(1:stimwindow),RCmeanTrace(1:stimwindow)',RCSDTrace(1:stimwindow),'r');    
shadedErrorBar(TimeX(1:stimwindow),(GCmeanTrace(1:stimwindow)' + 2.5),GCSDTrace(1:stimwindow),'g');
plot([5 5],[-1 5], 'k--','LineWidth', 0.5)
x1 = -0.5;
y1 = 0;
y2 = 2.5;

txt1 = 'Neurons';
txt2 = 'Astrocytes';

text(x1,y1,txt1,'HorizontalAlignment','right')
text(x1,y2,txt2,'HorizontalAlignment','right')

%ROITypes
figure ('name', 'ROITypes- stim window')
hold on
axis off
plot([5 13],[-1 -1], 'k','LineWidth', 2)
shadedErrorBar(TimeX(1:stimwindow),NmeanTrace(1:stimwindow)',NSDTrace(1:stimwindow),'r');    
shadedErrorBar(TimeX(1:stimwindow),(NPmeanTrace(1:stimwindow)' + 3),NPSDTrace(1:stimwindow),'b');
shadedErrorBar(TimeX(1:stimwindow),(SmeanTrace(1:stimwindow)' + 6),SSDTrace(1:stimwindow),'g');
shadedErrorBar(TimeX(1:stimwindow),(EmeanTrace(1:stimwindow)' + 9),ESDTrace(1:stimwindow),'y');
shadedErrorBar(TimeX(1:stimwindow),(PmeanTrace(1:stimwindow)' + 12),PSDTrace(1:stimwindow),'m');
plot([5 5],[-1 15], 'k--','LineWidth', 0.5)
x1 = -0.5;
y1 = 0;
y2 = 3;
y3 = 6;
y4 = 9;
y5 = 12;

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


%% Plot only the responding neurons and astrocytes from the same field of view


%% Plot only responding neurons and responding astrocytes from the same field of view


