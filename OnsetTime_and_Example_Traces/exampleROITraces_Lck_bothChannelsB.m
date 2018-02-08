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

load('D:\Data\GCaMP_RCaMP\Revision\cyto_GCaMP\FilesforMatlab\Traces_allMice_cyto_nostim_vs_longstim_12_2017.mat');
All_traces(:,7)=[];
All_traces(:,14)=[];
All_traces(:,14)=[];
All_traces(:,14)=[];

Shortstim=All_traces;

% Extract info for the required mice
for iROI=1:length(Shortstim)
    str4(iROI)= ~isempty(strfind(Shortstim{iROI,5},'WT_LR4'));
end

WT_LR4=Shortstim(str4',:);

TimeX(1:nframes) = (1:nframes)/WT_LR4{1,13};
BL_time=WT_LR4{1, 8}/WT_LR4{1,13};
baselineCorrectedTime=TimeX-BL_time;

%% Control

for iROI=1:length(WT_LR4)
    spt3B_str(iROI)= ~isempty(strfind(WT_LR4{iROI,4},'spot3'));
end
WT_LR4_spot3=WT_LR4(spt3B_str',:);

for iROI=1:length(WT_LR4_spot3)
    rcamp_str(iROI)= ~isempty(strfind(WT_LR4_spot3{iROI,3},'RCaMP'));
    gcamp_str(iROI)= ~isempty(strfind(WT_LR4_spot3{iROI,3},'GCaMP'));
end
WT_LR4_spot3_RC=WT_LR4_spot3(rcamp_str',:);
WT_LR4_spot3_GC=WT_LR4_spot3(gcamp_str',:);

% trial 1
for iROI=1:length(WT_LR4_spot3_RC)
    trial_str(iROI)= ~isempty(strfind(WT_LR4_spot3_RC{iROI,2},'trial01'));
end
trial8=WT_LR4_spot3_RC(trial_str',:);

figure('name', 'responding RCaMP ROIs- WT_LR4, Spot3 Trial1')
hold on
axis off
for ii=1:length(trial8)
    %Responding(ii)=~isempty(find((trial8{ii,14}>0 && trial8{ii,14}<=Fast_OnsetWindow),1));
    tempY1=smooth(trial8{ii,9},3);
    tempY1=tempY1(1:nframes);
    %if Responding(ii)
    plot(baselineCorrectedTime(1:nframes),tempY1'+(3*(ii-1)),'r')%'LineWidth',1);
    %end
    
end

for iROI=1:length(WT_LR4_spot3_GC)
    trial_str2(iROI)= ~isempty(strfind(WT_LR4_spot3_GC{iROI,2},'trial01'));
end
trial1_GC=WT_LR4_spot3_GC(trial_str2',:);

figure('name', 'responding GCaMP ROIs- WT_LR4, Spot3 Trial1')
hold on
axis off
for ii=1:length(trial1_GC)
    %Responding(ii)=~isempty(find((trial8{ii,14}>0 && trial8{ii,14}<=Fast_OnsetWindow),1));
    tempY1=smooth(trial1_GC{ii,9},3);
    tempY1=tempY1(1:nframes);
    %if Responding(ii)
    plot(baselineCorrectedTime(1:nframes),tempY1'+(3*(ii-1)),'g')%'LineWidth',1);
    %end
    
end


%% stim neurons
Maps=zeros(128,128);
figure('name','Lck RCaMP examples')
hold on
axis off
ROI=[29,59];%,71];
for ii=1:size(ROI,2)
    ROInum=ROI(ii);
    tempY1=smooth(trial8{ROInum,9},3);
    tempY1=tempY1(1:nframes);
    plot(baselineCorrectedTime(1:nframes),tempY1'+(5*(ii-1)),'r')%'LineWidth',1);
    
    % mask
    mask=trial8{ROInum,10};
    if islogical(mask)
        Mask2=double(mask);
    else
        Image1=zeros(127,128);
        Image1(mask)=1;
        Image2=zeros(1,128);
        Mask2=[Image1;Image2];
    end
    Maps=Maps+Mask2;
end
rectangle('Position', [0 -0.3 8 8])
plot([-10.1 -10.1],[0 1], 'k','LineWidth', 1)


% ROI mask plot
DendMask=im2bw(Maps);
Dend_B=bwboundaries(DendMask);

figure();
imshow(zeros(128,128)); hold on
for k=1:length(Dend_B)
    border=Dend_B{k};
    plot(border(:,2),border(:,1),'m','linewidth',1.5);
end


%%  no stim neurons
Maps=zeros(128,128);
figure('name','Lck RCaMP no stim examples')
hold on
axis off
ROI=[1,17];%,71];
for ii=1:size(ROI,2)
    ROInum=ROI(ii);
    tempY1=smooth(trial8{ROInum,9},3);
    tempY1=tempY1(1:nframes);
    plot(baselineCorrectedTime(1:nframes),tempY1'+(5*(ii-1)),'r')%'LineWidth',1);
    
    % mask
    mask=trial8{ROInum,10};
    if islogical(mask)
        Mask2=double(mask);
    else
        Image1=zeros(127,128);
        Image1(mask)=1;
        Image2=zeros(1,128);
        Mask2=[Image1;Image2];
    end
    Maps=Maps+Mask2;
end
plot([-10.1 -10.1],[0 1], 'k','LineWidth', 1)


% ROI mask plot
DendMask=im2bw(Maps);
Dend_B=bwboundaries(DendMask);

figure();
imshow(zeros(128,128)); hold on
for k=1:length(Dend_B)
    border=Dend_B{k};
    plot(border(:,2),border(:,1),'m','linewidth',1.5);
end


%%  no stim astrocytes
Maps=zeros(128,128);
figure('name','Lck GCaMP no stim examples')
hold on
axis off
ROI=[9,2];%,71];
for ii=1:size(ROI,2)
    ROInum=ROI(ii);
    tempY1=smooth(trial1_GC{ROInum,9},3);
    tempY1=tempY1(1:nframes);
    plot(baselineCorrectedTime(1:nframes),tempY1'+(5*(ii-1)),'r')%'LineWidth',1);
    
    % mask
    mask=trial1_GC{ROInum,10};
    if islogical(mask)
        Mask2=double(mask);
    else
        Image1=zeros(127,128);
        Image1(mask)=1;
        Image2=zeros(1,128);
        Mask2=[Image1;Image2];
    end
    Maps=Maps+Mask2;
end
plot([-10.1 -10.1],[0 1], 'k','LineWidth', 1)


% ROI mask plot
DendMask=im2bw(Maps);
Dend_B=bwboundaries(DendMask);

figure();
imshow(zeros(128,128)); hold on
for k=1:length(Dend_B)
    border=Dend_B{k};
    plot(border(:,2),border(:,1),'g','linewidth',1.5);
end
