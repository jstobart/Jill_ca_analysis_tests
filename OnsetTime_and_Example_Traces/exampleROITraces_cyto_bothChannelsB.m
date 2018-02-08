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

BL_time=5.5;

%% load control trace data

load('D:\Data\GCaMP_RCaMP\Revision\cytoGCaMP\FilesforMatlab\Traces_allMice_cyto_nostim_vs_longstim_12_2017.mat');
% All_traces(:,7)=[];
% All_traces(:,14)=[];
% All_traces(:,14)=[];
% All_traces(:,14)=[];

Shortstim=All_traces;

% Extract info for the required mice
for iROI=1:length(Shortstim)
    str4(iROI)= ~isempty(strfind(Shortstim{iROI,5},'RG14'));
end

RG14=Shortstim(str4',:);

TimeX(1:nframes) = (1:nframes)/RG14{1,14};
baselineCorrectedTime=TimeX-BL_time;

%% Control

for iROI=1:length(RG14)
    spt3B_str(iROI)= ~isempty(strfind(RG14{iROI,4},'16_02_24_spot1'));
end
RG14_spot1=RG14(spt3B_str',:);

for iROI=1:length(RG14_spot1)
    rcamp_str(iROI)= ~isempty(strfind(RG14_spot1{iROI,3},'RCaMP'));
    gcamp_str(iROI)= ~isempty(strfind(RG14_spot1{iROI,3},'GCaMP'));
end
RG14_spot1_RC=RG14_spot1(rcamp_str',:);
RG14_spot1_GC=RG14_spot1(gcamp_str',:);

% trial 1
for iROI=1:length(RG14_spot1_RC)
    trial_str(iROI)= ~isempty(strfind(RG14_spot1_RC{iROI,2},'trial04'));
end
trial8=RG14_spot1_RC(trial_str',:);

figure('name', 'responding RCaMP ROIs- RG14, Spot1 Trial4')
hold on
axis off
for ii=1:length(trial8)
    %Responding(ii)=~isempty(find((trial8{ii,14}>0 && trial8{ii,14}<=Fast_OnsetWindow),1));
    tempY1=smooth(trial8{ii,10},3);
    tempY1=tempY1(1:nframes);
    %if Responding(ii)
    plot(baselineCorrectedTime(1:nframes),tempY1'+(3*(ii-1)),'r')%'LineWidth',1);
    %end
    
end
rectangle('Position', [0 -0.3 8 140])

for iROI=1:length(RG14_spot1_GC)
    trial_str2(iROI)= ~isempty(strfind(RG14_spot1_GC{iROI,2},'trial01'));
end
trial1_GC=RG14_spot1_GC(trial_str2',:);

figure('name', 'responding GCaMP ROIs- RG14, Spot3 Trial4')
hold on
axis off
for ii=1:length(trial1_GC)
    %Responding(ii)=~isempty(find((trial8{ii,14}>0 && trial8{ii,14}<=Fast_OnsetWindow),1));
    tempY1=smooth(trial1_GC{ii,10},3);
    tempY1=tempY1(1:nframes);
    %if Responding(ii)
    plot(baselineCorrectedTime(1:nframes),tempY1'+(3*(ii-1)),'g')%'LineWidth',1);
    %end
    
end
rectangle('Position', [0 -0.3 8 100])

%% stim neurons
Maps=zeros(128,128);
figure('name','Lck RCaMP examples')
hold on
axis off
ROI=[40];
for ii=1:size(ROI,2)
    ROInum=ROI(ii);
    tempY1=smooth(trial8{ROInum,10},5);
    tempY1=tempY1(1:nframes);
    plot(baselineCorrectedTime(1:nframes),tempY1'+(5*(ii-1)),'r')%'LineWidth',1);
    
    % mask
    mask=trial8{ROInum,11};
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
rectangle('Position', [0 -0.3 8 12])
plot([-5.1 -5.1],[0 2], 'k','LineWidth', 1)


% ROI mask plot
DendMask=im2bw(Maps);
Dend_B=bwboundaries(DendMask);

figure();
imshow(zeros(128,128)); hold on
for k=1:length(Dend_B)
    border=Dend_B{k};
    plot(border(:,2),border(:,1),'m','linewidth',1.5);
end


%%  stim astrocytes
Maps=zeros(128,128);
figure('name','Lck GCaMP stim examples')
hold on
axis off
ROI=[10,14];%,71];
for ii=1:size(ROI,2)
    ROInum=ROI(ii);
    tempY1=smooth(trial1_GC{ROInum,10},5);
    tempY1=tempY1(1:nframes);
    plot(baselineCorrectedTime(1:nframes),tempY1'+(5*(ii-1)),'g')%'LineWidth',1);
    
    % mask
    mask=trial1_GC{ROInum,11};
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
rectangle('Position', [0 -0.3 8 12])
plot([-5.1 -5.1],[0 2], 'k','LineWidth', 1)


% ROI mask plot
DendMask=im2bw(Maps);
Dend_B=bwboundaries(DendMask);

figure();
imshow(zeros(128,128)); hold on
for k=1:length(Dend_B)
    border=Dend_B{k};
    plot(border(:,2),border(:,1),'g','linewidth',1.5);
end

%%  no stim astrocytes
Maps=zeros(128,128);
figure('name','Lck GCaMP nostim examples')
hold on
axis off
ROI=[1,7];%,71];
for ii=1:size(ROI,2)
    ROInum=ROI(ii);
    tempY1=smooth(trial1_GC{ROInum,10},5);
    tempY1=tempY1(1:nframes);
    plot(baselineCorrectedTime(1:nframes),tempY1'+(4*(ii-1)),'g')%'LineWidth',1);
    
    % mask
    mask=trial1_GC{ROInum,11};
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
plot([-5.1 -5.1],[0 2], 'k','LineWidth', 1)


% ROI mask plot
DendMask=im2bw(Maps);
Dend_B=bwboundaries(DendMask);

figure();
imshow(zeros(128,128)); hold on
for k=1:length(Dend_B)
    border=Dend_B{k};
    plot(border(:,2),border(:,1),'g','linewidth',1.5);
end

%%
clearvars trial8 trial1_GC trial_str trial_str2

for iROI=1:length(RG14_spot1_RC)
    trial_str(iROI)= ~isempty(strfind(RG14_spot1_RC{iROI,2},'trial05'));
end
trial8=RG14_spot1_RC(trial_str',:);


figure('name', 'responding RCaMP ROIs- RG14, Spot1 Trial4')
hold on
axis off
for ii=1:length(trial8)
    %Responding(ii)=~isempty(find((trial8{ii,14}>0 && trial8{ii,14}<=Fast_OnsetWindow),1));
    tempY1=smooth(trial8{ii,10},3);
    tempY1=tempY1(1:nframes);
    %if Responding(ii)
    plot(baselineCorrectedTime(1:nframes),tempY1'+(3*(ii-1)),'r')%'LineWidth',1);
    %end
    
end
rectangle('Position', [0 -0.3 8 140])

for iROI=1:length(RG14_spot1_GC)
    trial_str2(iROI)= ~isempty(strfind(RG14_spot1_GC{iROI,2},'trial05'));
end
trial1_GC=RG14_spot1_GC(trial_str2',:);

figure('name', 'responding GCaMP ROIs- RG14, Spot3 Trial4')
hold on
axis off
for ii=1:length(trial1_GC)
    %Responding(ii)=~isempty(find((trial8{ii,14}>0 && trial8{ii,14}<=Fast_OnsetWindow),1));
    tempY1=smooth(trial1_GC{ii,10},3);
    tempY1=tempY1(1:nframes);
    %if Responding(ii)
    plot(baselineCorrectedTime(1:nframes),tempY1'+(3*(ii-1)),'g')%'LineWidth',1);
    %end
    
end
rectangle('Position', [0 -0.3 8 100])

%%  no stim neurons
Maps=zeros(128,128);
figure('name','Lck RCaMP no stim examples')
hold on
axis off
ROI=[3,11];%,71];
for ii=1:size(ROI,2)
    ROInum=ROI(ii);
    tempY1=smooth(trial8{ROInum,10},5);
    tempY1=tempY1(1:nframes);
    plot(baselineCorrectedTime(1:nframes),tempY1'+(3*(ii-1)),'r')%'LineWidth',1);
    
    % mask
    mask=trial8{ROInum,11};
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
plot([-5.1 -5.1],[0 2], 'k','LineWidth', 1)


% ROI mask plot
DendMask=im2bw(Maps);
Dend_B=bwboundaries(DendMask);

figure();
imshow(zeros(128,128)); hold on
for k=1:length(Dend_B)
    border=Dend_B{k};
    plot(border(:,2),border(:,1),'m','linewidth',1.5);
end

%%  stim 1 neurons
Maps=zeros(128,128);
figure('name','Lck RCaMP stim examples')
hold on
axis off
ROI=[42];%,71];
for ii=1:size(ROI,2)
    ROInum=ROI(ii);
    tempY1=smooth(trial8{ROInum,10},5);
    tempY1=tempY1(1:nframes);
    plot(baselineCorrectedTime(1:nframes),tempY1'+(3*(ii-1)),'r')%'LineWidth',1);
    
    % mask
    mask=trial8{ROInum,11};
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
plot([-5.1 -5.1],[0 2], 'k','LineWidth', 1)
rectangle('Position', [0 -0.3 8 5])

% ROI mask plot
DendMask=im2bw(Maps);
Dend_B=bwboundaries(DendMask);

figure();
imshow(zeros(128,128)); hold on
for k=1:length(Dend_B)
    border=Dend_B{k};
    plot(border(:,2),border(:,1),'m','linewidth',1.5);
end
