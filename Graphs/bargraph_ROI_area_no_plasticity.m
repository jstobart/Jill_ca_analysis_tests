%% Graphs of AP ROI Area- Not for Plasticity Study
% Version 3
%  16.12.2013  JS

%open data files
DirName = 'E:\Data\Two_Photon_Data\GCaMP6\Whisker_Stim\Analysis_Results\activity_ROIs_2013_12_16';
cd(DirName)
flist=dir(fullfile(DirName, '*.mat')); % will give you list of specifed files
for ii=1:length(flist)
     load(flist(ii).name) % load files in workspAPe
end

%% Cell of Average signal Frequencies for eAPh ROI type (with both layers)

RowInfo = ['NS_AP';'S_AP '];
Area_Mean =  cellstr(RowInfo);
Area_SD =  cellstr(RowInfo);
Area_groups = cellstr(RowInfo);
L1_A_Mean =  cellstr(RowInfo);
L1_A_SD =  cellstr(RowInfo);
L2_A_Mean =  cellstr(RowInfo);
L2_A_SD =  cellstr(RowInfo);
L1_A_groups = cellstr(RowInfo);
L2_A_groups = cellstr(RowInfo);


%% determine average AP ROI area (for both layers)
% AP ROIs 
%NO STIM 
%Both Layers
     for iROI = 1:length(Nostim.AP)
           A{iROI} = (Nostim.AP{iROI}{9,2});
           Area =cell2mat(A);  
     end
    Area_Mean{1,2} = nanmean(Area);
    Area_Mean{1,3} = nanstd(Area); 
    Area_groups{1,2} = {Area};

 % Layer 1
     for iROI = 1:length(Layer1.NS_AP)
          L1_A{iROI} = (Layer1.NS_AP{iROI}{9,2});
         L1_Area =cell2mat(L1_A);
    end
    L1_A_groups{1,2} = {L1_Area};
    L1_A_Mean{1,2} = nanmean(L1_Area);
    L1_A_Mean{1,3} = nanstd(L1_Area);
    
    %Layer 2
    for iROI = 1:length(Layer2.NS_AP)
          L2_A{iROI} = (Layer2.NS_AP{iROI}{9,2});
         L2_Area =cell2mat(L2_A);
    end
    L2_A_groups{1,2} = {L2_Area};
    L2_A_Mean{1,2} = nanmean(L2_Area);
    L2_A_Mean{1,3} = nanstd(L2_Area); 
 
   clear 'Area' 'A' 'L1_A' 'L1_Area' 'L2_A' 'L2_Area'
   
   %% STIM 
%Both Layers
     for iROI = 1:length(Stim8s.AP)
           A{iROI} = (Stim8s.AP{iROI}{9,2});
           Area =cell2mat(A);  
     end
    Area_Mean{2,2} = nanmean(Area);
    Area_Mean{2,3} = nanstd(Area); 
    Area_groups{2,2} = {Area};

 % Layer 1
     for iROI = 1:length(Layer1.S_AP)
          L1_A{iROI} = (Layer1.S_AP{iROI}{9,2});
         L1_Area =cell2mat(L1_A);
    end
    L1_A_groups{2,2} = {L1_Area};
    L1_A_Mean{2,2} = nanmean(L1_Area);
    L1_A_Mean{2,3} = nanstd(L1_Area);
    
    %Layer 2
    for iROI = 1:length(Layer2.S_AP)
          L2_A{iROI} = (Layer2.S_AP{iROI}{9,2});
         L2_Area =cell2mat(L2_A);
    end
    L2_A_groups{2,2} = {L2_Area};
    L2_A_Mean{2,2} = nanmean(L2_Area);
    L2_A_Mean{2,3} = nanstd(L2_Area) ;
    
    
   
   
%%  Bar Graphs
%BOTH LAYERS
No_stim = cell2mat(Area_Mean(1,2)');
Stim =cell2mat(Area_Mean(2,2)');
No_stim_SD = cell2mat(Area_Mean(2,3)');
Stim_SD =cell2mat(Area_Mean(2,3)');

Both = [No_stim' Stim'];
Both_SD = [No_stim_SD' Stim_SD'];

%% Layer 1
L1_No_stim = cell2mat(L1_A_Mean(1,2)');
L1_Stim =cell2mat(L1_A_Mean(2,2)');
L1_No_stim_SD = cell2mat(L1_A_Mean(1,3)');
L1_Stim_SD =cell2mat(L1_A_Mean(2,3)');

L1 = [L1_No_stim' L1_Stim'];
L1_SD = [L1_No_stim_SD' L1_Stim_SD'];

%Layer 2
L2_No_stim = cell2mat(L2_A_Mean(1,2)');
L2_Stim =cell2mat(L2_A_Mean(2,2)');
L2_No_stim_SD = cell2mat(L2_A_Mean(1,3)');
L2_Stim_SD =cell2mat(L2_A_Mean(2,3)');

L2 = [L2_No_stim' L2_Stim'];
L2_SD = [L2_No_stim_SD' L2_Stim_SD'];

%Figures
figure ('name','AP ROI Area','numbertitle','off')
subplot (3,1,1)
hold on
barwitherr(Both, Both_SD)
set(gca,'XTickLabel', '')
set(gca,'YLim', [0 200],'YTick',0:50:200)
ylabel('AP ROI area')
title('Both Layers')
legend('NoStim','Stim') 

subplot (3,1,2)
hold on
barwitherr(L1,L1_SD)
set(gca,'XTickLabel', '')
set(gca,'YLim', [0 200],'YTick',0:50:200)
ylabel('AP ROI area')
title('Layer1')
legend('NoStim','Stim') 

subplot (3,1,3)
hold on
barwitherr(L2,L2_SD)
set(gca,'XTickLabel', '')
set(gca,'YLim', [0 200],'YTick',0:50:200)
ylabel('AP ROI area')
title('Layer2')
legend('NoStim','Stim')




%% Statistics- T-test for stim vs. nostim
%Both Layers

%AP ROIs
AP_NS = cell2mat(Area_groups{1,2});
AP_S = cell2mat(Area_groups{2,2});
AP = padcat(AP_NS, AP_S);
AP1= AP(1,:);
AP2= AP(2,:);

[h(1,1), p(1,1)] =ttest(AP1, AP2)


%% Layer 1

%AP ROIs
L1_AP_NS = cell2mat(L1_A_groups{1,2});
L1_AP_S = cell2mat(L1_A_groups{2,2});
L1_AP = padcat(L1_AP_NS, L1_AP_S);
L1_AP1= L1_AP(1,:);
L1_AP2= L1_AP(2,:);

[h(1,2), p(1,2)] =ttest(L1_AP1, L1_AP2)

%% Layer 2

%AP ROIs
L2_AP_NS = cell2mat(L2_A_groups{1,2});
L2_AP_S = cell2mat(L2_A_groups{2,2});
L2_AP = padcat(L2_AP_NS, L2_AP_S);
L2_AP1= L2_AP(1,:);
L2_AP2= L2_AP(2,:);

[h(1,3), p(1,3)] =ttest(L2_AP1, L2_AP2)


%% Create Output folder- groups
file = 'E:\Data\Two_Photon_Data\GCaMP6\Whisker_Stim\Analysis_Results\activity_ROIs_2013_12_16';
outputFolder = fullfile(file, 'Data_Groups');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

cd(outputFolder)

filename1= 'Both_Layers_AP_Area';
save(filename1, 'Area_groups');

filename2 ='Layer1_AP_Area';
save(filename2, 'L1_A_groups');

filename3 = 'Layer2_AP_Area';
save(filename3, 'L2_A_groups');

%% Create Output folder- Means
file = 'E:\Data\Two_Photon_Data\GCaMP6\Whisker_Stim\Analysis_Results\activity_ROIs_2013_12_16';
outputFolder = fullfile(file, 'Means');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

cd(outputFolder)

filename4= 'Area_BothLMean&SD';
save(filename4, 'Area_Mean');

filename5 ='Area_Lay1_Mean&SD';
save(filename5, 'L1_A_Mean');

filename6 = 'Area_Lay2_Mean&SD';
save(filename6, 'L2_A_Mean');

stats.hypothesis_test = h;
stats.pvalues = p;
filename7 = 'Area_ttest_results';
save(filename7, 'stats');