%% Graphs of AUC for each spike vs. duration
% Version 2
%  19.12.2013  JS

%open data files
DirName = 'E:\Data\Two_Photon_Data\GCaMP6\Whisker_Stim\Analysis_Results\matched_ROIs_2013_12_16\Data_Groups';
cd(DirName)
flist=dir(fullfile(DirName, '*.mat')); % will give you list of specifed files
for ii=1:length(flist)
     load(flist(ii).name) % load files in workspace
end

%% determine peak values and the duration
% BR
S1= repmat(20,1,length(PA_groups{3,2}));  % circle size on the plot
S2= repmat(20,1,length(PA_groups{7,2}));  % circle size on the plot

%Both Layers- Peak AUC vs Duration
Both_BR_AUC_NS = cell2mat(PA_groups{3,2});
Both_BR_Dur_NS = cell2mat(Dur_groups{3,2});
Both_BR_AUC_S = cell2mat(PA_groups{7,2});
Both_BR_Dur_S = cell2mat(Dur_groups{7,2});


figure ('name','BR ROI-Peak AUC vs. Duration','numbertitle','off')
subplot (3,2,1)
hold on
scatter(Both_BR_Dur_NS,Both_BR_AUC_NS,S1,'b')
scatter(Both_BR_Dur_S,Both_BR_AUC_S,S2,'r')
ylabel('AUC')
xlabel('Duration (s)')
title('Both Layers')
set(gca,'XLim',[0 180])
set(gca,'YLim', [0 10000],'YTick',0:1000:10000)
legend('Nostim','Stim')  

% Peak number
NS = length(Both_BR_Dur_NS);
S= length(Both_BR_Dur_S); 

subplot (3,2,2)
hold on
bar(1,NS,'b')
hold on
bar(2,S,'r')
set(gca,'XTickLabel', '')
set(gca,'YLim', [0 2000],'YTick',0:500:2000)
ylabel('# of Peaks')
xlabel('ROI Type')
title('Both Layers')
legend('NoStim','Stim')

%Layer 1
L1_BR_AUC_NS = cell2mat(L1_PA_groups{3,2});
L1_BR_Dur_NS = cell2mat(L1_Dur_groups{3,2});
L1_BR_AUC_S = cell2mat(L1_PA_groups{7,2});
L1_BR_Dur_S = cell2mat(L1_Dur_groups{7,2});


subplot (3,2,3)
hold on
scatter(L1_BR_Dur_NS,L1_BR_AUC_NS,S1,'b')
scatter(L1_BR_Dur_S,L1_BR_AUC_S,S2,'r')
ylabel('AUC')
xlabel('Duration (s)')
title('Layer1')
set(gca,'XLim',[0 180])
set(gca,'YLim', [0 10000],'YTick',0:1000:10000)
legend('Nostim','Stim')  

% Peak number
L1_NS = length(L1_BR_Dur_NS);
L1_S= length(L1_BR_Dur_S); 

subplot (3,2,4)
hold on
bar(1,L1_NS,'b')
hold on
bar(2,L1_S,'r')
set(gca,'XTickLabel', '')
set(gca,'YLim', [0 2000],'YTick',0:500:2000)
ylabel('# of Peaks')
xlabel('ROI Type')
title('Layer1')
legend('NoStim','Stim')

%Layer 2
 L2_BR_AUC_NS = cell2mat(L2_PA_groups{3,2});
L2_BR_Dur_NS = cell2mat(L2_Dur_groups{3,2});
L2_BR_AUC_S = cell2mat(L2_PA_groups{7,2});
L2_BR_Dur_S = cell2mat(L2_Dur_groups{7,2});


subplot (3,2,5)
hold on
scatter(L2_BR_Dur_NS,L2_BR_AUC_NS,S1,'b')
scatter(L2_BR_Dur_S,L2_BR_AUC_S,S2,'r')
ylabel('AUC')
xlabel('Duration (s)')
title('Layer2')
set(gca,'XLim',[0 180])
set(gca,'YLim', [0 10000],'YTick',0:1000:10000)
legend('Nostim','Stim')  

% Peak number
L2_NS = length(L2_BR_Dur_NS);
L2_S= length(L2_BR_Dur_S); 

subplot (3,2,6)
hold on
bar(1,L2_NS,'b')
hold on
bar(2,L2_S,'r')
set(gca,'XTickLabel', '')
set(gca,'YLim', [0 2000],'YTick',0:500:2000)
ylabel('# of Peaks')
xlabel('ROI Type')
title('Layer2')
legend('NoStim','Stim')









  