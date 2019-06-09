%% Graphs of AP ROIs-
%Signal Frequency (new calculation)- one number per trial
% Peak Amp, Duration 
% Version 5
%  30.10.2014  JS

open data files
DirName = 'E:\Data\Two_Photon_Data\GFAP_GCaMP6\Whisker_Stim\LongStim\Analysis_Results_Old\activity_ROIs_2013_12_16';
cd(DirName)
flist=dir(fullfile(DirName, '*.mat')); % will give you list of specifed files
for ii=1:length(flist)
     load(flist(ii).name) % load files in workspace
end
% DirName = 'E:\Data\Two_Photon_Data\GCaMP6\Whisker_Stim\Analysis_Results\activity_ROIs_2013_12_16';
% cd(DirName)
% flist=dir(fullfile(DirName, '*.mat')); % will give you list of specifed files
% for ii=1:length(flist)
%      load(flist(ii).name) % load files in workspace
% end


%% determine average Duration over time 
for iROI = 1:length(Layer1.NS_AP)
    for itrial = 1:length(Layer1.NS_AP{iROI}{11,2})
        L1_NSDur1{iROI}= cell2mat(Layer1.NS_AP{iROI}{11,2});
        L1_NSDur = cell2mat(L1_NSDur1);
        
        L1_NS_Freq1{iROI}{itrial} = length(Layer1.NS_AP{iROI}{11,2}{itrial})/3; % calculation of the frequency per trial
        L1_NS_Freq2{iROI} = cell2mat(L1_NS_Freq1{iROI});    
        L1_NSFreq = cell2mat(L1_NS_Freq2);
        
        L1_NSDur1{iROI}= cell2mat(Layer1.NS_AP{iROI}{11,2});
        L1_NSDur = cell2mat(L1_NSDur1);
        
        L1_NSDurA1{iROI}{itrial}= nanmean(Layer1.NS_AP{iROI}{11,2}{itrial});
        L1_NSDurA2{iROI}= nanmean(cell2mat(L1_NSDurA1{iROI}));
    end
end
L1_MeanNSDur = nanmean(L1_NSDur);
L1_SDNSDur = nanstd(L1_NSDur);
L1_MeanNSFreq = nanmean(L1_NSFreq);
L1_SDNSFreq = nanstd(L1_NSFreq);
        
for iROI = 1:length(Layer1.S_AP)
    for itrial = 1:length(Layer1.S_AP{iROI}{11,2})
        L1_SDur1{iROI}= cell2mat(Layer1.S_AP{iROI}{11,2});
        L1_SDur = cell2mat(L1_SDur1);
                
        L1_SDurA1{iROI}{itrial}= nanmean(Layer1.S_AP{iROI}{11,2}{itrial});
        L1_SDurA2{iROI}= nanmean(cell2mat(L1_SDurA1{iROI}));
    end
end
L1_MeanSDur = nanmean(L1_SDur);
L1_SDSDur = nanstd(L1_SDur);
        

for iROI = 1:length(Layer2.NS_AP)
    for itrial = 1:length(Layer2.NS_AP{iROI}{11,2})
        L2_NSDur1{iROI}= cell2mat(Layer2.NS_AP{iROI}{11,2});
        L2_NSDur = cell2mat(L2_NSDur1);
                
        L2_NSDurA1{iROI}{itrial}= nanmean(Layer2.NS_AP{iROI}{11,2}{itrial});
        L2_NSDurA2{iROI}= nanmean(cell2mat(L2_NSDurA1{iROI}));
    end
end
L2_MeanNSDur = nanmean(L2_NSDur);
L2_SDNSDur = nanstd(L2_NSDur);
        
for iROI = 1:length(Layer2.S_AP)
    for itrial = 1:length(Layer2.S_AP{iROI}{11,2})
        L2_SDur1{iROI}= cell2mat(Layer2.S_AP{iROI}{11,2});
        L2_SDur = cell2mat(L2_SDur1);
                
        L2_SDurA1{iROI}{itrial}= nanmean(Layer2.S_AP{iROI}{11,2}{itrial});
        L2_SDurA2{iROI}= nanmean(cell2mat(L2_SDurA1{iROI}));
    end
end
L2_MeanSDur = nanmean(L2_SDur);
L2_SDSDur = nanstd(L2_SDur);

%%
for ii = 1:4
L1_NS_SF = mean(cell2mat(L1_SF_groups{ii,2}));
L1_NS_SF_SD = std(cell2mat(L1_SF_groups{ii,2}));
L2_NS_SF = mean(cell2mat(L2_SF_groups{ii,2}));
L2_NS_SF_SD = std(cell2mat(L2_SF_groups{ii,2}));
end
for xx = 5:8
    L1_S_SF = mean(cell2mat(L1_SF_groups{xx,2}));
    L1_S_SF_SD = std(cell2mat(L1_SF_groups{xx,2}));
    L2_S_SF = mean(cell2mat(L2_SF_groups{xx,2}));
    L2_S_SF_SD = std(cell2mat(L2_SF_groups{xx,2}));
end
SF_NS = [L1_NS_SF, L2_NS_SF];
SF_NS_SD = [L1_NS_SF_SD, L2_NS_SF_SD];

SF_S = [L1_S_SF, L2_S_SF];
SF_S_SD = [L1_S_SF_SD, L2_S_SF_SD];

SF = [SF_NS' SF_S'];
SF_SD = [SF_NS_SD' SF_S_SD'];

L1_SF = [L1_NS_SF, L1_S_SF];
L2_SF = [L2_NS_SF, L2_S_SF];
L1_SF_SD = [L1_NS_SF_SD, L1_S_SF_SD];
L2_SF_SD = [L2_NS_SF_SD, L2_S_SF_SD];
SF_L = [L1_SF' L2_SF'];
SF_L_SD = [L1_SF_SD' L2_SF_SD'];

%% determine average duration over time (ALL ROIs from each layer)

for ii = 1:4
L1_NS_Dur = mean(cell2mat(L1_Dur_groups{ii,2}));
L1_NS_Dur_SD = std(cell2mat(L1_Dur_groups{ii,2}));
L2_NS_Dur = mean(cell2mat(L2_Dur_groups{ii,2}));
L2_NS_Dur_SD = std(cell2mat(L2_Dur_groups{ii,2}));
end
for xx = 5:8
    L1_S_Dur = mean(cell2mat(L1_Dur_groups{xx,2}));
    L1_S_Dur_SD = std(cell2mat(L1_Dur_groups{xx,2}));
    L2_S_Dur = mean(cell2mat(L2_Dur_groups{xx,2}));
    L2_S_Dur_SD = std(cell2mat(L2_Dur_groups{xx,2}));
end
Dur_NS = [L1_NS_Dur, L2_NS_Dur];
Dur_NS_SD = [L1_NS_Dur_SD, L2_NS_Dur_SD];

Dur_S = [L1_S_Dur, L2_S_Dur];
Dur_S_SD = [L1_S_Dur_SD, L2_S_Dur_SD];

Dur = [Dur_NS' Dur_S'];
Dur_SD = [Dur_NS_SD' Dur_S_SD'];

L1_Dur = [L1_NS_Dur, L1_S_Dur];
L2_Dur = [L2_NS_Dur, L2_S_Dur];
L1_Dur_SD = [L1_NS_Dur_SD, L1_S_Dur_SD];
L2_Dur_SD = [L2_NS_Dur_SD, L2_S_Dur_SD];
Dur_L = [L1_Dur' L2_Dur'];
Dur_L_SD = [L1_Dur_SD' L2_Dur_SD'];

%% determine average Peak AUCover time (ALL ROIs from each layer)

for ii = 1:4
L1_NS_PA = mean(cell2mat(L1_PA_groups{ii,2}));
L1_NS_PA_SD = std(cell2mat(L1_PA_groups{ii,2}));
L2_NS_PA = mean(cell2mat(L2_PA_groups{ii,2}));
L2_NS_PA_SD = std(cell2mat(L2_PA_groups{ii,2}));
end
for xx = 5:8
    L1_S_PA = mean(cell2mat(L1_PA_groups{xx,2}));
    L1_S_PA_SD = std(cell2mat(L1_PA_groups{xx,2}));
    L2_S_PA = mean(cell2mat(L2_PA_groups{xx,2}));
    L2_S_PA_SD = std(cell2mat(L2_PA_groups{xx,2}));
end
PA_NS = [L1_NS_PA, L2_NS_PA];
PA_NS_SD = [L1_NS_PA_SD, L2_NS_PA_SD];

PA_S = [L1_S_PA, L2_S_PA];
PA_S_SD = [L1_S_PA_SD, L2_S_PA_SD];

L1_PA = [L1_NS_PA, L1_S_PA];
L2_PA = [L2_NS_PA, L2_S_PA];
L1_PA_SD = [L1_NS_PA_SD, L1_S_PA_SD];
L2_PA_SD = [L2_NS_PA_SD, L2_S_PA_SD];
PA_L = [L1_PA' L2_PA'];
PA_L_SD = [L1_PA_SD' L2_PA_SD'];

PA = [PA_NS' PA_S'];
PA_SD = [PA_NS_SD' PA_S_SD'];



%%  Bar Graphs
%Figures
figure ('name','Bar Graphs-All ROIs','numbertitle','off')

subplot (3,1,1)
hold on
barwitherr(PA,PA_SD)
set(gca,'XTickLabel', '')
set(gca,'YLim', [0 1200],'YTick',0:200:1200)
xlabel('Layer 1                                           Layer 2')
ylabel('Signal AUC')
%title('Signal AUC')
legend('No-Stim','Whisker-Stim')

subplot (3,1,2)
hold on
barwitherr(Dur,Dur_SD)
set(gca,'XTickLabel', '')
set(gca,'YLim', [0 30],'YTick',0:5:30)
xlabel('Layer 1                                                                       Layer 2')
ylabel('Duration (s)')
%title('Signal Duration')
legend('No-Stim','Whisker-Stim')

subplot (3,1,3)
hold on
barwitherr(SF,SF_SD)
set(gca,'XTickLabel', '')
set(gca,'YLim', [0 2.5],'YTick',0:0.5:2.5)
xlabel('Layer 1                                                                       Layer 2')
ylabel('Frequency (signals/min)')
%title('Signal Frequency')
legend('No-Stim','Whisker-Stim')

%%  Bar Graphs
%Figures
figure ('name','Bar Graphs-NS vs Stim','numbertitle','off')

subplot (3,1,1)
hold on
barwitherr(PA_L,PA_L_SD)
set(gca,'XTickLabel', '')
set(gca,'YLim', [0 1200],'YTick',0:200:1200)
xlabel('No Stim                                           Whisker Stim')
ylabel('Signal AUC')
%title('Signal AUC')
legend('Layer 1','Layer 2/3')

subplot (3,1,2)
hold on
barwitherr(Dur_L,Dur_L_SD)
set(gca,'XTickLabel', '')
set(gca,'YLim', [0 30],'YTick',0:5:30)
xlabel('No Stim                                           Whisker Stim')
ylabel('Duration (s)')
%title('Signal Duration')
legend('Layer 1','Layer 2/3')

subplot (3,1,3)
hold on
barwitherr(SF_L,SF_L_SD)
set(gca,'XTickLabel', '')
set(gca,'YLim', [0 2.5],'YTick',0:0.5:2.5)
xlabel('No Stim                                           Whisker Stim')
ylabel('Frequency (signals/min)')
%title('Signal Frequency')
legend('Layer 1','Layer 2/3')


%%
Dur_mat = padcat(B1,B2,B4,B5);

%Dur_matrix(:,1) = B1';
%Dur_matrix(:,2) = B2';
%Dur_matrix(:,3) = B4';
%Dur_matrix(:,4) = B5';
Dur_matrix = Dur_mat';
[p_v,a,s] = anova1(Dur_matrix);
[c,m,f,nms] = multcompare(s);

SF_mat = padcat(A1,A2,A4,A5);
SF_matrix = SF_mat';
[p_v2,a2,s2] = anova1(SF_matrix);
[c2,m2,f2,nms2] = multcompare(s2);

AUC_mat = padcat(C1,C2,C4,C5);
AUC_matrix = AUC_mat';
[p_v3,a3,s3] = anova1(AUC_matrix);
[c3,m3,f3,nms3] = multcompare(s3);



