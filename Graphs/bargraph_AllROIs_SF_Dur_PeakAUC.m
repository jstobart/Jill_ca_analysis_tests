%% Graphs of All ROIs together-
%Signal Frequency,Peak AUC, Duration 
% Version 4
%  20.01.2014  JS

%open data files
DirName = 'E:\Data\Two_Photon_Data\GCaMP6\Whisker_Stim\Analysis_Results\activity_ROIs_2013_12_16\Data_Groups';
cd(DirName)
flist=dir(fullfile(DirName, '*.mat')); % will give you list of specifed files
for ii=1:length(flist)
     load(flist(ii).name) % load files in workspace
end
DirName = 'E:\Data\Two_Photon_Data\GCaMP6\Whisker_Stim\Analysis_Results\activity_ROIs_2013_12_16';
cd(DirName)
flist=dir(fullfile(DirName, '*.mat')); % will give you list of specifed files
for ii=1:length(flist)
     load(flist(ii).name) % load files in workspace
end


%% determine average signal frequency over time (ALL ROIs from each layer)

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

%% number of ROIs/ cell  Needs revising to be accurate for each CELL!
%num_L1_NS = length(PA_L1_NS)/length(Layer1.NS_AC);
%num_L1_S = length(PA_L1_S)/length(Layer1.S_AC);
%num_L2_NS = length(PA_L2_NS)/length(Layer2.NS_AC);
%num_L2_S = length(PA_L2_S)/length(Layer2.S_AC);

%num_NS = [num_L1_NS, num_L2_NS];
%num_S = [num_L1_S, num_L2_S];
%num = [num_NS' num_S'];
   

%% Statistics- T-test for stim vs. nostim
SF_L1_NS = [];
SF_L1_S = [];
SF_L2_NS = [];
SF_L2_S = [];
%Layer 1, SF
for ii = 1:4
SF_L1_NS = horzcat(SF_L1_NS,cell2mat(L1_SF_groups{ii,2}));
end
for xx = 5:8
    SF_L1_S = horzcat(SF_L1_S,cell2mat(L1_SF_groups{xx,2}));
end

A = padcat(SF_L1_NS, SF_L1_S);
A1= A(1,:);
A2= A(2,:);

[h(1,1), p(1,1)] =ttest(A1, A2);

%Layer 2, SF
for ii = 1:4
SF_L2_NS = horzcat(SF_L2_NS,cell2mat(L2_SF_groups{ii,2}));
end
for xx = 5:8
    SF_L2_S = horzcat(SF_L2_S,cell2mat(L2_SF_groups{xx,2}));
end

A3 = padcat(SF_L2_NS, SF_L2_S);
A4= A3(1,:);
A5= A3(2,:);

[h(1,2), p(1,2)] =ttest(A4, A5);

%Duration
Dur_L1_NS = [];
Dur_L1_S = [];
Dur_L2_NS = [];
Dur_L2_S = [];
%Layer 1, Dur
for ii = 1:4
Dur_L1_NS = horzcat(Dur_L1_NS,cell2mat(L1_Dur_groups{ii,2}));
end
for xx = 5:8
    Dur_L1_S = horzcat(Dur_L1_S,cell2mat(L1_Dur_groups{xx,2}));
end

B = padcat(Dur_L1_NS, Dur_L1_S);
B1= B(1,:);
B2= B(2,:);

[h(2,1), p(2,1)] =ttest(B1, B2);

%Layer 2, Dur
for ii = 1:4
Dur_L2_NS = horzcat(Dur_L2_NS,cell2mat(L2_Dur_groups{ii,2}));
end
for xx = 5:8
    Dur_L2_S = horzcat(Dur_L2_S,cell2mat(L2_Dur_groups{xx,2}));
end

B3 = padcat(Dur_L2_NS, Dur_L2_S);
B4= B3(1,:);
B5= B3(2,:);

[h(2,2), p(2,2)] =ttest(B4, B5);

%Peak Area
PA_L1_NS = [];
PA_L1_S = [];
PA_L2_NS = [];
PA_L2_S = [];
%Layer 1, PA
for ii = 1:4
PA_L1_NS = horzcat(PA_L1_NS,cell2mat(L1_PA_groups{ii,2}));
end
for xx = 5:8
    PA_L1_S = horzcat(PA_L1_S,cell2mat(L1_PA_groups{xx,2}));
end

C = padcat(PA_L1_NS, PA_L1_S);
C1= C(1,:);
C2= C(2,:);

[h(3,1), p(3,1)] =ttest(C1, C2);

%Layer 2, PA
for ii = 1:4
PA_L2_NS = horzcat(PA_L2_NS,cell2mat(L2_PA_groups{ii,2}));
end
for xx = 5:8
    PA_L2_S = horzcat(PA_L2_S,cell2mat(L2_PA_groups{xx,2}));
end

C3 = padcat(PA_L2_NS, PA_L2_S);
C4= C3(1,:);
C5= C3(2,:);

[h(3,2), p(3,2)] =ttest(C4, C5);



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



