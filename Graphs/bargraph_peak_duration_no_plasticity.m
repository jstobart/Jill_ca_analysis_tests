%% Graphs of Peak Duration- Not for Plasticity Study
% Version 3
%  17.12.2013  JS

%open data files
DirName = 'E:\Data\Two_Photon_Data\GCaMP6\Whisker_Stim\Analysis_Results\matched_ROIs_2013_12_16';
cd(DirName)
flist=dir(fullfile(DirName, '*.mat')); % will give you list of specifed files
for ii=1:length(flist)
     load(flist(ii).name) % load files in workspace
end

%% Cell of Average Duration for each ROI type (with both layers)

 RowInfo = ['NS_AC';'NS_AP';'NS_BR';'NS_EF';'S8_AC';'S8_AP';'S8_BR';'S8_EF'];
Dur_Mean =  cellstr(RowInfo);
Dur_groups = cellstr(RowInfo);
L1_Dur_Mean =  cellstr(RowInfo);
L2_Dur_Mean =  cellstr(RowInfo);
L1_Dur_groups = cellstr(RowInfo);
L2_Dur_groups = cellstr(RowInfo);


%% determine average Duration over time (for both layers)
% AC ROIs 
%NO STIM 
%Both Layers
     for iROI = 1:length(Nostim.AC)
           y{iROI} = cell2mat(Nostim.AC{iROI}{11,2});
           y2 =cell2mat(y);  
     end
    Dur_Mean{1,2} = nanmean(y2);
    Dur_Mean{1,3} = nanstd(y2); 
    Dur_groups{1,2} = {y2};

    %Layer 1
     for iROI = 1:length(Layer1.NS_AC)
          k{iROI} = cell2mat(Layer1.NS_AC{iROI}{11,2});
         k2 =cell2mat(k);
    end
    L1_Dur_groups{1,2} = {k2};
    L1_Dur_Mean{1,2} = nanmean(k2);
    L1_Dur_Mean{1,3} = nanstd(k2);
    
    %Layer 2
    for iROI = 1:length(Layer2.NS_AC)
          q{iROI} = cell2mat(Layer2.NS_AC{iROI}{11,2});
         q2 =cell2mat(q);
    end
    L2_Dur_groups{1,2} = {q2};
    L2_Dur_Mean{1,2} = nanmean(q2);
    L2_Dur_Mean{1,3} = nanstd(q2) 
 
   clear 'y' 'y2' 'k' 'k2' 'q' 'q2'
   
   %STIM 
%Both Layers
     for iROI = 1:length(Stim8s.AC)
           y{iROI} = cell2mat(Stim8s.AC{iROI}{11,2});
           y2 =cell2mat(y);  
     end
    Dur_Mean{5,2} = nanmean(y2);
    Dur_Mean{5,3} = nanstd(y2); 
    Dur_groups{5,2} = {y2};

    %Layer 1
     for iROI = 1:length(Layer1.S_AC)
          k{iROI} = cell2mat(Layer1.S_AC{iROI}{11,2});
         k2 =cell2mat(k);
    end
    L1_Dur_groups{5,2} = {k2};
    L1_Dur_Mean{5,2} = nanmean(k2);
    L1_Dur_Mean{5,3} = nanstd(k2);
    
    %Layer 2
    for iROI = 1:length(Layer2.S_AC)
          q{iROI} = cell2mat(Layer2.S_AC{iROI}{11,2});
         q2 =cell2mat(q);
    end
    L2_Dur_groups{5,2} = {q2};
    L2_Dur_Mean{5,2} = nanmean(q2);
    L2_Dur_Mean{5,3} = nanstd(q2) 
 
   clear 'y' 'y2' 'k' 'k2' 'q' 'q2'


   %%  AP ROIs 
%NO STIM 
%Both Layers
     for iROI = 1:length(Nostim.AP)
           y{iROI} = cell2mat(Nostim.AP{iROI}{11,2});
           y2 =cell2mat(y);  
     end
    Dur_Mean{2,2} = nanmean(y2);
    Dur_Mean{2,3} = nanstd(y2); 
    Dur_groups{2,2} = {y2};
    
    %Layer 1
     for iROI = 1:length(Layer1.NS_AP)
          k{iROI} = cell2mat(Layer1.NS_AP{iROI}{11,2});
         k2 =cell2mat(k);
    end
    L1_Dur_groups{2,2} = {k2};
    L1_Dur_Mean{2,2} = nanmean(k2);
    L1_Dur_Mean{2,3} = nanstd(k2);
    
    %Layer 2
    for iROI = 1:length(Layer2.NS_AP)
          q{iROI} = cell2mat(Layer2.NS_AP{iROI}{11,2});
         q2 =cell2mat(q);
    end
    L2_Dur_groups{2,2} = {q2};
    L2_Dur_Mean{2,2} = nanmean(q2);
    L2_Dur_Mean{2,3} = nanstd(q2); 
 
   clear 'y' 'y2' 'k' 'k2' 'q' 'q2'
   
   %STIM 
%Both Layers
     for iROI = 1:length(Stim8s.AP)
           y{iROI} = cell2mat(Stim8s.AP{iROI}{11,2});
           y2 =cell2mat(y);  
     end
    Dur_Mean{6,2} = nanmean(y2);
    Dur_Mean{6,3} = nanstd(y2); 
    Dur_groups{6,2} = {y2};

    %Layer 1
     for iROI = 1:length(Layer1.S_AP)
          k{iROI} = cell2mat(Layer1.S_AP{iROI}{11,2});
         k2 =cell2mat(k);
    end
    L1_Dur_groups{6,2} = {k2};
    L1_Dur_Mean{6,2} = nanmean(k2);
    L1_Dur_Mean{6,3} = nanstd(k2);
    
    %Layer 2
    for iROI = 1:length(Layer2.S_AP)
          q{iROI} = cell2mat(Layer2.S_AP{iROI}{11,2});
         q2 =cell2mat(q);
    end
    L2_Dur_groups{6,2} = {q2};
    L2_Dur_Mean{6,2} = nanmean(q2);
    L2_Dur_Mean{6,3} = nanstd(q2); 
 
   clear 'y' 'y2' 'k' 'k2' 'q' 'q2'
   
      %%  BR ROIs 
%NO STIM 
%Both Layers
     for iROI = 1:length(Nostim.BR)
           y{iROI} = cell2mat(Nostim.BR{iROI}{11,2});
           y2 =cell2mat(y);  
     end
    Dur_Mean{3,2} = nanmean(y2);
    Dur_Mean{3,3} = nanstd(y2); 
    Dur_groups{3,2} = {y2};

    %Layer 1
     for iROI = 1:length(Layer1.NS_BR)
          k{iROI} = cell2mat(Layer1.NS_BR{iROI}{11,2});
         k2 =cell2mat(k);
    end
    L1_Dur_groups{3,2} = {k2};
    L1_Dur_Mean{3,2} = nanmean(k2);
    L1_Dur_Mean{3,3} = nanstd(k2);
    
    %Layer 2
    for iROI = 1:length(Layer2.NS_BR)
          q{iROI} = cell2mat(Layer2.NS_BR{iROI}{11,2});
         q2 =cell2mat(q);
    end
    L2_Dur_groups{3,2} = {q2};
    L2_Dur_Mean{3,2} = nanmean(q2);
    L2_Dur_Mean{3,3} = nanstd(q2); 
 
   clear 'y' 'y2' 'k' 'k2' 'q' 'q2'
   
   %STIM 
%Both Layers
     for iROI = 1:length(Stim8s.BR)
           y{iROI} = cell2mat(Stim8s.BR{iROI}{11,2});
           y2 =cell2mat(y);  
     end
    Dur_Mean{7,2} = nanmean(y2);
    Dur_Mean{7,3} = nanstd(y2); 
    Dur_groups{7,2} = {y2};

    %Layer 1
     for iROI = 1:length(Layer1.S_BR)
          k{iROI} = cell2mat(Layer1.S_BR{iROI}{11,2});
         k2 =cell2mat(k);
    end
    L1_Dur_groups{7,2} = {k2};
    L1_Dur_Mean{7,2} = nanmean(k2);
    L1_Dur_Mean{7,3} = nanstd(k2);
    
    %Layer 2
    for iROI = 1:length(Layer2.S_BR)
          q{iROI} = cell2mat(Layer2.S_BR{iROI}{11,2});
         q2 =cell2mat(q);
    end
    L2_Dur_groups{7,2} = {q2};
    L2_Dur_Mean{7,2} = nanmean(q2);
    L2_Dur_Mean{7,3} = nanstd(q2) 
 
   clear 'y' 'y2' 'k' 'k2' 'q' 'q2'

   
     %%  EF ROIs 
%NO STIM 
%Both Layers
     for iROI = 1:length(Nostim.EF)
           y{iROI} = cell2mat(Nostim.EF{iROI}{11,2});
           y2 =cell2mat(y);  
                    
         end
         
    Dur_Mean{4,2} = nanmean(y2);
    Dur_Mean{4,3} = nanstd(y2); 
    Dur_groups{4,2} = {y2};

    %Layer 1
     for iROI = 1:length(Layer1.NS_EF)
          k{iROI} = cell2mat(Layer1.NS_EF{iROI}{11,2});
         k2 =cell2mat(k);
    end
    L1_Dur_groups{4,2} = {k2};
    L1_Dur_Mean{4,2} = nanmean(k2);
    L1_Dur_Mean{4,3} = nanstd(k2);
    
    %Layer 2
    for iROI = 1:length(Layer2.NS_EF)
          q{iROI} = cell2mat(Layer2.NS_EF{iROI}{11,2});
         q2 =cell2mat(q);
    end
    L2_Dur_groups{4,2} = {q2};
    L2_Dur_Mean{4,2} = nanmean(q2);
    L2_Dur_Mean{4,3} = nanstd(q2) 
 
   clear 'y' 'y2' 'k' 'k2' 'q' 'q2'
   
   %STIM 
%Both Layers
     for iROI = 1:length(Stim8s.EF)
           y{iROI} = cell2mat(Stim8s.EF{iROI}{11,2});
           y2 =cell2mat(y);  
     end
    Dur_Mean{8,2} = nanmean(y2);
    Dur_Mean{8,3} = nanstd(y2); 
    Dur_groups{8,2} = {y2};

    %Layer 1
     for iROI = 1:length(Layer1.S_EF)
          k{iROI} = cell2mat(Layer1.S_EF{iROI}{11,2});
         k2 =cell2mat(k);
    end
    L1_Dur_groups{8,2} = {k2};
    L1_Dur_Mean{8,2} = nanmean(k2);
    L1_Dur_Mean{8,3} = nanstd(k2);
    
    %Layer 2
    for iROI = 1:length(Layer2.S_EF)
          q{iROI} = cell2mat(Layer2.S_EF{iROI}{11,2});
         q2 =cell2mat(q);
    end
    L2_Dur_groups{8,2} = {q2};
    L2_Dur_Mean{8,2} = nanmean(q2);
    L2_Dur_Mean{8,3} = nanstd(q2) 
 
   clear 'y' 'y2' 'k' 'k2' 'q' 'q2'
   
%%  Bar Graphs
%BOTH LAYERS
No_stim = cell2mat(Dur_Mean(1:4,2)');
Stim =cell2mat(Dur_Mean(5:8,2)');
No_stim_SD = cell2mat(Dur_Mean(1:4,3)');
Stim_SD =cell2mat(Dur_Mean(5:8,3)');

Both = [No_stim' Stim'];
Both_SD = [No_stim_SD' Stim_SD'];

%Layer 1
L1_No_stim = cell2mat(L1_Dur_Mean(1:4,2)');
L1_Stim =cell2mat(L1_Dur_Mean(5:8,2)');
L1_No_stim_SD = cell2mat(L1_Dur_Mean(1:4,3)');
L1_Stim_SD =cell2mat(L1_Dur_Mean(5:8,3)');

L1 = [L1_No_stim' L1_Stim'];
L1_SD = [L1_No_stim_SD' L1_Stim_SD'];

%Layer 2
L2_No_stim = cell2mat(L2_Dur_Mean(1:4,2)');
L2_Stim =cell2mat(L2_Dur_Mean(5:8,2)');
L2_No_stim_SD = cell2mat(L2_Dur_Mean(1:4,3)');
L2_Stim_SD =cell2mat(L2_Dur_Mean(5:8,3)');

L2 = [L2_No_stim' L2_Stim'];
L2_SD = [L2_No_stim_SD' L2_Stim_SD'];

%Figures
figure ('name','Peak Duration','numbertitle','off')
subplot (3,1,1)
hold on
barwitherr(Both, Both_SD)
set(gca,'XTickLabel', '')
set(gca,'YLim', [0 40],'YTick',0:5:40)
ylabel('Duration (s)')
xlabel('ROI Type')
title('Both Layers')
legend('NoStim','Stim') 

subplot (3,1,2)
hold on
barwitherr(L1,L1_SD)
set(gca,'XTickLabel', '')
set(gca,'YLim', [0 40],'YTick',0:5:40)
ylabel('Duration (s)')
xlabel('ROI Type')
title('Layer1')
legend('NoStim','Stim') 

subplot (3,1,3)
hold on
barwitherr(L2,L2_SD)
set(gca,'XTickLabel', '')
set(gca,'YLim', [0 40],'YTick',0:5:40)
ylabel('Duration (s)')
xlabel('ROI Type')
title('Layer2')
legend('NoStim','Stim')




%% Statistics- T-test for stim vs. nostim
%Both Layers

%AC ROIs
AC_NS = cell2mat(Dur_groups{1,2});
AC_S = cell2mat(Dur_groups{5,2});
AC = padcat(AC_NS, AC_S);
AC1= AC(1,:);
AC2= AC(2,:);

[h(1,1), p(1,1)] =ttest(AC1, AC2)

%AP ROIs
AP_NS = cell2mat(Dur_groups{2,2});
AP_S = cell2mat(Dur_groups{6,2});
AP = padcat(AP_NS, AP_S);
AP1= AP(1,:);
AP2= AP(2,:);

[h(2,1), p(2,1)] =ttest(AP1, AP2)

%BR ROIs
BR_NS = cell2mat(Dur_groups{3,2});
BR_S = cell2mat(Dur_groups{7,2});
BR = padcat(BR_NS, BR_S);
BR1= BR(1,:);
BR2= BR(2,:);

[h(3,1), p(3,1)] =ttest(BR1, BR2)

%EF ROIs
EF_NS = cell2mat(Dur_groups{4,2});
EF_S = cell2mat(Dur_groups{8,2});
EF = padcat(EF_NS, EF_S);
EF1= EF(1,:);
EF2= EF(2,:);

[h(4,1), p(4,1)] =ttest(EF1, EF2)

%% Layer 1

%AC ROIs
L1_AC_NS = cell2mat(L1_Dur_groups{1,2});
L1_AC_S = cell2mat(L1_Dur_groups{5,2});
L1_AC = padcat(L1_AC_NS, L1_AC_S);
L1_AC1= L1_AC(1,:);
L1_AC2= L1_AC(2,:);

[h(1,2), p(1,2)] =ttest(L1_AC1, L1_AC2)

%AP ROIs
L1_AP_NS = cell2mat(L1_Dur_groups{2,2});
L1_AP_S = cell2mat(L1_Dur_groups{6,2});
L1_AP = padcat(L1_AP_NS, L1_AP_S);
L1_AP1= L1_AP(1,:);
L1_AP2= L1_AP(2,:);

[h(2,2), p(2,2)] =ttest(L1_AP1, L1_AP2)

%BR ROIs
L1_BR_NS = cell2mat(L1_Dur_groups{3,2});
L1_BR_S = cell2mat(L1_Dur_groups{7,2});
L1_BR = padcat(L1_BR_NS, L1_BR_S);
L1_BR1= L1_BR(1,:);
L1_BR2= L1_BR(2,:);

[h(3,2), p(3,2)] =ttest(L1_BR1, L1_BR2)

%EF ROIs
L1_EF_NS = cell2mat(L1_Dur_groups{4,2});
L1_EF_S = cell2mat(L1_Dur_groups{8,2});
L1_EF = padcat(L1_EF_NS, L1_EF_S);
L1_EF1= L1_EF(1,:);
L1_EF2= L1_EF(2,:);

[h(4,2), p(4,2)] =ttest(L1_EF1, L1_EF2)

%% Layer 2

%AC ROIs
L2_AC_NS = cell2mat(L2_Dur_groups{1,2});
L2_AC_S = cell2mat(L2_Dur_groups{5,2});
L2_AC = padcat(L2_AC_NS, L2_AC_S);
L2_AC1= L2_AC(1,:);
L2_AC2= L2_AC(2,:);

[h(1,3), p(1,3)] =ttest(L2_AC1, L2_AC2)

%AP ROIs
L2_AP_NS = cell2mat(L2_Dur_groups{2,2});
L2_AP_S = cell2mat(L2_Dur_groups{6,2});
L2_AP = padcat(L2_AP_NS, L2_AP_S);
L2_AP1= L2_AP(1,:);
L2_AP2= L2_AP(2,:);

[h(2,3), p(2,3)] =ttest(L2_AP1, L2_AP2)

%BR ROIs
L2_BR_NS = cell2mat(L2_Dur_groups{3,2});
L2_BR_S = cell2mat(L2_Dur_groups{7,2});
L2_BR = padcat(L2_BR_NS, L2_BR_S);
L2_BR1= L2_BR(1,:);
L2_BR2= L2_BR(2,:);

[h(3,3), p(3,3)] =ttest(L2_BR1, L2_BR2)

%EF ROIs
L2_EF_NS = cell2mat(L2_Dur_groups{4,2});
L2_EF_S = cell2mat(L2_Dur_groups{8,2});
L2_EF = padcat(L2_EF_NS, L2_EF_S);
L2_EF1= L2_EF(1,:);
L2_EF2= L2_EF(2,:);

[h(4,3), p(4,3)] =ttest(L2_EF1, L2_EF2)

%% histogram graph
% AC
s = 0:2:124;

figure ('name', 'AC_number of signals vs. duration', 'numbertitle', 'off')
subplot (3,1,1)
hist(AC_NS,s)
r = findobj(gca,'Type','patch');
set(r,'FaceColor','r','EdgeColor','w','facealpha',0.75)
hold on;
hist(AC_S,s)

r1 = findobj(gca,'Type','patch');
set(r1,'facealpha',0.75);

set(gca,'Xlim',[0 125])
ylabel('Number of Signals')
xlabel('Duration (s)')
legend('NoStim','8secStim')
title('Both Layers')

subplot (3,1,2)
hist(L1_AC_NS,s)
r = findobj(gca,'Type','patch');
set(r,'FaceColor','r','EdgeColor','w','facealpha',0.75)
hold on;
hist(L1_AC_S,s)

r1 = findobj(gca,'Type','patch');
set(r1,'facealpha',0.75);

set(gca,'Xlim',[0 125])
ylabel('Number of Signals')
xlabel('Duration (s)')
legend('NoStim','8secStim')
title('Layer1')

subplot (3,1,3)
hist(L2_AC_NS,s)
r = findobj(gca,'Type','patch');
set(r,'FaceColor','r','EdgeColor','w','facealpha',0.75)
hold on;
hist(L2_AC_S,s)

r1 = findobj(gca,'Type','patch');
set(r1,'facealpha',0.75);

set(gca,'Xlim',[0 125])
ylabel('Number of Signals')
xlabel('Duration (s)')
legend('NoStim','8secStim')
title('Layer2')

%% AP
s = 0:2:124;

figure ('name', 'AP_number of signals vs. duration', 'numbertitle', 'off')
subplot (3,1,1)
hist(AP_NS,s)
r = findobj(gca,'Type','patch');
set(r,'FaceColor','r','EdgeColor','w','facealpha',0.75)
hold on;
hist(AP_S,s)

r1 = findobj(gca,'Type','patch');
set(r1,'facealpha',0.75);

set(gca,'Xlim',[0 125])
ylabel('Number of Signals')
xlabel('Duration (s)')
legend('NoStim','8secStim')
title('Both Layers')

subplot (3,1,2)
hist(L1_AP_NS,s)
r = findobj(gca,'Type','patch');
set(r,'FaceColor','r','EdgeColor','w','facealpha',0.75)
hold on;
hist(L1_AP_S,s)

r1 = findobj(gca,'Type','patch');
set(r1,'facealpha',0.75);

set(gca,'Xlim',[0 125])
ylabel('Number of Signals')
xlabel('Duration (s)')
legend('NoStim','8secStim')
title('Layer1')

subplot (3,1,3)
hist(L2_AP_NS,s)
r = findobj(gca,'Type','patch');
set(r,'FaceColor','r','EdgeColor','w','facealpha',0.75)
hold on;
hist(L2_AP_S,s)

r1 = findobj(gca,'Type','patch');
set(r1,'facealpha',0.75);

set(gca,'Xlim',[0 125])
ylabel('Number of Signals')
xlabel('Duration (s)')
legend('NoStim','8secStim')
title('Layer2')

%% BR

s = 0:2:124;

figure ('name', 'BR_number of signals vs. duration', 'numbertitle', 'off')
subplot (3,1,1)
hist(BR_NS,s)
r = findobj(gca,'Type','patch');
set(r,'FaceColor','r','EdgeColor','w','facealpha',0.75)
hold on;
hist(BR_S,s)

r1 = findobj(gca,'Type','patch');
set(r1,'facealpha',0.75);

set(gca,'Xlim',[0 125])
ylabel('Number of Signals')
xlabel('Duration (s)')
legend('NoStim','8secStim')
title('Both Layers')

subplot (3,1,2)
hist(L1_BR_NS,s)
r = findobj(gca,'Type','patch');
set(r,'FaceColor','r','EdgeColor','w','facealpha',0.75)
hold on;
hist(L1_BR_S,s)

r1 = findobj(gca,'Type','patch');
set(r1,'facealpha',0.75);

set(gca,'Xlim',[0 125])
ylabel('Number of Signals')
xlabel('Duration (s)')
legend('NoStim','8secStim')
title('Layer1')

subplot (3,1,3)
hist(L2_BR_NS,s)
r = findobj(gca,'Type','patch');
set(r,'FaceColor','r','EdgeColor','w','facealpha',0.75)
hold on;
hist(L2_BR_S,s)

r1 = findobj(gca,'Type','patch');
set(r1,'facealpha',0.75);

set(gca,'Xlim',[0 125])
ylabel('Number of Signals')
xlabel('Duration (s)')
legend('NoStim','8secStim')
title('Layer2')

%% EF

s = 0:2:124;

figure ('name', 'EF_number of signals vs. duration', 'numbertitle', 'off')
subplot (3,1,1)
hist(EF_NS,s)
r = findobj(gca,'Type','patch');
set(r,'FaceColor','r','EdgeColor','w','facealpha',0.75)
hold on;
hist(EF_S,s)

r1 = findobj(gca,'Type','patch');
set(r1,'facealpha',0.75);

set(gca,'Xlim',[0 125])
ylabel('Number of Signals')
xlabel('Duration (s)')
legend('NoStim','8secStim')
title('Both Layers')

subplot (3,1,2)
hist(L1_EF_NS,s)
r = findobj(gca,'Type','patch');
set(r,'FaceColor','r','EdgeColor','w','facealpha',0.75)
hold on;
hist(L1_EF_S,s)

r1 = findobj(gca,'Type','patch');
set(r1,'facealpha',0.75);

set(gca,'Xlim',[0 125])
ylabel('Number of Signals')
xlabel('Duration (s)')
legend('NoStim','8secStim')
title('Layer1')

subplot (3,1,3)
hist(L2_EF_NS,s)
r = findobj(gca,'Type','patch');
set(r,'FaceColor','r','EdgeColor','w','facealpha',0.75)
hold on;
hist(L2_EF_S,s)

r1 = findobj(gca,'Type','patch');
set(r1,'facealpha',0.75);

set(gca,'Xlim',[0 125])
ylabel('Number of Signals')
xlabel('Duration (s)')
legend('NoStim','8secStim')
title('Layer2')

%% Create Output folder- groups
file = 'E:\Data\Two_Photon_Data\GCaMP6\Whisker_Stim\Analysis_Results\matched_ROIs_2013_12_16';
outputFolder = fullfile(file, 'Data_Groups');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

cd(outputFolder)

filename1= 'Both_Layers_Dur';
save(filename1, 'Dur_groups');

filename2 ='Layer1_Dur';
save(filename2, 'L1_Dur_groups');

filename3 = 'Layer2_Dur';
save(filename3, 'L2_Dur_groups');

%% Create Output folder- Means
file = 'E:\Data\Two_Photon_Data\GCaMP6\Whisker_Stim\Analysis_Results\matched_ROIs_2013_12_16';
outputFolder = fullfile(file, 'Means');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

cd(outputFolder)

filename4= 'Dur_BothL_Mean&SD';
save(filename4, 'Dur_Mean');

filename5 ='Dur_Lay1_Mean&SD';
save(filename5, 'L1_Dur_Mean');

filename6 = 'Dur_Lay2_Mean&SD';
save(filename6, 'L2_Dur_Mean');

stats.hypothesis_test = h;
stats.pvalues = p;
filename7 = 'Dur_ttest_results';
save(filename7, 'stats');

