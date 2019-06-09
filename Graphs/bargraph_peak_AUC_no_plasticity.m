%% Graphs of Peak AUC- Not for Plasticity Study
% Version 3
%  17.12.2013  JS

%open data files
DirName = 'E:\Data\Two_Photon_Data\GCaMP6\Whisker_Stim\Analysis_Results\activity_ROIs_2013_12_16';
cd(DirName)
flist=dir(fullfile(DirName, '*.mat')); % will give you list of specifed files
for ii=1:length(flist)
     load(flist(ii).name) % load files in workspace
end

%% Cell of Average signal Frequencies for each ROI type (with both layers)

 RowInfo = ['NS_AC';'NS_AP';'NS_BR';'NS_EF';'S8_AC';'S8_AP';'S8_BR';'S8_EF'];
PA_Mean =  cellstr(RowInfo);
PA_groups = cellstr(RowInfo);
L1_PA_Mean =  cellstr(RowInfo);
L2_PA_Mean =  cellstr(RowInfo);
L1_PA_groups = cellstr(RowInfo);
L2_PA_groups = cellstr(RowInfo);


%% determine average peak frequency over time (for both layers)
% AC ROIs 
%NO STIM 
%Both Layers
     for iROI = 1:length(Nostim.AC)
         for itrial = 1:length(Nostim.AC{iROI}{14,2})
         for ipeak = 1:length(Nostim.AC{iROI}{14,2}{itrial})
           y{iROI}{itrial} = cell2mat(Nostim.AC{iROI}{14,2}{itrial});
           y2{iROI} =cell2mat(y{iROI});
           y3 = cell2mat(y2);
         end
         end
     end
    
     
    PA_Mean{1,2} = nanmean(y3);
    PA_Mean{1,3} = nanstd(y3); 
    PA_groups{1,2} = {y3};

    %Layer 1
     for iROI = 1:length(Layer1.NS_AC)
         for itrial = 1:length(Layer1.NS_AC{iROI}{14,2})
         for ipeak = 1:length(Layer1.NS_AC{iROI}{14,2}{itrial})
          k{iROI}{itrial} =cell2mat(Layer1.NS_AC{iROI}{14,2}{itrial});
          k2{iROI} =cell2mat(k{iROI});
          k3 = cell2mat(k2);
         end
         end
    end
    L1_PA_groups{1,2} = {k3};
    L1_PA_Mean{1,2} = nanmean(k3);
    L1_PA_Mean{1,3} = nanstd(k3);
    
    %Layer 2
    for iROI = 1:length(Layer2.NS_AC)
        for itrial = 1:length(Layer2.NS_AC{iROI}{14,2})
         for ipeak = 1:length(Layer2.NS_AC{iROI}{14,2}{itrial})
          q{iROI}{itrial} =cell2mat(Layer2.NS_AC{iROI}{14,2}{itrial});
         q2{iROI} =cell2mat(q{iROI});
         q3 = cell2mat(q2);
         end
        end
    end
    L2_PA_groups{1,2} = {q3};
    L2_PA_Mean{1,2} = nanmean(q3);
    L2_PA_Mean{1,3} = nanstd(q3) 
 
   clear 'y' 'y2' 'y3' 'k' 'k2' 'k3' 'q' 'q2' 'q3'
   
   %STIM 
%Both Layers
     for iROI = 1:length(Stim8s.AC)
         for itrial = 1:length(Stim8s.AC{iROI}{14,2})
         for ipeak = 1:length(Stim8s.AC{iROI}{14,2}{itrial})
           y{iROI}{itrial} = cell2mat(Stim8s.AC{iROI}{14,2}{itrial});
           y2{iROI} =cell2mat(y{iROI});
           y3 = cell2mat(y2);
         end
         end
     end
    PA_Mean{5,2} = nanmean(y3);
    PA_Mean{5,3} = nanstd(y3); 
    PA_groups{5,2} = {y3};

    %Layer 1
     for iROI = 1:length(Layer1.S_AC)
         for itrial = 1:length(Layer1.S_AC{iROI}{14,2})
         for ipeak = 1:length(Layer1.S_AC{iROI}{14,2}{itrial})
          k{iROI}{itrial} =cell2mat(Layer1.S_AC{iROI}{14,2}{itrial});
          k2{iROI} =cell2mat(k{iROI});
          k3 = cell2mat(k2);
         end
         end
    end
    L1_PA_groups{5,2} = {k3};
    L1_PA_Mean{5,2} = nanmean(k3);
    L1_PA_Mean{5,3} = nanstd(k3);
    
    %Layer 2
     for iROI = 1:length(Layer2.S_AC)
        for itrial = 1:length(Layer2.S_AC{iROI}{14,2})
         for ipeak = 1:length(Layer2.S_AC{iROI}{14,2}{itrial})
          q{iROI}{itrial} =cell2mat(Layer2.S_AC{iROI}{14,2}{itrial});
         q2{iROI} =cell2mat(q{iROI});
         q3 = cell2mat(q2);
         end
        end
    end
    L2_PA_groups{5,2} = {q3};
    L2_PA_Mean{5,2} = nanmean(q3);
    L2_PA_Mean{5,3} = nanstd(q3) 
 
    clear 'y' 'y2' 'y3' 'k' 'k2' 'k3' 'q' 'q2' 'q3'


   %%  AP ROIs 
%NO STIM 
%Both Layers
      for iROI = 1:length(Nostim.AP)
         for itrial = 1:length(Nostim.AP{iROI}{14,2})
         for ipeak = 1:length(Nostim.AP{iROI}{14,2}{itrial})
           y{iROI}{itrial} = cell2mat(Nostim.AP{iROI}{14,2}{itrial});
           y2{iROI} =cell2mat(y{iROI});
           y3 = cell2mat(y2);
         end
         end
     end
    
     
    PA_Mean{2,2} = nanmean(y3);
    PA_Mean{2,3} = nanstd(y3); 
    PA_groups{2,2} = {y3};

    %Layer 1
     for iROI = 1:length(Layer1.NS_AP)
         for itrial = 1:length(Layer1.NS_AP{iROI}{14,2})
         for ipeak = 1:length(Layer1.NS_AP{iROI}{14,2}{itrial})
          k{iROI}{itrial} =cell2mat(Layer1.NS_AP{iROI}{14,2}{itrial});
          k2{iROI} =cell2mat(k{iROI});
          k3 = cell2mat(k2);
         end
         end
    end
    L1_PA_groups{2,2} = {k3};
    L1_PA_Mean{2,2} = nanmean(k3);
    L1_PA_Mean{2,3} = nanstd(k3);
    
    %Layer 2
    for iROI = 1:length(Layer2.NS_AP)
        for itrial = 1:length(Layer2.NS_AP{iROI}{14,2})
         for ipeak = 1:length(Layer2.NS_AP{iROI}{14,2}{itrial})
          q{iROI}{itrial} =cell2mat(Layer2.NS_AP{iROI}{14,2}{itrial});
         q2{iROI} =cell2mat(q{iROI});
         q3 = cell2mat(q2);
         end
        end
    end
    L2_PA_groups{2,2} = {q3};
    L2_PA_Mean{2,2} = nanmean(q3);
    L2_PA_Mean{2,3} = nanstd(q3) 
 
   clear 'y' 'y2' 'y3' 'k' 'k2' 'k3' 'q' 'q2' 'q3'
   
   %STIM 
%Both Layers
     for iROI = 1:length(Stim8s.AP)
         for itrial = 1:length(Stim8s.AP{iROI}{14,2})
         for ipeak = 1:length(Stim8s.AP{iROI}{14,2}{itrial})
           y{iROI}{itrial} = cell2mat(Stim8s.AP{iROI}{14,2}{itrial});
           y2{iROI} =cell2mat(y{iROI});
           y3 = cell2mat(y2);
         end
         end
     end
    PA_Mean{6,2} = nanmean(y3);
    PA_Mean{6,3} = nanstd(y3); 
    PA_groups{6,2} = {y3};

    %Layer 1
     for iROI = 1:length(Layer1.S_AP)
         for itrial = 1:length(Layer1.S_AP{iROI}{14,2})
         for ipeak = 1:length(Layer1.S_AP{iROI}{14,2}{itrial})
          k{iROI}{itrial} =cell2mat(Layer1.S_AP{iROI}{14,2}{itrial});
          k2{iROI} =cell2mat(k{iROI});
          k3 = cell2mat(k2);
         end
         end
    end
    L1_PA_groups{6,2} = {k3};
    L1_PA_Mean{6,2} = nanmean(k3);
    L1_PA_Mean{6,3} = nanstd(k3);
    
    %Layer 2
     for iROI = 1:length(Layer2.S_AP)
        for itrial = 1:length(Layer2.S_AP{iROI}{14,2})
         for ipeak = 1:length(Layer2.S_AP{iROI}{14,2}{itrial})
          q{iROI}{itrial} =cell2mat(Layer2.S_AP{iROI}{14,2}{itrial});
         q2{iROI} =cell2mat(q{iROI});
         q3 = cell2mat(q2);
         end
        end
    end
    L2_PA_groups{6,2} = {q3};
    L2_PA_Mean{6,2} = nanmean(q3);
    L2_PA_Mean{6,3} = nanstd(q3) 
 
    clear 'y' 'y2' 'y3' 'k' 'k2' 'k3' 'q' 'q2' 'q3'
   
      %%  BR ROIs 
%NO STIM 
%Both Layers
     for iROI = 1:length(Nostim.BR)
         for itrial = 1:length(Nostim.BR{iROI}{14,2})
         for ipeak = 1:length(Nostim.BR{iROI}{14,2}{itrial})
           y{iROI}{itrial} = cell2mat(Nostim.BR{iROI}{14,2}{itrial});
           y2{iROI} =cell2mat(y{iROI});
           y3 = cell2mat(y2);
         end
         end
     end
    
     
    PA_Mean{3,2} = nanmean(y3);
    PA_Mean{3,3} = nanstd(y3); 
    PA_groups{3,2} = {y3};

    %Layer 1
     for iROI = 1:length(Layer1.NS_BR)
         for itrial = 1:length(Layer1.NS_BR{iROI}{14,2})
         for ipeak = 1:length(Layer1.NS_BR{iROI}{14,2}{itrial})
          k{iROI}{itrial} =cell2mat(Layer1.NS_BR{iROI}{14,2}{itrial});
          k2{iROI} =cell2mat(k{iROI});
          k3 = cell2mat(k2);
         end
         end
    end
    L1_PA_groups{3,2} = {k3};
    L1_PA_Mean{3,2} = nanmean(k3);
    L1_PA_Mean{3,3} = nanstd(k3);
    
    %Layer 2
    for iROI = 1:length(Layer2.NS_BR)
        for itrial = 1:length(Layer2.NS_BR{iROI}{14,2})
         for ipeak = 1:length(Layer2.NS_BR{iROI}{14,2}{itrial})
          q{iROI}{itrial} =cell2mat(Layer2.NS_BR{iROI}{14,2}{itrial});
         q2{iROI} =cell2mat(q{iROI});
         q3 = cell2mat(q2);
         end
        end
    end
    L2_PA_groups{3,2} = {q3};
    L2_PA_Mean{3,2} = nanmean(q3);
    L2_PA_Mean{3,3} = nanstd(q3) 
 
   clear 'y' 'y2' 'y3' 'k' 'k2' 'k3' 'q' 'q2' 'q3'
   
   %STIM 
%Both Layers
     for iROI = 1:length(Stim8s.BR)
         for itrial = 1:length(Stim8s.BR{iROI}{14,2})
         for ipeak = 1:length(Stim8s.BR{iROI}{14,2}{itrial})
           y{iROI}{itrial} = cell2mat(Stim8s.BR{iROI}{14,2}{itrial});
           y2{iROI} =cell2mat(y{iROI});
           y3 = cell2mat(y2);
         end
         end
     end
    PA_Mean{7,2} = nanmean(y3);
    PA_Mean{7,3} = nanstd(y3); 
    PA_groups{7,2} = {y3};

    %Layer 1
     for iROI = 1:length(Layer1.S_BR)
         for itrial = 1:length(Layer1.S_BR{iROI}{14,2})
         for ipeak = 1:length(Layer1.S_BR{iROI}{14,2}{itrial})
          k{iROI}{itrial} =cell2mat(Layer1.S_BR{iROI}{14,2}{itrial});
          k2{iROI} =cell2mat(k{iROI});
          k3 = cell2mat(k2);
         end
         end
    end
    L1_PA_groups{7,2} = {k3};
    L1_PA_Mean{7,2} = nanmean(k3);
    L1_PA_Mean{7,3} = nanstd(k3);
    
    %Layer 2
     for iROI = 1:length(Layer2.S_BR)
        for itrial = 1:length(Layer2.S_BR{iROI}{14,2})
         for ipeak = 1:length(Layer2.S_BR{iROI}{14,2}{itrial})
          q{iROI}{itrial} =cell2mat(Layer2.S_BR{iROI}{14,2}{itrial});
         q2{iROI} =cell2mat(q{iROI});
         q3 = cell2mat(q2);
         end
        end
    end
    L2_PA_groups{7,2} = {q3};
    L2_PA_Mean{7,2} = nanmean(q3);
    L2_PA_Mean{7,3} = nanstd(q3) 
 
    clear 'y' 'y2' 'y3' 'k' 'k2' 'k3' 'q' 'q2' 'q3'

   
     %%  EF ROIs 
%NO STIM 
%Both Layers
     for iROI = 1:length(Nostim.EF)
         for itrial = 1:length(Nostim.EF{iROI}{14,2})
         for ipeak = 1:length(Nostim.EF{iROI}{14,2}{itrial})
           y{iROI}{itrial} = cell2mat(Nostim.EF{iROI}{14,2}{itrial});
           y2{iROI} =cell2mat(y{iROI});
           y3 = cell2mat(y2);
         end
         end
     end
    
     
    PA_Mean{4,2} = nanmean(y3);
    PA_Mean{4,3} = nanstd(y3); 
    PA_groups{4,2} = {y3};

    %Layer 1
     for iROI = 1:length(Layer1.NS_EF)
         for itrial = 1:length(Layer1.NS_EF{iROI}{14,2})
         for ipeak = 1:length(Layer1.NS_EF{iROI}{14,2}{itrial})
          k{iROI}{itrial} =cell2mat(Layer1.NS_EF{iROI}{14,2}{itrial});
          k2{iROI} =cell2mat(k{iROI});
          k3 = cell2mat(k2);
         end
         end
    end
    L1_PA_groups{4,2} = {k3};
    L1_PA_Mean{4,2} = nanmean(k3);
    L1_PA_Mean{4,3} = nanstd(k3);
    
    %Layer 2
    for iROI = 1:length(Layer2.NS_EF)
        for itrial = 1:length(Layer2.NS_EF{iROI}{14,2})
         for ipeak = 1:length(Layer2.NS_EF{iROI}{14,2}{itrial})
          q{iROI}{itrial} =cell2mat(Layer2.NS_EF{iROI}{14,2}{itrial});
         q2{iROI} =cell2mat(q{iROI});
         q3 = cell2mat(q2);
         end
        end
    end
    L2_PA_groups{4,2} = {q3};
    L2_PA_Mean{4,2} = nanmean(q3);
    L2_PA_Mean{4,3} = nanstd(q3) 
 
   clear 'y' 'y2' 'y3' 'k' 'k2' 'k3' 'q' 'q2' 'q3'
   
   %STIM 
%Both Layers
     for iROI = 1:length(Stim8s.EF)
         for itrial = 1:length(Stim8s.EF{iROI}{14,2})
         for ipeak = 1:length(Stim8s.EF{iROI}{14,2}{itrial})
           y{iROI}{itrial} = cell2mat(Stim8s.EF{iROI}{14,2}{itrial});
           y2{iROI} =cell2mat(y{iROI});
           y3 = cell2mat(y2);
         end
         end
     end
    PA_Mean{8,2} = nanmean(y3);
    PA_Mean{8,3} = nanstd(y3); 
    PA_groups{8,2} = {y3};

    %Layer 1
     for iROI = 1:length(Layer1.S_EF)
         for itrial = 1:length(Layer1.S_EF{iROI}{14,2})
         for ipeak = 1:length(Layer1.S_EF{iROI}{14,2}{itrial})
          k{iROI}{itrial} =cell2mat(Layer1.S_EF{iROI}{14,2}{itrial});
          k2{iROI} =cell2mat(k{iROI});
          k3 = cell2mat(k2);
         end
         end
    end
    L1_PA_groups{8,2} = {k3};
    L1_PA_Mean{8,2} = nanmean(k3);
    L1_PA_Mean{8,3} = nanstd(k3);
    
    %Layer 2
     for iROI = 1:length(Layer2.S_EF)
        for itrial = 1:length(Layer2.S_EF{iROI}{14,2})
         for ipeak = 1:length(Layer2.S_EF{iROI}{14,2}{itrial})
          q{iROI}{itrial} =cell2mat(Layer2.S_EF{iROI}{14,2}{itrial});
         q2{iROI} =cell2mat(q{iROI});
         q3 = cell2mat(q2);
         end
        end
    end
    L2_PA_groups{8,2} = {q3};
    L2_PA_Mean{8,2} = nanmean(q3);
    L2_PA_Mean{8,3} = nanstd(q3) 
 
    clear 'y' 'y2' 'y3' 'k' 'k2' 'k3' 'q' 'q2' 'q3'
   
%%  Bar Graphs
%BOTH LAYERS
No_stim = cell2mat(PA_Mean(1:4,2)');
Stim =cell2mat(PA_Mean(5:8,2)');
No_stim_SD = cell2mat(PA_Mean(1:4,3)');
Stim_SD =cell2mat(PA_Mean(5:8,3)');

Both = [No_stim' Stim'];
Both_SD = [No_stim_SD' Stim_SD'];

%Layer 1
L1_No_stim = cell2mat(L1_PA_Mean(1:4,2)');
L1_Stim =cell2mat(L1_PA_Mean(5:8,2)');
L1_No_stim_SD = cell2mat(L1_PA_Mean(1:4,3)');
L1_Stim_SD =cell2mat(L1_PA_Mean(5:8,3)');

L1 = [L1_No_stim' L1_Stim'];
L1_SD = [L1_No_stim_SD' L1_Stim_SD'];

%Layer 2
L2_No_stim = cell2mat(L2_PA_Mean(1:4,2)');
L2_Stim =cell2mat(L2_PA_Mean(5:8,2)');
L2_No_stim_SD = cell2mat(L2_PA_Mean(1:4,3)');
L2_Stim_SD =cell2mat(L2_PA_Mean(5:8,3)');

L2 = [L2_No_stim' L2_Stim'];
L2_SD = [L2_No_stim_SD' L2_Stim_SD'];

%Figures
figure ('name','Peak AUC','numbertitle','off')
subplot (3,1,1)
hold on
barwitherr(Both, Both_SD)
set(gca,'XTickLabel', '')
set(gca,'YLim', [0 2000],'YTick',0:200:2000)
ylabel('Area Under Curve')
title('Both Layers')
legend('NoStim','Stim') 

subplot (3,1,2)
hold on
barwitherr(L1,L1_SD)
set(gca,'XTickLabel', '')
set(gca,'YLim', [0 2000],'YTick',0:200:2000)
ylabel('Area Under Curve')
title('Layer1')
legend('NoStim','Stim') 

subplot (3,1,3)
hold on
barwitherr(L2,L2_SD)
set(gca,'XTickLabel', '')
set(gca,'YLim', [0 2000],'YTick',0:200:2000)
ylabel('Area Under Curve')
title('Layer2')
legend('NoStim','Stim')




%% Statistics- T-test for stim vs. nostim
%Both Layers

%AC ROIs
AC_NS = cell2mat(PA_groups{1,2});
AC_S = cell2mat(PA_groups{5,2});
AC = padcat(AC_NS, AC_S);
AC1= AC(1,:);
AC2= AC(2,:);

[h(1,1), p(1,1)] =ttest(AC1, AC2)

%AP ROIs
AP_NS = cell2mat(PA_groups{2,2});
AP_S = cell2mat(PA_groups{6,2});
AP = padcat(AP_NS, AP_S);
AP1= AP(1,:);
AP2= AP(2,:);

[h(2,1), p(2,1)] =ttest(AP1, AP2)

%BR ROIs
BR_NS = cell2mat(PA_groups{3,2});
BR_S = cell2mat(PA_groups{7,2});
BR = padcat(BR_NS, BR_S);
BR1= BR(1,:);
BR2= BR(2,:);

[h(3,1), p(3,1)] =ttest(BR1, BR2)

%EF ROIs
EF_NS = cell2mat(PA_groups{4,2});
EF_S = cell2mat(PA_groups{8,2});
EF = padcat(EF_NS, EF_S);
EF1= EF(1,:);
EF2= EF(2,:);

[h(4,1), p(4,1)] =ttest(EF1, EF2)

%% Layer 1

%AC ROIs
L1_AC_NS = cell2mat(L1_PA_groups{1,2});
L1_AC_S = cell2mat(L1_PA_groups{5,2});
L1_AC = padcat(L1_AC_NS, L1_AC_S);
L1_AC1= L1_AC(1,:);
L1_AC2= L1_AC(2,:);

[h(1,2), p(1,2)] =ttest(L1_AC1, L1_AC2)

%AP ROIs
L1_AP_NS = cell2mat(L1_PA_groups{2,2});
L1_AP_S = cell2mat(L1_PA_groups{6,2});
L1_AP = padcat(L1_AP_NS, L1_AP_S);
L1_AP1= L1_AP(1,:);
L1_AP2= L1_AP(2,:);

[h(2,2), p(2,2)] =ttest(L1_AP1, L1_AP2)

%BR ROIs
L1_BR_NS = cell2mat(L1_PA_groups{3,2});
L1_BR_S = cell2mat(L1_PA_groups{7,2});
L1_BR = padcat(L1_BR_NS, L1_BR_S);
L1_BR1= L1_BR(1,:);
L1_BR2= L1_BR(2,:);

[h(3,2), p(3,2)] =ttest(L1_BR1, L1_BR2)

%EF ROIs
L1_EF_NS = cell2mat(L1_PA_groups{4,2});
L1_EF_S = cell2mat(L1_PA_groups{8,2});
L1_EF = padcat(L1_EF_NS, L1_EF_S);
L1_EF1= L1_EF(1,:);
L1_EF2= L1_EF(2,:);

[h(4,2), p(4,2)] =ttest(L1_EF1, L1_EF2)

%% Layer 2

%AC ROIs
L2_AC_NS = cell2mat(L2_PA_groups{1,2});
L2_AC_S = cell2mat(L2_PA_groups{5,2});
L2_AC = padcat(L2_AC_NS, L2_AC_S);
L2_AC1= L2_AC(1,:);
L2_AC2= L2_AC(2,:);

[h(1,3), p(1,3)] =ttest(L2_AC1, L2_AC2)

%AP ROIs
L2_AP_NS = cell2mat(L2_PA_groups{2,2});
L2_AP_S = cell2mat(L2_PA_groups{6,2});
L2_AP = padcat(L2_AP_NS, L2_AP_S);
L2_AP1= L2_AP(1,:);
L2_AP2= L2_AP(2,:);

[h(2,3), p(2,3)] =ttest(L2_AP1, L2_AP2)

%BR ROIs
L2_BR_NS = cell2mat(L2_PA_groups{3,2});
L2_BR_S = cell2mat(L2_PA_groups{7,2});
L2_BR = padcat(L2_BR_NS, L2_BR_S);
L2_BR1= L2_BR(1,:);
L2_BR2= L2_BR(2,:);

[h(3,3), p(3,3)] =ttest(L2_BR1, L2_BR2)

%EF ROIs
L2_EF_NS = cell2mat(L2_PA_groups{4,2});
L2_EF_S = cell2mat(L2_PA_groups{8,2});
L2_EF = padcat(L2_EF_NS, L2_EF_S);
L2_EF1= L2_EF(1,:);
L2_EF2= L2_EF(2,:);

[h(4,3), p(4,3)] =ttest(L2_EF1, L2_EF2)

%% Create Output folder- groups
file = 'E:\Data\Two_Photon_Data\GCaMP6\Whisker_Stim\Analysis_Results\matched_ROIs_2013_12_16';
outputFolder = fullfile(file, 'Data_Groups');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

cd(outputFolder)

filename1= 'Peak_AUC_BothL';
save(filename1, 'PA_groups');

filename2 ='Peak_AUC_Layer1';
save(filename2, 'L1_PA_groups');

filename3 = 'Peak_AUC_Layer2';
save(filename3, 'L2_PA_groups');

%% Create Output folder- Means
file = 'E:\Data\Two_Photon_Data\GCaMP6\Whisker_Stim\Analysis_Results\matched_ROIs_2013_12_16';
outputFolder = fullfile(file, 'Means');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

cd(outputFolder)

filename4= 'PA_BothLMean&SD';
save(filename4, 'PA_Mean');

filename5 ='PA_Lay1_Mean&SD';
save(filename5, 'L1_PA_Mean');

filename6 = 'PA_Lay2_Mean&SD';
save(filename6, 'L2_PA_Mean');

stats.hypothesis_test = h;
stats.pvalues = p;
filename7 = 'PA_ttest_results';
save(filename7, 'stats');