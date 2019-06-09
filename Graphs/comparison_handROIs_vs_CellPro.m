%% Graphs of All ROIs together-
%Comparing ROI areas
%  01.014.2014  JS

%open CellProfiler data files
DirName = 'E:\Data\Two_Photon_Data\GCaMP6\Whisker_Stim\Analysis_Results\CellProfiler_Results\';
cd(DirName)
flist=dir(fullfile(DirName, '*.mat')); % will give you list of specifed files
for ii=1:length(flist)
     load(flist(ii).name) % load files in workspace
end

load('E:\Data\Two_Photon_Data\GCaMP6\Whisker_Stim\Analysis_Results\activity_ROIs_2013_12_16\Layer1.mat')
load('E:\Data\Two_Photon_Data\GCaMP6\Whisker_Stim\Analysis_Results\activity_ROIs_2013_12_16\Layer2.mat')

%% determine average ROI area (ALL ROIs from each layer)
% CellProfiler ROIs
ColumnInfo = ['L1_NS';'L2_NS';'L1_S ';'L2_S '];
Column =  cellstr(ColumnInfo);
CP_Area_Mean = Column';
CP_Area_groups = Column';

%Nostim
for x = [1,2]
    t = eval(['Layer' num2str(x(:), '%01d'), '_NS']);
    name = ['L' num2str(x(:), '%01d'),'_NS'];
       
    for iexpt = 1:length(t)
            Signal1{iexpt} = cell2mat(t{iexpt}{8,2});
            Signal2= cell2mat(Signal1);            
      

 CP_Area_Mean{2,x} = nanmean(Signal2);
 CP_Area_Mean{3,x} = nanstd(Signal2);
 CP_Area_groups{2,x} = Signal2;
 clear 'Signal2'
  end
    end
%% Stim

for y = [1, 2]
    q = eval(['Layer' num2str(y(:), '%01d'), '_S']);
    name2 = ['L' num2str(y(:), '%01d'),'_S'];
         
    for xexpt = 1:length(q)
            Signal3{xexpt} = cell2mat(q{xexpt}{8,2});
            Signal4 = cell2mat(Signal3);
    
 CP_Area_Mean{2,y+2} = nanmean(Signal4);
 CP_Area_Mean{3,y+2} = nanstd(Signal4);
 CP_Area_groups{2,y+2} = Signal4;
 clear 'Signal4'
    end
end

 %% Hand-clicked ROIs
 ColumnInfo = ['L1_NS';'L2_NS';'L1_S ';'L2_S '];
Column =  cellstr(ColumnInfo);
HC_Area_Mean = Column';
HC_Area_groups = Column';

%Nostim
for z = [1, 2]
    v = eval(['Layer' num2str(z(:), '%01d'), '.NS_AP']);
for iexpt = 1:length(v)
     Signal5{iexpt}= cell2mat(v{1,iexpt}(9,2));
     Signal6 = cell2mat(Signal5);
 
 HC_Area_Mean{2,z} = nanmean(Signal6);
 HC_Area_Mean{3,z} = nanstd(Signal6);
 HC_Area_groups{2,z} = Signal6;
 
end
clear 'Signal6' 'Signal5'
end

%% 
  %Stim
for zz = [1, 2]
    vv = eval(['Layer' num2str(zz(:), '%01d'), '.S_AP']);
for iexpt = 1:length(vv)
     Signal7{iexpt}= cell2mat(vv{1,iexpt}(9,2));
     Signal8 = cell2mat(Signal7);
 
 HC_Area_Mean{2,zz+2} = nanmean(Signal8);
 HC_Area_Mean{3,zz+2} = nanstd(Signal8);
 HC_Area_groups{2,zz+2} = Signal8;
 
end
clear 'Signal8' 'Signal7'
end        
      

%% histogram
figure ('name', 'CellProfiler- ROI area histogram', 'numbertitle', 'off')
subplot (2,2,1)
CP_L1_NS= cell2mat(CP_Area_groups(2,1));
binranges = 0:700;
[CP_L1_NS_counts] = histc(CP_L1_NS,binranges);
bar(binranges,CP_L1_NS_counts, 'b')
set(gca,'Xlim',[0 700])
ylabel('Number of Signals')
xlabel('ROI area')
title('Layer1 Nostim')
hold on 
HC_L1_NS= cell2mat(HC_Area_groups(2,1));
[HC_L1_NS_counts] = histc(HC_L1_NS,binranges);
bar(binranges,HC_L1_NS_counts,'r')
%set(gca, colour, 'r')
legend('CellPro','HandClick') 

subplot (2,2,2)
CP_L2_NS= cell2mat(CP_Area_groups(2,2));
[CP_L2_NS_counts] = histc(CP_L2_NS,binranges)
bar(binranges,CP_L2_NS_counts,'b')
set(gca,'Xlim',[0 700])
ylabel('Number of Signals')
xlabel('ROI area')
title('Layer2 Nostim')
hold on 
HC_L2_NS= cell2mat(HC_Area_groups(2,1));
[HC_L2_NS_counts] = histc(HC_L2_NS,binranges);
bar(binranges,HC_L2_NS_counts,'r')
legend('CellPro','HandClick') 

subplot (2,2,3)
CP_L1_S= cell2mat(CP_Area_groups(2,3));
[CP_L1_S_counts] = histc(CP_L1_S,binranges);
bar(binranges,CP_L1_S_counts,'b')
set(gca,'Xlim',[0 700])
ylabel('Number of Signals')
xlabel('ROI area')
title('Layer1 Stim')
hold on 
HC_L1_S= cell2mat(HC_Area_groups(2,1));
[HC_L1_S_counts] = histc(HC_L1_S,binranges);
bar(binranges,HC_L1_S_counts,'r')
legend('CellPro','HandClick') 

subplot (2,2,4)
CP_L2_S= cell2mat(CP_Area_groups(2,4));
[CP_L2_S_counts] = histc(CP_L2_S,binranges);
bar(binranges,CP_L2_S_counts,'b')
set(gca,'Xlim',[0 700])
ylabel('Number of Signals')
xlabel('ROI area')
title('Layer2 Stim')
hold on 
HC_L2_S= cell2mat(HC_Area_groups(2,1));
[HC_L2_S_counts] = histc(HC_L2_S,binranges);
bar(binranges,HC_L2_S_counts,'r')
legend('CellPro','HandClick') 
 

%%  Bar Graphs
%Figures

CP(1,1:4) = cell2mat(CP_Area_Mean(2,1:4));
CP_SD(1:4) = cell2mat(CP_Area_Mean(3,1:4));
HC(1,1:4) = cell2mat(HC_Area_Mean(2,1:4));
HC_SD(1:4) = cell2mat(HC_Area_Mean(3,1:4));

L1_NS = [CP(1,1), HC(1,1)]; 
L1_NS_SD = [CP_SD(1,1), HC_SD(1,1)]; 

L2_NS = [CP(1,2), HC(1,2)]; 
L2_NS_SD = [CP_SD(1,2), HC_SD(1,2)]; 

L1_S = [CP(1,3), HC(1,3)]; 
L1_S_SD = [CP_SD(1,3), HC_SD(1,3)]; 

L2_S = [CP(1,4), HC(1,4)]; 
L2_S_SD = [CP_SD(1,4), HC_SD(1,4)]; 
%%
figure ('name','Average ROI Area- CellPro vs HandClick','numbertitle','off')
subplot (2,2,1)
hold on
%bar(L1_NS_SD)
barwitherr(L1_NS_SD, L1_NS)
set(gca,'XTickLabel', '')
%set(gca,'YLim', [0 100],'YTick',0:25:100)
ylabel('Average ROI Area')
title('Layer1-NoStim')
legend('CellPro','HandClick')

subplot (2,2,2)
hold on
barwitherr(L2_NS_SD, L2_NS)
set(gca,'XTickLabel', '')
%set(gca,'YLim', [0 100],'YTick',0:25:100)
ylabel('Average ROI Area')
title('Layer2-NoStim')
legend('CellPro','HandClick')


subplot (2,2,3)
hold on
barwitherr(L1_S_SD,L1_S)
set(gca,'XTickLabel', '')
%set(gca,'YLim', [0 100],'YTick',0:25:100)
ylabel('Average ROI Area')
title('Layer1-Stim')
legend('CellPro','HandClick')

subplot (2,2,4)
hold on
barwitherr(L2_S_SD,L2_S)
set(gca,'XTickLabel', '')
%set(gca,'YLim', [0 100],'YTick',0:25:100)
ylabel('Average ROI Area')
title('Layer2-Stim')
legend('CellPro','HandClick')


%% Comparitive Scatterplot

BR = 1:length(binranges);

figure ('name','Counts per Bin- CellPro vs HandClick','numbertitle','off')
subplot (2,2,1)
hold on
scatter(BR,CP_L1_NS_counts, 'b')
scatter(BR, HC_L1_NS_counts,'r')

ylabel('Number of ROIs')
title('Layer1-NoStim')
legend('CellPro','HandClick')

subplot (2,2,2)
hold on
scatter(BR,CP_L2_NS_counts, 'b')
scatter(BR, HC_L2_NS_counts,'r')
xlabel('ROI area')
ylabel('Number of ROIs')
title('Layer2-NoStim')
legend('CellPro','HandClick')

subplot (2,2,3)
hold on
scatter(BR,CP_L1_S_counts, 'b')
scatter(BR, HC_L1_S_counts,'r')
xlabel('ROI area')
ylabel('Number of ROIs')
title('Layer1-Stim')
legend('CellPro','HandClick')

subplot (2,2,4)
hold on
scatter(BR,CP_L2_S_counts, 'b')
scatter(BR, HC_L2_S_counts,'r')
xlabel('ROI area')
ylabel('Number of ROIs')
title('Layer2-Stim')
legend('CellPro','HandClick')





