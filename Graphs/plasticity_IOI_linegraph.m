%% Line graph of IOI data over time
%For GC603, GC604, GC605
%  5.9.2013  JS
%Note: there is not a complete data set for each point


%% Input the Ratio Data (change in intensity compared to baseline)
%trimmed
T_Day0 = [0,0,0];
T_Day0_M = mean(T_Day0);
T_Day0_SD = std(T_Day0);
T_Day4 = [-12.60878564,-102.4483616,49.26657367];
T_Day4_M = mean(T_Day4);
T_Day4_SD = std(T_Day4);
T_Day7 = [8.705808143,-80.18941341];
T_Day7_M = mean(T_Day7);
T_Day7_SD = std(T_Day7);
T_Day18 = [-11.95434699,-49.82715767,-104.6633874];
T_Day18_M = mean(T_Day18);
T_Day18_SD = std(T_Day18);

%spared
S_Day0 = [0,0,0];
S_Day0_M = mean(S_Day0);
S_Day0_SD = std(S_Day0);
S_Day4 = [-63.79280095,-5.816005099,55.34486817];
S_Day4_M = mean(S_Day4);
S_Day4_SD = std(S_Day4);
S_Day7 = [18.15904162,23.97716573];
S_Day7_M = mean(S_Day7);
S_Day7_SD = std(S_Day7);
S_Day18 = [11.98742377,-28.8786454,-55.64383239];
S_Day18_M = mean(S_Day18);
S_Day18_SD = std(S_Day18);

%Spared-Trimmed
S_T0 =(S_Day0)-(T_Day0);
S_T0_M = mean(S_T0);
S_T0_SD = std(S_T0);
S_T4 =(S_Day4)-(T_Day4);
S_T4_M = mean(S_T4);
S_T4_SD = std(S_T4);
S_T7 =(S_Day7)-(T_Day7);
S_T7_M = mean(S_T7);
S_T7_SD = std(S_T7);
S_T18 =(S_Day18)-(T_Day18);
S_T18_M = mean(S_T18);
S_T18_SD = std(S_T18);


%% Line Graphs

%Spared and trimmed change in amplitude
Sp = [S_Day0_M, S_Day4_M,S_Day7_M,S_Day18_M];
Tr = [T_Day0_M, T_Day4_M,T_Day7_M,T_Day18_M];
time = [0,4,7,18];
Sp_SD = [S_Day0_SD, S_Day4_SD,S_Day7_SD,S_Day18_SD];
Tr_SD = [T_Day0_SD, T_Day4_SD,T_Day7_SD,T_Day18_SD];

figure ('name','IOI Results','numbertitle','off')
subplot (2,1,1)
hold on
errorbar(time,Sp,Sp_SD,'b','LineWidth', 2,'Marker', 'square')
errorbar(time,Tr,Tr_SD,'r','LineWidth', 2,'Marker', 'square')

set(gca,'XLim', [0 20],'XTick',0:1:20,'YLim', [-120 80],'YTick',-120:20:80)
ylabel('Change in Amp (-%DR/R_0)')
xlabel('Time (Days after start of trimming)')
title('Change in Amplitude')
legend('Spared','Trimmed')  

% Spared-Trimmed

S_T = [S_T0_M,S_T4_M,S_T7_M,S_T18_M];
S_T_SD = [S_T0_SD,S_T4_SD,S_T7_SD,S_T18_SD];

subplot (2,1,2)
hold on
errorbar(time,S_T,S_T_SD, 'k','LineWidth', 2,'Marker', 'square')

set(gca,'XLim', [0 20],'XTick',0:1:20,'YLim', [-80 140],'YTick',-80:20:140)
ylabel('Spared-Trimmed (-%DR/R_0)')
xlabel('Time (Days after start of trimming)')
title('Spared-Trimmed')
  


%% Statistics
%Two way ANOVA- 2 factors (time after whisker cutting and spared vs. trimmed)

groups1 = padcat(S_Day0', S_Day4',S_Day7',S_Day18',T_Day0', T_Day4',T_Day7',T_Day18');
groups2 = padcat(S_T0',S_T4',S_T7',S_T18');
groups3 = padcat(S_T0',[S_T4,S_T7,S_T18]');

names = {'D00','D04','D07','D18'};
names2 = [0,4,7,18,0,4,7,18];

[p,anovatab,stats] = anova1(groups2, names); 
[c,m,h,nms] = multcompare(stats,'ctype','bonferroni');

[p2,table2,stats2]= anova1(groups1,names2);
[c2,m2,h2,nms2] = multcompare(stats,'ctype','bonferroni');
%set(gca, 'Xticklabel',)

[h,p] = ttest2(S_T0, [S_T4,S_T7,S_T18])





