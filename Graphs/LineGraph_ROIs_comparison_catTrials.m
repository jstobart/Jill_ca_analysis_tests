%% test ROIs example plot

% Stimulation

%AC ROIs
figure('name','short_trial_108sec__8secstim_example_ROIs','numbertitle','off')

 %trace over all trials

    AC = cell2mat(GC604_spot1D1_shtrial.AC{1,1}{6,2});
    AP1 = cell2mat(GC604_spot1D1_shtrial.AP{1,1}{6,2});
    AP2 = cell2mat(GC604_spot1D1_shtrial.AP{1,2}{6,2});
    AP3 = cell2mat(GC604_spot1D1_shtrial.AP{1,3}{6,2});
    AP4 = cell2mat(GC604_spot1D1_shtrial.AP{1,4}{6,2});
    BR = cell2mat(GC604_spot1D1_shtrial.BR{1,1}{6,2});
    EF = cell2mat(GC604_spot1D1_shtrial.EF{1,1}{6,2});
    tvec= (1:1:length(AC))/GC604_spot1D1_shtrial.AC{1,1}{4,2}; 
    
    hold on
   
    plot(tvec,AC, 'k','LineWidth', 2)
        
    plot(tvec,EF+100, 'g','LineWidth', 2)
    plot(tvec,BR+200, 'b','LineWidth', 2)
    plot(tvec,AP1+300, 'r','LineWidth', 2)
    plot(tvec,AP2+400, 'y','LineWidth', 2)
    plot(tvec,AP3+500, 'm','LineWidth', 2)
    plot(tvec,AP4+600, 'c','LineWidth', 2)
           
   
    plot([-30 -30],[0 50],'LineWidth', 2)
    %stimulation line
    for x = 0:4
        hold on
        plot([(98*x) (98*x+8)],[-50 -50], 'm','LineWidth', 4)
    end
    set(gca,'XTick',0:98:490)
    grid on
    
  
  legend('AC','EF','BR','AP1','AP2','AP3','AP4')

%% Stimulation

%AC ROIs
figure('name','8sec_stim_300sec_trial_example_ROIs','numbertitle','off')

 %trace over all trials

    AC = cell2mat(GC605_spot3_8sec.AC{1,1}{6,2});
    AP4 = cell2mat(GC605_spot3_8sec.AP{1,1}{6,2});
    AP5 = cell2mat(GC605_spot3_8sec.AP{1,2}{6,2});
    BR = cell2mat(GC605_spot3_8sec.BR{1,1}{6,2});
    EF = cell2mat(GC605_spot3_8sec.EF{1,1}{6,2});
    tvec= (1:1:length(AC))/GC605_spot3_8sec.AC{1,1}{4,2}; 
    
    hold on
   
    plot(tvec,AC, 'k','LineWidth', 2)
        
    plot(tvec,EF+100, 'g','LineWidth', 2)
    plot(tvec,BR+200, 'b','LineWidth', 2)
    plot(tvec,AP4+300, 'r','LineWidth', 2)
    plot(tvec,AP5+400, 'y','LineWidth', 2)       
   
    plot([-30 -30],[0 50],'LineWidth', 2)
    %stimulation line
    for x = 0:9
       hold on
        plot([(180*x) (180*x+8)],[-50 -50], 'm','LineWidth', 4)
    end
    set(gca,'XTick',0:180:1800)
    grid on
    
  
  legend('AC','EF','BR','AP4', 'AP5')

  %% 1 sec stim
  figure('name','1sec_stim_300sec_trial_example_ROIs','numbertitle','off')

 %trace over all trials

    AC = cell2mat(GC605_spot3_1sec.AC{1,1}{6,2});
    AP2 = cell2mat(GC605_spot3_1sec.AP{1,1}{6,2});
    AP3 = cell2mat(GC605_spot3_1sec.AP{1,2}{6,2});
    BR = cell2mat(GC605_spot3_1sec.BR{1,1}{6,2});
    EF = cell2mat(GC605_spot3_1sec.EF{1,1}{6,2});
    tvec= (1:1:length(AC))/GC605_spot3_1sec.AC{1,1}{4,2}; 
    
    hold on
   
    plot(tvec,AC, 'k','LineWidth', 2)
        
    plot(tvec,EF+100, 'g','LineWidth', 2)
    plot(tvec,BR+200, 'b','LineWidth', 2)
    plot(tvec,AP2+300, 'r','LineWidth', 2)
    plot(tvec,AP3+400, 'y','LineWidth', 2)       
   
    plot([-30 -30],[0 50],'LineWidth', 2)
    %stimulation line
    for x = 0:9
        hold on
       plot([(180*x) (180*x+1)],[-50 -50], 'm','LineWidth', 4)
    end
    set(gca,'XTick',0:180:1800)
    set(gca,'yticklabel',[]);
    grid on
    
  
  legend('AC','EF','BR','AP2', 'AP3')
