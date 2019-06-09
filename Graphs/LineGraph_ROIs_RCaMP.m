%% test ROIs example plot

% Stimulation

%AC ROIs
figure('name','RCaMP_1sec_stim_5sec_trial','numbertitle','off')

 %trace over all trials

    tvec= (1:1:length(ROI1))/11.84; 
    
    hold on
   
    plot(tvec,ROI3, 'k','LineWidth', 2)
        
    plot(tvec,ROI2+600, 'g','LineWidth', 2)
    plot(tvec,ROI1+1400, 'b','LineWidth', 2)
   
    
    plot([-2 -2],[200 700],'LineWidth', 2)
    %stimulation line
    for x = 0:17
        hold on
        plot([(5*x+1) (5*x+2)],[100 100], 'm','LineWidth', 4)
            end
    set(gca,'XLim',[-5 90], 'XTick', -5:5:90)
    set(gca,'yticklabel',[]);
    grid on
    
  legend('Neuron1','Neuron2','Neuron3')

