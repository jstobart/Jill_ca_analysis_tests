%% shaded error bar with
green=[(0/255) (158/255) (115/255)];
purple=[(204/255) (121/255) (167/255)];

time=time+1;

figure('name', 'Ret ret vs ret plus diamox average- plus SEM')
hold on
%axis off
%xlim([-3 25]);
lineProps.width = 1;
lineProps.edgestyle = ':';

%ylim([-0.2 3]);
lineProps.col = {green};
mseb(time',ret_plus_meanTrace_Cortex',ret_plus_SEMTrace_Cortex',lineProps)
lineProps.col = {purple};
mseb(time',ret_ret_meanTrace_Cortex',ret_ret_SEMTrace_Cortex',lineProps)

% H2=shadedErrorBar(TimeX,(fastAC_mean+0.75),fastAC_SEM)%,green,1);
% H3=shadedErrorBar(TimeX,(OT_RCaMP_mean+2),RC_SEM)%,purple,1);
%rectangle('Position', [0 -0.3 8 4])
%plot([-1 -1],[-5 40], 'k','LineWidth', 1)
