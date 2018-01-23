
clearvars
close all

%% Load data

% save files names

%traces data
load('E:\Data\Two_Photon_Data\Acta2-RCaMP\Results\FilesforMatlab\Traces_VSMC_01_2018.mat');
%load('D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_2D_longstim_28_04_2017.mat');


VSMC_traces=All_traces;

for x=1:length(VSMC_traces)
    number_pos=regexp(VSMC_traces{x,2},'[0-9]');
    if length(number_pos)<2
        trialname=VSMC_traces{x,2};
        number=str2double(trialname(number_pos));
        VSMC_traces{x,2}=strcat('trial',num2str(number,'%02d'));
    end
end

% trace info- ROIType etc.
% get info for plotting
FrameRate=11.84;
nframes=length(VSMC_traces{1,8});
TimeX(1:nframes) = (1:nframes)/FrameRate;


for iROI=1:length(VSMC_traces)
    %find ROITypes
    N_str= strfind(VSMC_traces{iROI, 1},'V');
    r_str= strfind(VSMC_traces{iROI, 1},'r'); %FLIKA ROIs
    
    if ~isempty(N_str)
        VSMC_traces{iROI,13}='VSMC';
    elseif ~isempty(r_str)
        VSMC_traces{iROI,13}='Dendrite';
    end
    
    % make new unique trial names
    VSMC_traces{iROI,14}=strcat(VSMC_traces{iROI,5},'_',VSMC_traces{iROI,4},'_',VSMC_traces{iROI,2});
    % unique ROI names
    VSMC_traces{iROI,15}=strcat(VSMC_traces{iROI,5},'_',VSMC_traces{iROI,4},'_',VSMC_traces{iROI,2}, '_',VSMC_traces{iROI,1});
end


%
for iROI=1:length(VSMC_traces)
    rc_str(iROI)= ~isempty(strfind(VSMC_traces{iROI,1},'VSMC'));
end

VSMC=VSMC_traces(rc_str',:);

%% Plot traces

figure('name','spontaneous VSMC examples')
hold on
axis off
for ii=[9,8,7]
    tempY1=smooth(VSMC{ii,8},5);
    %tempY1=tempY1(1:stimwindow);
    plot(TimeX,tempY1'+(10*(ii-1)),'r')%'LineWidth',1);
    
end
plot([0 0],[0 5], 'k','LineWidth', 1)
plot([0 10],[0 0], 'k','LineWidth', 1)

% near ROI4
Diameter1=[2.228279576
0.455516555
-0.443117384
-1.311262906
-0.929415842
1.575637952
3.453951314
1.736284724
3.226852497
4.598173769
6.254839827
6.673897694
5.354367187
4.695842836
6.963518839
2.546332509
2.353283807
-1.07456214
-3.232413765
-5.872255889
-2.705593726
-2.045829634
-2.065026097
-1.344788409
-0.854858086
1.177648705
-1.121732223
1.260884644
1.099368939
3.039642843
3.5382007
3.678699829
1.523285864
5.550664864
2.430640272
3.462901795
1.981510797
3.492114227
3.529253007
3.4
3.922787924
6.401150123
2.548907712
2.607640189
2.942740018
3.221897252
3.0962462
5.213476341
6.799394939
3.806504163
5.215232795
2.548912358
-1.483533871
5.583464552
0.954631552
6.667495
-0.506805106
-0.472010588
1.870828173
5.194339821
7.946042509
2.301166378
5.860593034
-1.253397504
1.425715308
5.406366125
0.203523952
0.821499772
1.07647463
-2.757035525
0.226884326
3.223353529
3.587655103
-0.342240233
4.881989324
5.490866433
1.806407666
2.347693358
1.590101282
2.030100449
3.064263867
1.195950212
5.398179098
7.761390401
-0.197387608
5.13878476
1.132466945
-1.164000499
0.133412162
1.942823736
-2.174454112
3.025375604
-4.28033979
4.879192937
-0.282545424
2.426158991
0.344732168
6.757917068
0.139788369
-0.922199976
3.885764847
1.505067068
-2.073868774
0.253081977
3.627121882
5.185109608
4.350969594
3.530248331
2.272327075
3.201799793
2.026038765
0.505707926
3.495130867
0.438639727
4.214846266
-1.433116206
2.249245039
4.270378094
0.115259349
2.368482246
3.234066039
2.885500063
1.64545003
4.425496682
5.534194158
3.957083851
2.918007008
6.056057557
2.388773472
1.513105401
2.339904088
2.673361773
-2.475739878
0.421029185
2.378540036
-0.586296267
2.480311212
5.158990951
3.885737897
1.376881704
-0.740399085
2.33806167
-0.84712086
-0.140522176
-2.757905388
-1.998070813
-1.660446894
-5.605521153
1.353473005
-2.287019116
-0.832006991
1.711544279
-1.952080063
4.245865802
7.569552631
1.018265372
0.559659874
-1.886537377
-4.932648246
-1.659216912
0.227004211
-1.802412987
-0.957595313
1.765517137
4.123162154
-1.2784674
-1.318201364
1.466201109
-1.751591988
-2.758424424
-3.498445919
-5.481546632
-4.392217547
-3.028141415
-6.022188291
-2.52528071
-2.89941378
1.315390015
-0.320798483
-1.521961647
3.627350499
1.862443196
1.162877763
-1.784900955
1.096574412
-1.428878412
-4.422136743
-5.540951022
-4.410767653
-2.069845658
-1.970994006
-2.496249499
2.151133607
1.970472553
2.266199932
0.675473389
4.366012368
-0.25391383
2.97511128
-1.373176516
3.538538515
1.945379423
-0.983667967
-4.564630702
-5.113755383
-7.066245936
-8.551670272
-5.811481186
-7.399820252
-6.330590019
-9.48241851
-2.197924612
-3.06127613];

figure('name','spontaneous Diameter example')
hold on
axis off
tempY3=smooth(Diameter1,5);
TimeX2=TimeX(5:5:end);
plot(TimeX2,tempY3','g')%'LineWidth',1);
plot([0 0],[10 15], 'k','LineWidth', 1)
plot([0 10],[10 10], 'k','LineWidth', 1)


%%  Ensheathing pericytes and endfeet
%traces data
load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\FilesforMatlab\Traces_allMice_Lck_nostim_vs_longstim_12_2017.mat');
%load('D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_2D_longstim_28_04_2017.mat');


EP_traces=All_traces;

for x=1:length(EP_traces)
    number_pos=regexp(EP_traces{x,2},'[0-9]');
    if length(number_pos)<2
        trialname=EP_traces{x,2};
        number=str2double(trialname(number_pos));
        EP_traces{x,2}=strcat('trial',num2str(number,'%02d'));
    end
end


for iROI=1:length(EP_traces)
    rc_str(iROI)= ~isempty(strfind(EP_traces{iROI,5},'ARG2'));
end

ARG2 =EP_traces(rc_str',:);

for iROI=1:length(ARG2)
    ep_str(iROI)= ~isempty(strfind(ARG2{iROI,1},'SM'));
    ef_str(iROI)= ~isempty(strfind(ARG2{iROI,1},'EF'));
end

EnsheathingPericytes= ARG2(ep_str',:);
Endfeet=ARG2(ef_str',:);


%% ensheathing pericyte graphs
FrameRate=11.84;
nframes=length(EnsheathingPericytes{1,9});
TimeX3(1:415) = (1:415)/FrameRate;

% stim vs no stim
for iROI=1:length(EnsheathingPericytes)
    ep_NS_str(iROI)= ~isempty(strfind(EnsheathingPericytes{iROI,6},'Nostim'));
    ep_S_str(iROI)= ~isempty(strfind(EnsheathingPericytes{iROI,6},'Stim'));
end

EP_NS= EnsheathingPericytes(ep_NS_str',:);
EP_S= EnsheathingPericytes(ep_S_str',:);

figure('name','spontaneous EP no stim examples')
hold on
axis off
tempY1=smooth(EP_NS{40,9},5);
tempY2=smooth(EP_NS{41,9},5);
%tempY1=tempY1(1:stimwindow);
plot(TimeX3,tempY1(1:415)','b')%'LineWidth',1);
plot(TimeX3,tempY2(1:415)'+2,'r')%'LineWidth',1);

plot([-1 -1],[-1 0], 'k','LineWidth', 1)
plot([-1 4],[-1 -1], 'k','LineWidth', 1)

allTraces=[];
figure('name','stim EP examples')
hold on
axis off
for ii=1:length(EP_S)
    tempY1=smooth(EP_S{ii,9},3);
    if length(tempY1)<415
        continue
    else
        allTraces=horzcat(allTraces,tempY1(1:415));
        plot(TimeX3,tempY1(1:415)','r')%'LineWidth',1);
    end
end
plot([0 0],[0 5], 'k','LineWidth', 1)
plot([0 10],[0 0], 'k','LineWidth', 1)

% mean trace for ensheathing pericytes during whisker stim
EPmeanStim=mean(allTraces');
EPstdStim=std(allTraces');
plot(TimeX3,smooth(EPmeanStim,5),'k')%'LineWidth',1);


%% shaded error bar with mean and sem during stim
red=[1 0 0];

% SEM calculations
EP_SEMStim=EPstdStim/sqrt(size(EPstdStim,2));


figure('name', 'EP RCaMP means- plus SEM during stim')
hold on
axis off
lineProps.width = 1;
lineProps.edgestyle = ':';

%ylim([-0.2 3]);
lineProps.col = {red};
mseb(TimeX3,smooth(EPmeanStim,5)',smooth(EP_SEMStim,5)',lineProps)

% H2=shadedErrorBar(TimeX,(fastAC_mean+0.75),fastAC_SEM)%,green,1);
% H3=shadedErrorBar(TimeX,(OT_RCaMP_mean+2),RC_SEM)%,purple,1);
rectangle('Position', [10 -0.5 8 1])
plot([-0.1 -0.1],[0 0.3], 'k','LineWidth', 1)


% calculate the onset time of the calcium drop in EP cells

% peak onsets and AUC in the first second after stim for each ROI
for iROI= 1:length(EP_S)
    trace=-(EP_S{iROI,9});
    FrameRate=EP_S{iROI,13};
    nframes=length(trace);
    TimeX(1:nframes) = (1:nframes)/FrameRate;
    BL_time=EP_S{iROI, 8}/EP_S{iROI,13};
    baselineCorrectedTime=TimeX-BL_time;
    %first 1 sec after stim onset
    x1=EP_S{iROI, 8};
    x2=round((BL_time+1)*FrameRate);
    x3= round(FrameRate*(BL_time+10));
    
    % onset time
    Onsets=find_first_onset_time(baselineCorrectedTime(10:end), trace(10:end),2.5,2);
    if isempty(Onsets)
        Onsets=nan(1,1);
    end
    EP_S{iROI, 17}= Onsets;
    % AUC
    EP_S{iROI,18}=trapz(trace(x1:x2));
    EP_S{iROI,19}=trapz(trace(x1:x3));
    
end

MeanActaOnset=nanmean(cell2mat(EP_S(:,17)));

%% Endfeet calcium graphs

% stim vs no stim
for iROI=1:length(Endfeet)
    ef_NS_str(iROI)= ~isempty(strfind(Endfeet{iROI,6},'Nostim'));
    ef_S_str(iROI)= ~isempty(strfind(Endfeet{iROI,6},'Stim'));
end

EF_NS= Endfeet(ef_NS_str',:);
EF_S= Endfeet(ef_S_str',:);

figure('name','spontaneous EF no stim examples- same trial as EP')
hold on
axis off
tempY1=smooth(EF_NS{86,9},5);
tempY2=smooth(EF_NS{87,9},5);
tempY3=smooth(EF_NS{88,9},5);
%tempY1=tempY1(1:stimwindow);
plot(TimeX3,tempY1(1:415)','g')%'LineWidth',1);
plot(TimeX3,tempY2(1:415)'+2,'g')%'LineWidth',1);
plot(TimeX3,tempY3(1:415)'+4,'g')%'LineWidth',1);
plot([-1 -1],[-1 0], 'k','LineWidth', 1)
plot([-1 4],[-1 -1], 'k','LineWidth', 1)


EF_StimTraces=[];
for ii=1:length(EF_S)
    tempY1=smooth(EF_S{ii,9},3);
    if length(tempY1)<415
        continue
    else
        EF_StimTraces=horzcat(EF_StimTraces,tempY1(1:415));
    end
end

% mean trace for ensheathing pericytes during whisker stim
EFmeanStim=mean(EF_StimTraces');
EFstdStim=std(EF_StimTraces');

figure('name','stim EF mean')
hold on
axis off
plot(TimeX3,smooth(EFmeanStim,5),'g')%'LineWidth',1);
plot([-1 -1],[-0.5 0], 'k','LineWidth', 1)
plot([10 18],[0.25 0.25], 'k','LineWidth', 1)


%% calculate the onset time of endfeet so we can separate into fast and delayed

% peak onsets and AUC in the first second after stim for each ROI
for iROI= 1:length(EF_S)
    trace=EF_S{iROI,9};
    FrameRate=EF_S{iROI,13};
    nframes=length(trace);
    TimeX(1:nframes) = (1:nframes)/FrameRate;
    BL_time=EF_S{iROI, 8}/EF_S{iROI,13};
    baselineCorrectedTime=TimeX-BL_time;
    %first 1 sec after stim onset
    x1=EF_S{iROI, 8};
    x2=round((BL_time+1)*FrameRate);
    x3= round(FrameRate*(BL_time+10));
    
    % onset time
    Onsets=find_first_onset_time(baselineCorrectedTime(10:end), trace(10:end),2.5,2);
    if isempty(Onsets)
        Onsets=nan(1,1);
    end
    EF_S{iROI, 17}= Onsets;
    % AUC
    EF_S{iROI,18}=trapz(trace(x1:x2));
    EF_S{iROI,19}=trapz(trace(x1:x3));
    
end

% fast and slow endfeet
for iROI=1:length(EF_S)
    fastIdx(iROI)=~isempty(find(EF_S{iROI,17}>0 && EF_S{iROI,17}<=2));
    slowIdx(iROI)=~isempty(find(EF_S{iROI,17}>2));
end
fastAC=EF_S(fastIdx',:);
slowAC=EF_S(slowIdx',:);

MeanFastOnset=nanmean(cell2mat(fastAC(:,17)));
MeanSlowOnset=nanmean(cell2mat(slowAC(:,17)));

% fast  mean
fastAC_StimTraces=[];
for ii=1:size(fastAC,1)
    tempY1=smooth(fastAC{ii,9},3);
    if length(tempY1)<415
        continue
    else
        fastAC_StimTraces=horzcat(fastAC_StimTraces,tempY1(1:415));
    end
end

% mean trace for ensheathing pericytes during whisker stim
fastACmeanStim=mean(fastAC_StimTraces');
fastACstdStim=std(fastAC_StimTraces');

% slow  mean
slowAC_StimTraces=[];
for ii=1:size(slowAC,1)
    tempY1=smooth(slowAC{ii,9},3);
    if length(tempY1)<415
        continue
    else
        slowAC_StimTraces=horzcat(slowAC_StimTraces,tempY1(1:415));
    end
end

% mean trace for ensheathing pericytes during whisker stim
slowACmeanStim=mean(slowAC_StimTraces');
slowACstdStim=std(slowAC_StimTraces');


% plot mean of fast and slow
figure('name','stim EF fast vs slow')
hold on
axis off
plot(TimeX3,smooth(fastACmeanStim,5),'g')%'LineWidth',1);
plot(TimeX3,smooth(slowACmeanStim,5),'r')%'LineWidth',1);
plot([-1 -1],[-0.5 0], 'k','LineWidth', 1)
plot([10 18],[0.25 0.25], 'k','LineWidth', 1)



red=[1 0 0];
blue=[0 0.8 0.8];
green= [0.3 0.8 0];  

% SEM calculations
EP_SEMStim=EPstdStim/sqrt(size(EPstdStim,2));
slowAC_SEMStim=slowACstdStim/sqrt(size(slowACstdStim,2));
fastAC_SEMStim=fastACstdStim/sqrt(size(fastACstdStim,2));

% plot of mean plus SEM
figure('name', 'EP RCaMP plus EF GCaMP- plus SEM during stim')
hold on
axis off
lineProps.width = 1;
lineProps.edgestyle = ':';

lineProps.col = {red};
mseb(TimeX3,smooth(EPmeanStim,5)',smooth(EP_SEMStim,5)',lineProps)
lineProps.col = {green};
mseb(TimeX3,smooth(slowACmeanStim,5)',smooth(slowAC_SEMStim,5)',lineProps)
lineProps.col = {blue};
mseb(TimeX3,smooth(fastACmeanStim,5)',smooth(fastAC_SEMStim,5)',lineProps)
rectangle('Position', [10 -0.5 8 3])
plot([-1 -1],[0 1], 'k','LineWidth', 1)

%%  Ensheathing pericytes and endfeet pharmacology
%traces data
load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\FilesforMatlab\Traces_pharmacology_Lck_nostim_vs_longstim_12_2017.mat');
%load('D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_2D_longstim_28_04_2017.mat');


EP_pharm_traces=All_traces;

for x=1:length(EP_pharm_traces)
    number_pos=regexp(EP_pharm_traces{x,2},'[0-9]');
    if length(number_pos)<2
        trialname=EP_pharm_traces{x,2};
        number=str2double(trialname(number_pos));
        EP_pharm_traces{x,2}=strcat('trial',num2str(number,'%02d'));
    end
end


for iROI=1:length(EP_pharm_traces)
    rc_str(iROI)= ~isempty(strfind(EP_pharm_traces{iROI,5},'ARG2'));
end

ARG2_pharm =EP_pharm_traces(rc_str',:);

for iROI=1:length(ARG2_pharm)
    atropine_str(iROI)= ~isempty(strfind(ARG2_pharm{iROI,7},'Atropine'));
    prazosin_str(iROI)= ~isempty(strfind(ARG2_pharm{iROI,7},'Prazosin'));
end

Atropine= ARG2_pharm(atropine_str',:);
Prazosin= ARG2_pharm(prazosin_str',:);


for iROI=1:length(Atropine)
    ep_str_Atropine(iROI)= ~isempty(strfind(Atropine{iROI,1},'SM'));
    ef_str_Atropine(iROI)= ~isempty(strfind(Atropine{iROI,1},'EF'));
end

Atropine_EP= Atropine(ep_str_Atropine',:);
Atropine_EF=Atropine(ef_str_Atropine',:);


for iROI=1:length(Prazosin)
    ep_str_Prazosin(iROI)= ~isempty(strfind(Prazosin{iROI,1},'SM'));
    ef_str_Prazosin(iROI)= ~isempty(strfind(Prazosin{iROI,1},'EF'));
end

Prazosin_EP= Prazosin(ep_str_Prazosin',:);
Prazosin_EF=Prazosin(ef_str_Prazosin',:);


%%
