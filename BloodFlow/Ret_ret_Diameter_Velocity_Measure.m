clear all
AllData=[];

%% Information about your images
Settings.MainDir = 'E:\Data\Pericyte_project\Two-photon-data\Ret_ret_Mice\Baseline_BloodFlow';

Settings.AnimalNames = {
    'RG14',...
    'RG16',...
    'RG17',...
    'RG18',...
    };
Settings.ScoreSheetNames = {
    'RG14_Scoresheet_AllTrials.xls',...
    'RG16_Scoresheet_AllTrials.xls',...
    'RG17_Scoresheet_AllTrials.xls',...
    'RG18_Scoresheet_AllTrials.xls',...
    };
Settings.Genotype = {'Ret+'};

channel = struct('blood_plasma',2);

% final data file name
SaveFiles{1,1} = fullfile(Settings.MainDir, 'Results', 'Ret+_diam_velocity_30_08_2017.csv');%'Control_Peaks_3Conds.csv'); % all data

%% Load calibration file
calibration ='E:\matlab\CalibrationFiles\calibration_20x.mat';
CalFile = CalibrationPixelSize.load(calibration);

%% Open image files and process flux
% % run this section multiple times
% VesselRadon=ImgGroup.from_files();
% %VesselRadon.children{1}.calcVelocity.config.windowTime=1000;
% VesselRadon.children{1}.calcVelocity.config.thresholdProm=0.001;
% 
% VesselRadon=VesselRadon.process();
% VesselRadon.plot()
% 
% % window time
% %VesselRadon.children{1}.calcVelocity.config.windowTime
% 
% % flux adjustment
% VesselRadon.children{1}.plot('windows','nPlotSQRT',3)
% VesselRadon.children{1}.calcVelocity.config.thresholdProm
% 



%% Open image files and process diameters
% run this section multiple times
VesselDiam=ImgGroup.from_files;

VesselDiam=VesselDiam.process();
VesselDiam.plot()

%VesselDiam.children{1,x}.calcDiameter.config.maxRate=5
%close all
%manually save data in the appropriate imaging folder

% Open image files and process velocity
% run this section multiple times
% VesselLSPIV=ImgGroup.from_files;
% 
% VesselLSPIV = VesselLSPIV.process(false);
% VesselLSPIV.plot()

%close all
%manually save data in the appropriate imaging folder


%% Extract diameters and append into one big cell array
% run this section multiple times
    filename=VesselDiam.children{1,1}.name; % file name
    %animal name
    expression1 = 'pr..';
    data{1,1} = regexp(filename,expression1,'match');
    
    %vessel number
    expression2 = 'v(\d+)';%'v..';
    data{1,2} = regexp(filename,expression2,'match');
    
    %depth
    expression3 = '-(\d+)-';
    depthString = regexp(filename,expression3,'match');
    expression4= '\d+';
    depthString2 = regexp(depthString,expression4,'match');
    depthString3 =depthString2{1,1}{1,1};
    data{1,3} = str2double(depthString3);
    
    %branch order
    expression5 = '(\d)[s,n,r,t]';
    order = regexp(filename,expression5,'match');
    expression6= '\d+';
    orderString2 = regexp(order,expression6,'match');
    orderString3 =orderString2{1,1}{1,1};
    data{1,4} = str2double(orderString3);
    
    %diameter mean
    data{1,5}=VesselDiam.children{1,1}.calcDiameter.data.means.raw.diameter;
    
    %diameter standard deviation
    data{1,6}=VesselDiam.children{1,1}.calcDiameter.data.stdevs.raw.diameter;
    
    %velocity mean LSPIV
    %Vmean = abs(VesselLSPIV.children{1,1}.calcVelocity.data.means.raw.velocity);
    %data{1,7}=Vmean;
    
    %velocity standard deviation
    %data{1,8}=abs(VesselLSPIV.children{1,1}.calcVelocity.data.stdevs.raw.velocity);
    
     %velocity mean Radon
    Vmean2 = abs(VesselRadon.children{1,1}.calcVelocity.data.means.raw.velocity);
    data{1,9}=Vmean2;
    
    %velocity standard deviation
    data{1,10}=abs(VesselRadon.children{1,1}.calcVelocity.data.stdevs.raw.velocity);
    
    %%
    % pulsatility index= (Vmax-Vmin)/Vmean (see Bouvy et al. 2014)
    %use the max 95th percentile and min 5th percentile
    Vmax = prctile(abs(VesselRadon.children{1,1}.calcVelocity.data.velocity),95);
    Vmin = prctile(abs(VesselRadon.children{1,1}.calcVelocity.data.velocity),5);
    %PI = Pulsatility index
    data{1,11}=(Vmax-Vmin)/Vmean2;
      
    % flux (# of RBCs/s)
    data{1,12} = VesselRadon.children{1,1}.calcVelocity.data.means.raw.flux;
    
    % linear density (# RBC/mm)
    data{1,13} = VesselRadon.children{1,1}.calcVelocity.data.means.raw.lineDensity;
    
    % hematocrit (volume of RBC/ total blood volume)
    %Hct= (#RBCs x Volume)/ (Pi x (radius)^2 x length)  RBC volume = 45um3 (Hawkey et al. Haemtology 1991).
    %Volume_radius= RBC volume / Pi*r^2
    Volume_radius = (45/(pi*(data{1,5}/2)^2))*0.001; %convert to mm
    
    % hematocrit= linear density* RBC volume radius above
    data{1,14}= data{1,13}*Volume_radius;
    
    %Reynold's number= (density*velocity*diameter)/viscosity
    % viscosity calculation (see Pries et al. 1992)
    %X = 1 + (1.7*exp(-0.35*data{x,4}))-(0.6*exp(-0.01*data{x,4})); %account of influence of diameter
    %HCD= -(X/(2-2*X)) + ((X/(2-2*X))^2 + (data{x,11}/(1-X)))^0.5; % discharage hematocrit
    
    %Blood_density= 



%%
% flux notation

data{1,15}=1; % 1 if it's good, 0 if it's bad
data{1,15}=0;
%%
AllData=vertcat(AllData,data);

%% make final cell array and save as CSV
 names={'animal','vessel','depth','branch_order','diameter','diameter_SD','vel_LSIV','vel_LSIV_SD',...
     'vel_radon','vel_radon_SD','PI','flux','linear_density','Hc','flux_note'};
AllData2=vertcat(names, AllData);

cd(fullfile(Settings.MainDir, 'Results'));
% write date to created file
cell2csv(SaveFiles{1,1}, AllData2);


%% Adjust process window time for more precise pulsatility index
%AllData=[];
FolderName='E:\Data\Pericyte project\Two-photon data\Max\Baseline_BloodFlow\Ret_Ret\Prr6';

D = dir([FolderName, '\*.mat']);
Num = length(D(not([D.isdir])));

cd(FolderName);

for ifile= 1:Num
    load(D(ifile).name)
    
    %reprocess velocities if the window time is not 1000
    if VesselRadon.children{1}.calcVelocity.config.windowTime~=1000;
        VesselRadon.children{1}.calcVelocity.config.windowTime=1000;
        VesselRadon=VesselRadon.process();
    end
    
    filename=VesselDiam.children{1,1}.name; % file name
    %animal name
    expression1 = 'pr..';
    data{1,1} = char(regexp(filename,expression1,'match'));
    
    %vessel number
    expression2 = 'v(\d+)';%'v..';
    data{1,2} = char(regexp(filename,expression2,'match'));
    
    %depth
    expression3 = '-(\d+)-';
    depthString = regexp(filename,expression3,'match');
    expression4= '\d+';
    depthString2 = regexp(depthString,expression4,'match');
    depthString3 =depthString2{1,1}{1,1};
    data{1,3} = str2double(depthString3);
    
    %branch order
    expression5 = '(\d)[s,n,r,t]';
    order = regexp(filename,expression5,'match');
    expression6= '\d+';
    orderString2 = regexp(order,expression6,'match');
    orderString3 =orderString2{1,1}{1,1};
    data{1,4} = str2double(orderString3);
    
    %diameter mean
    data{1,5}=VesselDiam.children{1,1}.calcDiameter.data.means.raw.diameter;
    
    %diameter standard deviation
    data{1,6}=VesselDiam.children{1,1}.calcDiameter.data.stdevs.raw.diameter;
    
     %velocity mean Radon
    Vmean2 = abs(VesselRadon.children{1,1}.calcVelocity.data.means.raw.velocity);
    data{1,7}=Vmean2;
    
    %velocity standard deviation
    data{1,8}=abs(VesselRadon.children{1,1}.calcVelocity.data.stdevs.raw.velocity);
    
    %pulsatility index= (Vmax-Vmin)/Vmean (see Bouvy et al. 2014)
    %use the max 95th percentile and min 5th percentile
    Vmax = max(abs(VesselRadon.children{1,1}.calcVelocity.data.velocity));
    Vmin = min(abs(VesselRadon.children{1,1}.calcVelocity.data.velocity));
    %PI = Pulsatility index
    data{1,9}=(Vmax-Vmin)/Vmean2;
      
       
    AllData=vertcat(AllData,data);
end


%% make final cell array and save as CSV
 names={'animal','vessel','depth','branch_order','diameter','diameter_SD',...
     'vel_radon','vel_radon_SD','PI'};
 All_Diam2=vertcat(names, AllData);

cd('E:\Data\Pericyte project\Two-photon data\Max\Baseline_BloodFlow\Results');
% write date to created file
cell2csv('PulsatilityIndex_ret_ret.csv', All_Diam2);



