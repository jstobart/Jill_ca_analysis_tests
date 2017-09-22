clear all
AllData=[];

%% Open image files and process diameters
% run this section multiple times
VesselDiam=ImgGroup.from_files;

VesselDiam=VesselDiam.process();
VesselDiam.plot()

%VesselDiam.children{1,x}.calcDiameter.config.maxRate=5
%close all
%manually save data in the appropriate imaging folder

%% Open image files and process velocity
% run this section multiple times
% VesselLSPIV=ImgGroup.from_files;
% 
% VesselLSPIV = VesselLSPIV.process(false);
% VesselLSPIV.plot()

%close all
%manually save data in the appropriate imaging folder
%% Open image files and process flux
% run this section multiple times
VesselRadon=ImgGroup.from_files;

VesselRadon=VesselRadon.process();
VesselRadon.plot()

% window time
%VesselRadon.children{1,x}.calcVelocity.config.windowTime

% flux adjustment
%VesselRadon.children{1,x}.plot('windows','nPlotSQRT',3)
%VesselRadon.children{1,x}.calcVelocity.config.thresholdProm

r= [3,4,9,11,17,19,22]; %position of number
for i=1:length(r)
    VesselRadon.children{1,r(i)}.calcVelocity.config.windowTime=100;
end

xx= [1,2,9,14,12,13,14,16,18,19,25]; %position of number
for i=1:length(xx)
    VesselRadon.children{1,xx(i)}.calcVelocity.config.thresholdProm=0.01;
end

xx= [4,5,7,8,17,22,23]; %position of number
for i=1:length(xx)
    VesselRadon.children{1,xx(i)}.calcVelocity.config.thresholdProm=0.1;
end
%close all
%manually save data in the appropriate imaging folder
%% Extract diameters and append into one big cell array
% run this section multiple times
for x = 1:length(VesselDiam.children)
    filename{x,1}=VesselDiam.children{1,x}.name; % file name
    %animal name
    expression1 = 'pr..';
    data{x,1} = regexp(filename{x,1},expression1,'match');
    
    %vessel number
    expression2 = 'v..';
    data{x,2} = regexp(filename{x,1},expression2,'match');
    
    %depth
    expression3 = '-(\d+)-';
    depthString{x,1} = regexp(filename{x,1},expression3,'match');
    expression4= '\d+';
    depthString2 = regexp(depthString{x,1},expression4,'match');
    depthString3 =depthString2{1,1}{1,1};
    data{x,3} = str2double(depthString3);
    
    %branch order
    expression5 = '(\d)[s,n,r,t]';
    order{x,1} = regexp(filename{x,1},expression5,'match');
    expression6= '\d+';
    orderString2 = regexp(order{x,1},expression6,'match');
    orderString3 =orderString2{1,1}{1,1};
    data{x,4} = str2double(orderString3);
    
    %diameter mean
    data{x,5}=VesselDiam.children{1,x}.calcDiameter.data.means.raw.diameter;
    
    %diameter standard deviation
    data{x,6}=VesselDiam.children{1,x}.calcDiameter.data.stdevs.raw.diameter;
    
    %velocity mean
    Vmean = abs(VesselLSPIV.children{1,x}.calcVelocity.data.means.raw.velocity);
    data{x,7}=Vmean;
    
    %velocity standard deviation
    data{x,8}=abs(VesselLSPIV.children{1,x}.calcVelocity.data.stdevs.raw.velocity);
    
    % pulsatility index= (Vmax-Vmin)/Vmean (see Bouvy et al. 2014)
    %use the max 95th percentile and min 5th percentile
    Vmax = prctile(abs(VesselRadon.children{1,x}.calcVelocity.data.velocity),95);
    Vmin = prctile(abs(VesselRadon.children{1,x}.calcVelocity.data.velocity),5);
    %PI = Pulsatility index
    data{x,9}=(Vmax-Vmin)/Vmean;
      
    % flux (# of RBCs/s)
    data{x,10} = VesselRadon.children{1,x}.calcVelocity.data.means.raw.flux;
    
    % linear density (# RBC/mm)
    data{x,11} = VesselRadon.children{1,x}.calcVelocity.data.means.raw.lineDensity;
    
    % hematocrit (volume of RBC/ total blood volume)
    %Hct= (#RBCs x Volume)/ (Pi x (radius)^2 x length)  RBC volume = 45um3 (Hawkey et al. Haemtology 1991).
    %Volume_radius= RBC volume / Pi*r^2
    Volume_radius = (45/(pi*(data{x,5}/2)^2))*0.001; %convert to mm
    
    % hematocrit= linear density*volume radius above
    data{x,12}= data{x,11}*Volume_radius;
    
    %Reynold's number= (density*velocity*diameter)/viscosity
    % viscosity calculation (see Pries et al. 1992)
    %X = 1 + (1.7*exp(-0.35*data{x,4}))-(0.6*exp(-0.01*data{x,4})); %account of influence of diameter
    %HCD= -(X/(2-2*X)) + ((X/(2-2*X))^2 + (data{x,11}/(1-X)))^0.5; % discharage hematocrit
    
    %Blood_density= 

    
end

%%
% flux notation

data{:,13}=1;
iv=[1]; % position of bad flux measurements
data{iv,13}=0;

AllData=vertcat(AllData,data);

%% make final cell array and save as CSV
% names={'animal','vessel','depth','branch_order','diameter','diameter_SD','velocity','velocity_SD',...
%     'PI','flux','linear_density','Hc'};
% All_Diam=vertcet(names, AllData);

%xlswrite('E:\Data\Pericyte project\Two-photon data\Max\Baseline_BloodFlow\AnalyzedDiameters.xls',All_Diam)

