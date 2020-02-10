%% LINE SCAN ANALYSIS

% load Prairie View Line Scan

MJ2= load_prairie_line();


% Unmixing is important
%MJ= MJ.unmix_chs();


%% Line Scan Velocity

lsv002 = LineScanVel([],MJ);
lsv002.process()
lsv002.plot()
lsv002.opt_config()
lsv002.output_data()

%% Line Scan Diameter

lsd001 = LineScanDiam([], MJ2);
lsd001.process()
lsd001.plot()
lsd001.opt_config()
lsd001.output_data()

%% measure calcium from Line scans