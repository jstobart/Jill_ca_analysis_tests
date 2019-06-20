%% LINE SCAN ANALYSIS

% load Prairie View Line Scan

MJ= load_prairie_line();


% Unmixing is important
%MJ= MJ.unmix_chs();


%% Line Scan Velocity

lsv001 = LineScanVel([],MJ);
lsv001.process()
lsv001.plot()
lsv001.opt_config()


%% Line Scan Diameter

lsd001 = LineScanDiam([], MJ_diam);
lsd001.process()
lsd001.plot()
lsd001.opt_config()


%% measure calcium from Line scans