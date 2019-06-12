%% LINE SCAN ANALYSIS

% load Prairie View Line Scan

MJ= load_prairie_line();


% Unmixing is important
MJ= MJ.unmix_chs();


%% Line Scan Velocity

lsv001 = LineScanVel([],test);
lsv001.process()
lsv001.plot()
lsv001.opt_config()


%% Line Scan Diameter

lsd001 = LineScanDiam();
lsd001.process()
lsd001.plot()
lsv001.opt_config()


%% measure calcium from Line scans