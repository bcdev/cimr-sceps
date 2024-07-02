
%=== adding pathts to software mscripts, etc
%    with a previous path cleaning to avoid 
%    naming conflicts and demonstarte stand alone
%    execution



clear
restoredefaultpath


%=  to be changed by user
UID              = '/obs/cjimenez';
codeSCEPSpath    = [ UID, '/Work/DataE/SCEPSOpenSF/SCEPSscd/SCEPScodes' ];
dataSCEPSpath    = [ UID, '/Work/DataE/SCEPSOpenSF/SCEPSscd/SCEPSdata' ];
wspaceSCEPSpath  = [ UID, '/Work/DataE/SCEPSOpenSF/MatlabWorkSpace' ];


%=== startup package

cd(codeSCEPSpath);
startup_matlab_SCEPScodes( codeSCEPSpath, dataSCEPSpath, wspaceSCEPSpath  ) 



%=== running the Observation Simulation Module to produce
%    sensor antenna temperatures from the Top-of-Atmosphere brightnes 
%    temperatures corresponding to the geophysical scene with fields
%    stored in SCEPSdata/InputData/GeoInputData/GeoCardScenes/SCEPS_ts2

cd( [ codeSCEPSpath, '/Sessions' ])

session_L1_Scene_Sensor_Simulation_GeoCard2_Day1

