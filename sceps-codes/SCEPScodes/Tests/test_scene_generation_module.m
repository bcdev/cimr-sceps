
%=== adding pathts to software mscripts, etc
%    with a previous path cleaning to avoid 
%    naming conflicts and demonstarte stand alone
%    execution


clear
restoredefaultpath


UID              = '/obs/cjimenez';
codeSCEPSpath    = [ UID, '/Work/Studies/SCEPS2OpenSF/SCEPScodes' ];
dataSCEPSpath    = [ UID, '/Work/DataE/SCEPS2OpenSF/SCEPSdata' ];
wspaceSCEPSpath  = [ UID, '/Work/DataE/SCEPS2OpenSF/MatlabWorkSpace' ];


%=== startup package

cd(codeSCEPSpath);
startup_matlab_SCEPScodes( codeSCEPSpath, dataSCEPSpath, wspaceSCEPSpath  ) 



%=== running the Scene Generation Module to produce
%    Top-of-Atmosphere brightness temperatures from
%    a geophysical scene with fields stored in
%    SCEPSdata/InputData/GeoInputData/GeoCardS

cd( [ codeSCEPSpath, '/Sessions' ])

session_L1_Scene_FMS_L2PAD_Svalbard_AM_test


