%-------------------------------------------------------------------------------
%
% SESSION   session_L1_Scene_Sensor_Simulation_GeoCard2_Day1
%
%    This session takes produced sensor-free scene radiance
%    fields at the CIMR channel center frequencies and transform them
%    with a simulated CIMR antenna and radiometer into observed radiances
%    at the locations specified by a CIMR orbit simulation. The antenna
%    integration is carried out using simulated antenna patterns. The 
%    simulated radiances caare then ready to be inverted by the CIMR L2
%    retrieval algorithm. 
%
%-------------------------------------------------------------------------------
% Project:	  SCEPS
% Package:	  OSS
% Developer:	  Estellus 
% Contact:	  carlos.jimenez@estellus.fr 
% Initiated:	  
%-------------------------------------------------------------------------------


global E2E_HOME
global SCEPS_WORK_SPACE
global SCEPS_CODES_PATH
global LOG

	
global SESSION_ID
SESSION_ID  = 'L1_Scene_Sensor_Simulation_GeoCard2_Day1';

LOG.info( [ SESSION_ID, ' ** Starting session' ] );


%= if session exists, it is saved with the saving time
%  and then removed to start the new session

sfo = [ SCEPS_WORK_SPACE, '/', SESSION_ID ];

if exist( sfo, 'dir' )

  aux  = char(datetime);
  aux  = strrep( aux, ' ', '-');
  aux  = strrep( aux, ':', '-');
  isfo = [ sfo, '.', aux ];
  eval( [ '!mv ', sfo, ' ', isfo])

end


LOG.info( [ SESSION_ID, ' ** Creating session folder ', sfo ]);
mkdir( sfo);


%= starting the simulation



%=== Extracting required geophysical fields from
%    netcdf database into a regular lat-lon grid
%    and storage in .mat files 


module_id  = 'Sensor_Apply_Antenna';

% global and local configuration files

gconf  = [ SCEPS_CODES_PATH, '/General/ConfigFiles/Global_Configuration.xml'];
lconf  = [ SCEPS_CODES_PATH, '/ObsSimMod/ConfigFiles/Sensor_Apply_Antenna/Sensor_Apply_Antenna_Local_Configuration_geocard2_day1.xml'];
configurationParameters = [ gconf, ',' lconf ];


%= input file with TOA-TB data

gSESSION_ID  = 'L1_Scene_Forward_Model_Simulation_GeoCard2_Day1';
gmodule_id   = 'Forward_Model';
SCENE_DATE   = '20161217';
aaa          = '000';

input{1}  = [ SCEPS_WORK_SPACE, '/', gSESSION_ID, '/', gmodule_id, '_Output/', gmodule_id, '_Output_', SCENE_DATE, '_BTS_aa_', aaa  ];


%= input file with orbit data

input{2}  = [ E2E_HOME, '/SCEPSdata/InputData/OrbitData/SCEPS/OSS_CMR_TEST_MPL_ORBSCT_20280101T180001_99999999T999999_0001.mat' ];


inputs = [ input{1}, ',', input{2} ];



%= folder for module outputs

outputs = [ SCEPS_WORK_SPACE, '/', SESSION_ID, '/', module_id, '_Output' ];



if 0

  %= matlab call
  eval([ module_id, '( configurationParameters, inputs, outputs )'])

else

  %= system call to test executable

  cfile  = [ '/obs/cjimenez/Work/DataE/SCEPSOpenSF/SCEPSscd/SCEPScodes/ObsSimMod/Modules/Sensor_Apply_Antenna.mcr' ];  
  isys = [ '!', cfile, ' ',  configurationParameters, ' ', inputs, ' ', outputs ];
  eval( isys )



end



