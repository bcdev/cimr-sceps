%-------------------------------------------------------------------------------
%
% SIMULATION   L1_Scene_Sensor_Simulation
%
%    This simulation converts a scene of ...  
%
%    These antenna temperatures ... 
%
%-------------------------------------------------------------------------------
% Project:	  SCEPS
% Package:	  OSS
% Developer:	  Estellus 
% Contact:	  carlos.jimenez@estellus.fr 
% Initiated:	  2023-10-01
%-------------------------------------------------------------------------------



global E2E_HOME
global SESSION_ID



%=== Extracting required geophysical fields from
%    netcdf database into a regular lat-lon grid
%    and storage in .mat files 

module_id  = 'Orbit_Geolocation_Extract';

% global and local configuration files
configurationParameters = [ 'Global_Configuration.xml,', module_id, '_Local_Configuration.xml' ];


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

eval([ module_id, '( configurationParameters, inputs, outputs )'])


%=== Forward modelling from data stored in the .mat files


module_id  = 'Sensor_Apply_Antenna';

% global and local configuration files
configurationParameters = [ 'Global_Configuration.xml,', module_id, '_Local_Configuration.xml' ];


% main input and output folders
inputs			= [ SCEPS_WORK_SPACE, '/', SESSION_ID, '/', module_id, '_Input' ];
outputs			= [ SCEPS_WORK_SPACE, '/', SESSION_ID, '/', module_id, '_Output' ];

eval([ module_id, '( configurationParameters, inputs, outputs )'])

