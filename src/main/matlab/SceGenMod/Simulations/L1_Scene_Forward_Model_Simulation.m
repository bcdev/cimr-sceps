%-------------------------------------------------------------------------------
%
% SIMULATION   L1_Scene_Forward_Model_Simulation
%
%    This simulation converts a scene of geophysical parameters
%    into radiances at the CIMR frequencies. Only the channel centre
%    frequencies depend on the CIMR sensor. No further sensor parameters
%    are used here, and the radiance fields are reproduced at the original
%    geophysical spatial resolution, without any further sensor simulation.   
%
%    These radiance fields will be used later to reproduce a CIMR observation
%    at a given location after transformation by a simulated CIMR antenna
%    and radiometer. 
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

module_id  = 'GeoInputs_Extract';

% global and local configuration files
configurationParameters = [ 'Global_Configuration.xml,', module_id, '_Local_Configuration.xml' ];

%= input file with geofields

input{1}  = [ E2E_HOME, '/SCEPSdata/InputData/GeoInputData/GeoCardScenes/SCEPS_Polar_Scene_1_v2/cimr_sceps_geo_card_devalgo_polarscene_1_20161217_harmonised_v2p0_atmosphere.nc' ];

input{2}  = [ E2E_HOME, '/SCEPSdata/InputData/GeoInputData/GeoCardScenes/SCEPS_Polar_Scene_1_v2/cimr_sceps_geo_card_devalgo_polarscene_1_20161217_harmonised_v2p0_observations.nc' ];

input{3}  = [ E2E_HOME, '/SCEPSdata/InputData/GeoInputData/GeoCardScenes/SCEPS_Polar_Scene_1_v2/cimr_sceps_geo_card_devalgo_polarscene_1_20161217_harmonised_v2p0_surface.nc' ];


inputs = [ input{1}, ',', input{2},',', input{3} ];



%= folder for module outputs

outputs = [ SCEPS_WORK_SPACE, '/', SESSION_ID, '/', module_id, '_Output' ];


if 1
  eval([ module_id, '( configurationParameters, inputs, outputs )'])
else
  disp('Geoinput commented out')
end


%=== Forward modelling from data stored in the .mat files

module_id  = 'Forward_Model';

% global and local configuration files
configurationParameters = [ 'Global_Configuration.xml,', module_id, '_Local_Configuration.xml' ];

% main input and output folders
inputs			= [ SCEPS_WORK_SPACE, '/', SESSION_ID, '/', module_id, '_Input' ];
outputs			= [ SCEPS_WORK_SPACE, '/', SESSION_ID, '/', module_id, '_Output' ];

eval([ module_id, '( configurationParameters, inputs, outputs )'])







   
  









