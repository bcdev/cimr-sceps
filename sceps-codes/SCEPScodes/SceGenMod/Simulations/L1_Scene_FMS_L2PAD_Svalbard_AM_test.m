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



global SESSION_ID

% Reading from system 
E2E_HOME = getenv( 'E2E_HOME' );


%=== Extracting required geophysical fields from
%    netcdf database into a regular lat-lon grid
%    and storage in .mat files 

module_ge  = 'GeoInputs_Extract';

% global and local configuration files
configurationParameters = [ 'Global_Configuration.xml,', module_ge, '_Local_Configuration.xml' ];

%= input file with geofields

input{1}  = [ E2E_HOME, '/InputData/GeoInputData/GeoScenesTest/cimr_geo_tds_svalbard_e2n01km_reftime_201712170700_asc_input_atmosphere_v1-1.nc' ];
input{2}  = [ E2E_HOME, '/InputData/GeoInputData/GeoScenesTest/cimr_geo_tds_svalbard_e2n01km_reftime_201712170700_asc_input_surface_v1-1.nc' ];
inputs = [ input{1}, ',', input{2} ];



%= folder for module outputs

outputs = [ SCEPS_WORK_SPACE, '/', SESSION_ID, '/', module_ge, '_Output' ];

eval([ module_ge, '( configurationParameters, inputs, outputs )'])


%=== Forward modelling from data stored in the .mat files
%    prepared from previous module

module_fm  = 'Forward_Model';

% global and local configuration files
configurationParameters = [ 'Global_Configuration.xml,', module_fm, '_Local_Configuration.xml' ];

% main input and output folders

%  input is output from previous module
inputs			= [ SCEPS_WORK_SPACE, '/', SESSION_ID, '/', module_ge, '_Output' ];
outputs			= [ SCEPS_WORK_SPACE, '/', SESSION_ID, '/', module_fm, '_Output' ];

eval([ module_fm, '( configurationParameters, inputs, outputs )'])







   
  









