%-------------------------------------------------------------------------------
%
% SESSION   session_L1_Scene_Forward_Model_Simulation_GeoCard2_Day1
%
%    This session converts a polar scene of geophysical parameters
%    into radiances at the CIMR frequencies. Only the channel centre
%    frequencies depend on the CIMR sensor. No further sensor parameters
%    are used here, and the radiance fields are reproduced at the original
%    geophysical spatial resolution, without any further sensor simulation.   
%
%    These radiance fields will be used later to reproduce a CIMR observation
%    at a given location after transformation by a simulated CIMR antenna
%    and radiometer.
%
%    If CIMR frequencies do not change, this session is independent of
%    CIMR instrumental configurations and therefore only need to be
%    run once.
%
%-------------------------------------------------------------------------------
% Project:	  SCEPS
% Package:	  OSS
% Developer:	  Estellus 
% Contact:	  carlos.jimenez@estellus.fr 
% Initiated:	  2023-10-20
%-------------------------------------------------------------------------------


global E2E_HOME
global SCEPS_WORK_SPACE
global SCEPS_CODES_PATH
global LOG

	
global SESSION_ID
SESSION_ID   = 'L1_Scene_Forward_Model_Simulation_GeoCard2_Day1';


LOG.info( [ SESSION_ID, ' ** Starting session' ] );



%= if session exists, it is saved with the saving time 
%  and then removed to start the new session

if 0

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
cd( [ sfo ] );
disp(sfo)

end


%=== Extracting required geophysical fields from
%    netcdf database into a regular lat-lon grid
%    and storage in .mat files 

module_id  = 'GeoInputs_Extract';

%= global and local configuration files
%  using original files instead of copying
%  to simulation folder

gconf  = [ SCEPS_CODES_PATH, '/General/ConfigFiles/Global_Configuration.xml'];
lconf  = [ SCEPS_CODES_PATH, '/SceGenMod/ConfigFiles/GeoInputs_Extract/GeoInputs_Extract_Local_Configuration_geocard2_day1_v2.xml'];
configurationParameters = [ gconf, ',' lconf ];

%= input file with geofields

input{1}  = [ E2E_HOME, '/SCEPSdata/InputData/GeoInputData/GeoCardScenes/SCEPS_Polar_Scene_1_v2/cimr_sceps_geo_card_devalgo_polarscene_1_20161217_harmonised_v2p0_atmosphere.nc' ];

input{2}  = [ E2E_HOME, '/SCEPSdata/InputData/GeoInputData/GeoCardScenes/SCEPS_Polar_Scene_1_v2/cimr_sceps_geo_card_devalgo_polarscene_1_20161217_harmonised_v2p0_surface.nc' ];

input{3}  = [ E2E_HOME, '/SCEPSdata/InputData/GeoInputData/GeoCardScenes/SCEPS_Polar_Scene_1_v2/cimr_sceps_geo_card_devalgo_polarscene_1_20161217_harmonised_v2p0_observations.nc' ];

inputs = [ input{1}, ',', input{2},',', input{3} ];




%= folder for module outputs

outputs = [ SCEPS_WORK_SPACE, '/', SESSION_ID, '/', module_id, '_Output' ];

if 1

  if 0

    %= original matlab call
    eval([ module_id, '( configurationParameters, inputs, outputs )'])

  else

    %= system call to test executable

    cfile  = [ '/obs/cjimenez/Work//DataE/SCEPSOpenSF/SCEPSscd/SCEPScodes/SceGenMod/Modules/GeoInputs_Extract.mcr' ];  
    isys = [ '!', cfile, ' ',  configurationParameters, ' ', inputs, ' ', outputs ];
    eval( isys )

  end

else

  %= recopying module output
  pout = [ SCEPS_WORK_SPACE, '/', SESSION_ID, '_complete/', module_id, '_Output' ]
  eval( [ '!cp -rfp ', pout, ' ', outputs ])


end


%= Output folder from previous simulations
%  is input to Forward_Model

predirinput =  [ outputs ];



module_id  = 'Forward_Model';

gconf  = [ SCEPS_CODES_PATH, '/General/ConfigFiles/Global_Configuration.xml'];
lconf  = [ SCEPS_CODES_PATH, '/SceGenMod/ConfigFiles/Forward_Model/Forward_Model_Local_Configuration_geocard2_day1.xml'], [ 'Forward_Model_Local_Configuration.xml'];
configurationParameters = [ gconf, ',' lconf ];



%= main input and output folders

inputs			= predirinput;
outputs			= [ SCEPS_WORK_SPACE, '/', SESSION_ID, '/', module_id, '_Output' ];

if 1

  %= original matlab call
  eval([ module_id, '( configurationParameters, inputs, outputs )'])

else

  %= system call to test executable

  cfile  = [ '/obs/cjimenez/Work//DataE/SCEPSOpenSF/SCEPSscd/SCEPScodes/SceGenMod/Modules/Forward_Model.mcr' ];  
  isys = [ '!', cfile, ' ',  configurationParameters, ' ', inputs, ' ', outputs ];

  eval( isys )

end





