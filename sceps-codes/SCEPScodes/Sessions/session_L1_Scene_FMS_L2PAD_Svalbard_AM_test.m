%-------------------------------------------------------------------------------
%
% SESSION   session_L1_Scene_FMS_SCEPS2_Polar1
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


global SCEPS_WORK_SPACE
global SCEPS_CODES_PATH

	
global SESSION_ID
SESSION_ID  = 'L1_Scene_FMS_L2PAD_Svalbard_AM_test';

LOG = Logger();
LOG.info( [ SESSION_ID, ' ** Starting session' ] );



%= if session exists, it is saved with the saving time 
%  and then removed to start the new session

sfo = [ SCEPS_WORK_SPACE, '/', SESSION_ID ];

if exist( sfo, 'dir' ) & 0

  aux  = char(datetime);
  aux  = strrep( aux, ' ', '-');
  aux  = strrep( aux, ':', '-');
  isfo = [ sfo, '.', aux ];
  eval( [ '!mv ', sfo, ' ', isfo])

end

LOG.info( [ SESSION_ID, ' ** Creating session folder ', sfo ]);
mkdir( sfo);
cd( [ sfo ] );




%= copying global configuration file from SWBforCIMR to session folder


[success] = copyfile( [ SCEPS_CODES_PATH, '/General/ConfigFiles/Global_Configuration.xml'], sfo );
if success == 0
  osfi_error('Problem with session folder, does it exist?')
end



%= copying global local configuration files to session folder
%  NOTE: using LOC_CIMR_PATH as file are there instead of usual
%        stand alone installation wirh files t SOF_CIMR_PATH

[success] = copyfile( [ SCEPS_CODES_PATH, '/SceGenMod/ConfigFiles/GeoInputs_Extract/GeoInputs_Extract_Local_Configuration_Sceps_Test.xml'], [ 'GeoInputs_Extract_Local_Configuration.xml'] );
if success == 0
  osfi_error([ SESSION_ID, ' ** Problem with local configuration file, does it exist?'])
end


[success] = copyfile( [ SCEPS_CODES_PATH, '/SceGenMod/ConfigFiles/Forward_Model/Forward_Model_Local_Configuration_Sceps_Test.xml'], [ 'Forward_Model_Local_Configuration.xml']   );
if success == 0
  osfi_error([ SESSION_ID, ' ** Problem with local configuration file, does it exist?'])
end


%= starting the simulation

eval( SESSION_ID );

LOG.info( [ SESSION_ID, ' ** Finishing session ', SESSION_ID ] );

   

