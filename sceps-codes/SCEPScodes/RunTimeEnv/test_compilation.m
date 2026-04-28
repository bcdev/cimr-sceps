
%=== testing the compilation

E2E_HOME       = getenv( 'E2E_HOME' );
E2E_CODE       = '/home/cjimenez/Work/Studies/SCEPS2OpenSF/SCEPSCodes';
E2E_WORKSPACE  = '/home/cjimenez/Work/DataE/SCEPS2OpenSF/MatlabWorkSpace';

if 0

  exec      = [ E2E_CODE, '/GeoInputs_Extract' ];

  cfile{1}  = [ E2E_CODE, '/General/ConfigFiles/Global_Configuration.xml' ];
  cfile{2}  = [ E2E_CODE, '/SceGenMod/ConfigFiles/GeoInputs_Extract/GeoInputs_Extract_Local_Configuration_Sceps_Test.xml' ];

  input{1}  = [ E2E_HOME, '/InputData/GeoInputData/GeoScenesTest/cimr_geo_tds_svalbard_e2n01km_reftime_201712170700_asc_input_atmosphere_v1-1.nc' ];
  input{2}  = [ E2E_HOME, '/InputData/GeoInputData/GeoScenesTest/cimr_geo_tds_svalbard_e2n01km_reftime_201712170700_asc_input_surface_v1-1.nc' ];

  output    = [ E2E_WORKSPACE, '/L1_Scene_FMS_L2PAD_Svalbard_AM_test/GeoInputs_Extract_Output' ];

  scol      = [ exec, ' ', cfile{1}, ',', cfile{2}, ' ', input{1}, ',', input{2}, ' ', output ];
  disp(scol)
  system( scol )

  clear exec input output cfile

end


if 1

  exec      = [ E2E_CODE, '/Forward_Model' ];

  cfile{1}  = [ E2E_CODE, '/General/ConfigFiles/Global_Configuration.xml' ];
  cfile{2}  = [ E2E_CODE, '/SceGenMod/ConfigFiles/Forward_Model/Forward_Model_Local_Configuration_Sceps_Test.xml' ];

  input     = [ E2E_WORKSPACE, '/L1_Scene_FMS_L2PAD_Svalbard_AM_test/GeoInputs_Extract_Output' ];
  output    = [ E2E_WORKSPACE, '/L1_Scene_FMS_L2PAD_Svalbard_AM_test/Forward_Model_Output' ];

  scol      = [ exec, ' ', cfile{1}, ',', cfile{2}, ' ', input, ' ', output ];
  disp(scol)
  system( scol )

  clear exec input output cfile

end