function Sensor_Apply_Antenna( configurationParameters, inputs, outputs)


global E2E_HOME

%= global variable to stored the dirout of the previous module
%global DIROUT_OGE


%= previous simulation output folder
global ORBITGEO_SIMULATION


%=  global variable that stored the name of the 
% TOA input file read by the previous module
global SCENE_TOA_FILE


idfunction = 'Sensor_Apply_Antenna';



%=== Handling module inputs


%= get log class saved as global variable

global LOG
LOG.setDebugMode (true);

%= initialize command line parsing class

clp = CLP (configurationParameters, inputs, outputs);

%= Get inputs, outputs and configuration files using

conf1 = clp.getConfFile(1);
conf2 = clp.getConfFile(2);

dirin  = clp.getInputFile(1);
dirout = clp.getOutputFile(1);

LOG.info([ idfunction, ' ** Input folder: ', dirin ])
LOG.info([ idfunction, ' ** Output folder: ', dirout ])

%= OD, 20240613: override global variable SCENE_TOA_FILE
%= with symbolic link set in Geolocation_Extract module.
SCENE_TOA_FILE   = [ dirin, '/scene_toa_file.nc' ];
LOG.info([ idfunction, ' ** SCENE_TOA_FILE: ', SCENE_TOA_FILE ])


%= Creating folder input and outputs if not existing already
%  Not in use any more

if 0
if ~exist( dirin, 'dir' )
  LOG.info( [ idfunction, ' ** Creating folder ', dirin ]);
  mkdir( dirin );
end  
end

%= Creating folder for inputs, link to previous module simulation

LOG.info( [ idfunction, ' ** Creating folder ', dirin ]);
mkdir(dirin)

% inputs to forward model from previous geo extraction
% creating a link to a subfolder

%=[ orbitgeo_dirout, orbitgeo_module ] = find_pre_folder( ORBITGEO_SIMULATION );
%=subdirin   = [ dirin, '/', orbitgeo_module, '_Output' ];

%=eval( [ '!ln -s ', orbitgeo_dirout, ' ', subdirin ] );

%= Creating folder for outputs 

% creating folder
LOG.info( [ idfunction, ' ** Creating folder ', dirout ]);
mkdir( dirout );


%= Parse configuration files 

cfm1 = ConFM(conf1);
cfm2 = ConFM(conf2);



%= Read global parameters

LOG.info( [ idfunction, ' ** Reading parameters from global configuration file']);

geodata_version = cfm1.getParameter('geodata_version').getValue;
LOG.info( [ idfunction, ' ** geodata_version has value ', geodata_version ]);

software_version = cfm1.getParameter('software_version').getValue;
LOG.info( [ idfunction, ' ** software_version has value ', software_version ]);

workers_number = cfm1.getParameter('workers_number').getValue;
LOG.info( [ idfunction, ' ** workers_number has value ', num2str(workers_number) ]);


%= Read local parameters

time_length_hours = cfm2.getParameter('time_length_hours').getValue;
LOG.info( [ idfunction, ' ** Length of simulation along the orbit is ', num2str(time_length_hours), ' hours' ]);
if time_length_hours > 0.5
  osfi_error('The length of simulation given exceeds the ones currently tested, check with developer' );
end

orbit_type = cfm2.getParameter('orbit_type').getValue;
if strcmp( orbit_type, 'asc' ) 
  LOG.info( [ idfunction, ' ** Simulation for an ascending orbit' ]);
elseif strcmp( orbit_type, 'des' ) 
  LOG.info( [ idfunction, ' ** Simulation for an descending orbit' ]);
else
  osfi_error('Orbit type not recognised, it has to be asc or des');
end

integration_sampling_step = cfm2.getParameter('integration_sampling_step').getValue;
if strcmp( integration_sampling_step, 'one-step' ) 
  LOG.info( [ idfunction, ' ** Antenna integrations with one step per nominal integration time' ]);
elseif strcmp( integration_sampling_step, 'five-steps' ) 
  LOG.info( [ idfunction, ' ** Antenna integrations with five steps per nominal integration time' ]);
else
  osfi_error('integration_sampling has to be one-step or five-steps');
end


antenna_beam_type = cfm2.getParameter('antenna_beam_type').getValue;
LOG.info( [ idfunction, ' ** Antenna beam selected with integration down to ', antenna_beam_type ]);
if ~strcmp( antenna_beam_type, '40db' ) &  ~strcmp( antenna_beam_type, '60db' )
  osfi_error('Antenna_beam_type not recognised, it has to be 40dB or 60dB');
end


band_selection = cfm2.getParameter('band_selection').getValue;
iband{1} = 'L ';
iband{2} = 'C ';
iband{3} = 'X ';
iband{4} = 'Ku ';
iband{5} = 'Ka';
pband    = [];
for a = 1:5
  if str2num(band_selection( 2*(a-1)+1 ))
    pband = [ pband iband{a}];
  end
end
LOG.info( [ idfunction, ' ** Performing sensor simulation for bands ', pband ]);
if length( band_selection ) ~= 9
  osfi_error('Band_selection seems wrong, it has to be .e.g. 1 1 0 1 1');
end
if length(pband)>2
  pband = strrep( pband, ' ', '-');
else
  pband = strrep( pband, ' ', '');
end


scan_start = cfm2.getParameter('scan_start').getValue;
LOG.info( [ idfunction, ' ** Scan starts at', num2str(scan_start) ]);
if scan_start < 0 
  osfi_error('Scan start cannot be negative');
end

scan_end = cfm2.getParameter('scan_end').getValue;
LOG.info( [ idfunction, ' ** Scan ends at', num2str(scan_end) ]);
if scan_end > 10000
  osfi_error('Scan end seems very large, check with developers');
end


nedt_150K_L = cfm2.getParameter('nedt_150K_L').getValue;
LOG.info( [ idfunction, ' ** NeDT for band L set at ', num2str(nedt_150K_L), ' K' ]);
if nedt_150K_L < 0 & nedt_150K_L > 2
  osfi_error('NeDT for band L is out of range');
end

nedt_150K_C = cfm2.getParameter('nedt_150K_C').getValue;
LOG.info( [ idfunction, ' ** NeDT for band C set at ', num2str(nedt_150K_C), ' K' ]);
if nedt_150K_C < 0 & nedt_150K_C > 2
  osfi_error('NeDT for band C is out of range');
end

nedt_150K_X = cfm2.getParameter('nedt_150K_X').getValue;
LOG.info( [ idfunction, ' ** NeDT for band X set at ', num2str(nedt_150K_X), ' K' ]);
if nedt_150K_X < 0 & nedt_150K_X > 2
  osfi_error('NeDT for band X is out of range');
end

nedt_150K_Ku = cfm2.getParameter('nedt_150K_Ku').getValue;
LOG.info( [ idfunction, ' ** NeDT for band Ku set at ', num2str(nedt_150K_Ku), ' K' ]);
if nedt_150K_Ku < 0 & nedt_150K_Ku > 2
  osfi_error('NeDT for band Ku is out of range');
end

nedt_150K_Ka = cfm2.getParameter('nedt_150K_Ka').getValue;
LOG.info( [ idfunction, ' ** NeDT for band Ka set at ', num2str(nedt_150K_Ka), ' K' ]);
if nedt_150K_Ka < 0 & nedt_150K_Ka > 2
  osfi_error('NeDT for band Ka is out of range');
end




%=== Create a folder form simulation internal tmp files

dirtmp = [ dirout, '/tmp' ];
if ~exist( dirtmp, 'dir')
  mkdir( dirtmp )
end


%=== Folder with executable cimrProject

dircimrProject = [ E2E_HOME, '/SCEPSdata/ExecFiles' ]; 


%=== Output file name

file_sim_output = [ dirout, '/cimr_sceps_l1b_', pband, '_beam_', antenna_beam_type, '_sampling_', integration_sampling_step  ];




%=== Some checks

if ~strcmp( geodata_version,'v1.3') 
  osfi_error('The geodata version given cannot be treated by the software');
end

if ~strcmp( software_version,'v1.3') 
  osfi_error('The software version given do not correspond to this software package');
end


%=== Loading center dates for ascending and descending orbit selection

%=odirout     = subdirin;
odfunction = 'Orbit_Geolocation_Extract';


%=filesave = [ odirout, '/', odfunction, '_Output_date_', orbit_type, '_orbit.asc' ];  
%=LOG.info( [ idfunction, ' ** Loading ascii file ', filesave ] );
%=date_center = load(filesave, '-ascii');

filesave = [ dirin, '/', odfunction, '_Output_date_', orbit_type, '_orbit.asc' ];
LOG.info( [ idfunction, ' ** Loading ascii file ', filesave ] );
date_center = load(filesave, '-ascii');








%=== Create the skeleton files for small 5-step or big one-step times

file_skeleton = [ dirtmp, '/cimr_scene_internal_sample_times_skeleton.mat' ];

if strcmp( integration_sampling_step, 'one-step' ) 

  create_cimr_sample_times_bigstep(date_center,   time_length_hours, file_skeleton);

elseif strcmp( integration_sampling_step, 'five-steps' ) 

  create_cimr_sample_times_smallstep(date_center, time_length_hours, file_skeleton);

end




%=== Reading the TOA TBs from the SceGenMod and converting into
%    an internal format. This is only valid for the geo-physical scenes
%    generated by SCEPS.
% 
%    NOTE: To be replaced by reading a TOA-Tb standard
%    file independent of scene provider and with the wind-related
%    output appended


% Create the scene Tb input file.
%ifn = 'cimr_sceps_toa_card_devalgo_polarscene_1_20161217_v2p0_aa_000.nc';
%ofn = 'cimr_sceps_toa_card_devalgo_polarscene_1_20161217_v2p0_aa_000.card.mat';

LOG.info( [ idfunction, ' ** here 1 ', filesave ] );
file_sim_input  = [ dirtmp, '/cimr_scene_internal_simulator_input.mat' ];
LOG.info( [ idfunction, ' ** here 2 ', filesave ] );

create_cimr_test_scene_from_sceps_toa_tbs( SCENE_TOA_FILE, file_sim_input );

LOG.info( [ idfunction, ' ** here 3 ', filesave ] );





%=== Create cimrProject control file

str.num_threads               = workers_number;
str.scan_start		      = scan_start;
str.scan_end		      = scan_end;
str.nedt_150K_L		      = nedt_150K_L;
str.nedt_150K_C		      = nedt_150K_C;
str.nedt_150K_X		      = nedt_150K_X;
str.nedt_150K_Ku	      = nedt_150K_Ku;
str.nedt_150K_Ka	      = nedt_150K_Ka;
str.antenna_beam_type	      = antenna_beam_type;
str.file_control	      = [ dirtmp, '/cimr_project_control_file.cfg' ];
str.file_skeleton	      = file_skeleton;
str.file_sim_input	      = file_sim_input;
str.file_sim_output           = file_sim_output;
str.do_add_nedt_to_netcdf_ta  = 1;
str.integration_sampling_step = integration_sampling_step;


create_simulator_control_file( str );




%=== Running the integrations by calling the simulator executable
%    with the control file


%= integration to produce internal mat file and netcdf with noise


cd( dircimrProject ) 
LOG.info( [ idfunction, ' ** Starting sensor simulation' ]);
irunsh = [ '!./cimrProject compute_l1b_cfg ', str.file_control, ' ', band_selection  ];
tic
eval( [ irunsh ] )

%= checking if output file has been generate
if ~exist( [ str.file_sim_output, '.mat' ], 'file' )
  osfi_error('Simulation was not completed, problem with the execution of cimrProject');
end


%= creating netcdf file no noise

irunsh = [ '!./cimrProject l1b_mat_to_netcdf ', file_sim_output, '.mat ', file_sim_output, '_no_nedt.nc 1 1 2 0 0 0 0 0' ];
eval( irunsh )


%= regenerating netcdf file with noise

irunsh = [ '!./cimrProject l1b_mat_to_netcdf ', file_sim_output, '.mat ', file_sim_output, '_with_nedt.nc 1 1 1' ];
eval( irunsh )
aux = toc/60;
LOG.info( [ idfunction, ' ** Sensor simulation completed in ', num2str(aux), ' minutes' ]);


%= removing tmp folder

eval(['!rm -r ', dirtmp ])



%= return



%========================================================================================

function create_simulator_control_file( str )

global E2E_HOME

%= creating tmp folder for control

a = 0;

if strcmp( str.integration_sampling_step, 'five-steps')

  a=a+1; txt(a).char = [ 'ifn_bs_hbs                               = ', E2E_HOME, '/SCEPSdata/InputData/SensorData/AntennaPaterns/apats_merged_hbs.ins.ds_001.mat']; 

  a=a+1; txt(a).char = [ 'ifn_bs_fhs                               = ', E2E_HOME, '/SCEPSdata/InputData/SensorData/AntennaPaterns/apats_merged_fhs.ins.ds_003.mat' ];

  a=a+1; txt(a).char = [ 'ifn_apats                                = ', E2E_HOME, '/SCEPSdata/InputData/SensorData/AntennaPaterns/cimr_hbs_apats_all_shifted_rescaled_final_renormalized_flipped_v.refinement1.mat' ];

elseif strcmp( str.integration_sampling_step, 'one-step')

  a=a+1; txt(a).char = [ 'ifn_bs_hbs                               = ', E2E_HOME, '/SCEPSdata/InputData/SensorData/AntennaPaterns/apats_merged_hbs.ave.ds_001.mat']; 

  a=a+1; txt(a).char = [ 'ifn_bs_fhs                               = ', E2E_HOME, '/SCEPSdata/InputData/SensorData/AntennaPaterns/apats_merged_fhs.ave.ds_003.mat' ];

  a=a+1; txt(a).char = [ 'ifn_apats                                = ', E2E_HOME, '/SCEPSdata/InputData/SensorData/AntennaPaterns/cimr_hbs_apats_ave_all_shifted_rescaled_final_renormalized_flipped_v.refinement1.mat' ];

end



if strcmp( str.antenna_beam_type, '60db')

  a=a+1; txt(a).char = [ 'ifn_apats_fhs                            = ', E2E_HOME, '/SCEPSdata/InputData/SensorData/AntennaPaterns/cimr_apat_fhs_all.hv.m60.mat' ];

elseif strcmp( str.antenna_beam_type, '40db')

  a=a+1; txt(a).char = 'ifn_apats_fhs                            = -';

end



a=a+1; txt(a).char = 'aptype                                   = actual';

a=a+1; txt(a).char = [ 'ifn_osf                                  = ', E2E_HOME, '/SCEPSdata/InputData/OrbitData/EOCFI/CMR_TEST_MPL_ORBSCT_20280101T180001_99999999T999999_0001.EOF' ];

a=a+1; txt(a).char = [ 'ifn_skel                                 = ', str.file_skeleton ];

a=a+1; txt(a).char = [ 'ifn_scene                                = ', str.file_sim_input ];

a=a+1; txt(a).char = [ 'scan_start                               = ', num2str( str.scan_start) ]; %0

a=a+1; txt(a).char = [ 'scan_end                                 = ', num2str( str.scan_end) ]; %5000

a=a+1; txt(a).char = [ 'ofn                                      = ', str.file_sim_output, '.mat' ];

a=a+1; txt(a).char = [ 'ofn_nc                                   = ', str.file_sim_output, '.nc' ];

a=a+1; txt(a).char = 'nc_output_extra_vars                     = 1';

a=a+1; txt(a).char = 'skip_done                                = 1';

a=a+1; txt(a).char = 'output_scanline_interval                 = -1';

a=a+1; txt(a).char = 'do_brightspots                           = 0';

a=a+1; txt(a).char = [ 'nedt_150K_L                              = ', num2str(str.nedt_150K_L) ];

a=a+1; txt(a).char = [ 'nedt_150K_C                              = ', num2str(str.nedt_150K_C) ];

a=a+1; txt(a).char = [ 'nedt_150K_X                              = ', num2str(str.nedt_150K_X) ];

a=a+1; txt(a).char = [ 'nedt_150K_Ku                             = ', num2str(str.nedt_150K_Ku) ];

a=a+1; txt(a).char = [ 'nedt_150K_Ka                             = ', num2str(str.nedt_150K_Ka) ];

a=a+1; txt(a).char = 'do_restart_from_output_file                = 0';

a=a+1; txt(a).char = [ 'add_nedt_to_netcdf_ta                    = ', num2str(str.do_add_nedt_to_netcdf_ta) ];

a=a+1; txt(a).char = [ 'num_threads                              = ', num2str( str.num_threads) ];

fileID = fopen( str.file_control,'w');   

for ia = 1:a
  fprintf(fileID,'%s\n',txt(ia).char);
end

fclose(fileID);


return

%========================================================================================


function isok = create_cimr_test_scene_from_sceps_toa_tbs(ifn, ofn)


%= testing that the file is valid

try

  ncread(ifn, 'toa_tbs_L_Hpo');
  isok = 1;

catch

  isok = 0;

end  

di = 1;

s.lon = double(npi2pi(ncread(ifn, 'longitude')));
s.lat = double(ncread(ifn, 'latitude'));
s.lat = s.lat(:,:);

s.lon = s.lon(end:-1:1,1:1:end);
s.lat = s.lat(end:-1:1,1:1:end);

s.gridName = 'EASE2_N01KM';
s.areaElement = ones(size(s.lon));

% Get the devalgo test card grid.
%o = load('../scenes/cimr_test_scene_devalgo1_unfiltered_brightspots.mat');
%s.gx_km = o.scene.gx_km;
%s.gy_km = o.scene.gy_km;
%s.lon   = o.scene.lon;
%s.lat   = o.scene.lat;

s.lsm = ones(size(s.lon));



%%%%%%%%%%

s.eia = ncread(ifn, 'incidence_angle');

neia = length(s.eia);

s.bands.L.ths =  ncread(ifn, 'toa_tbs_L_Hpo');
s.bands.L.tvs =  ncread(ifn, 'toa_tbs_L_Vpo');
s.bands.L.t3s =  ncread(ifn, 'toa_tbs_L_3rd');
s.bands.L.t4s =  ncread(ifn, 'toa_tbs_L_4th');

s.bands.C.ths =  ncread(ifn, 'toa_tbs_C_Hpo');
s.bands.C.tvs =  ncread(ifn, 'toa_tbs_C_Vpo');
s.bands.C.t3s =  ncread(ifn, 'toa_tbs_C_3rd');
s.bands.C.t4s =  ncread(ifn, 'toa_tbs_C_4th');

s.bands.X.ths =  ncread(ifn, 'toa_tbs_X_Hpo');
s.bands.X.tvs =  ncread(ifn, 'toa_tbs_X_Vpo');
s.bands.X.t3s =  ncread(ifn, 'toa_tbs_X_3rd');
s.bands.X.t4s =  ncread(ifn, 'toa_tbs_X_4th');

s.bands.Ku.ths =  ncread(ifn, 'toa_tbs_Ku_Hpo');
s.bands.Ku.tvs =  ncread(ifn, 'toa_tbs_Ku_Vpo');
s.bands.Ku.t3s =  ncread(ifn, 'toa_tbs_Ku_3rd');
s.bands.Ku.t4s =  ncread(ifn, 'toa_tbs_Ku_4th');

s.bands.Ka.ths =  ncread(ifn, 'toa_tbs_Ka_Hpo');
s.bands.Ka.tvs =  ncread(ifn, 'toa_tbs_Ka_Vpo');
s.bands.Ka.t3s =  ncread(ifn, 'toa_tbs_Ka_3rd');
s.bands.Ka.t4s =  ncread(ifn, 'toa_tbs_Ka_4th');

s.bands.L.ths =  s.bands.L.ths(end:-1:1,end:-1:1,:);
s.bands.L.tvs =  s.bands.L.tvs(end:-1:1,end:-1:1,:);
s.bands.L.t3s =  s.bands.L.t3s(end:-1:1,end:-1:1,:);
s.bands.L.t4s =  s.bands.L.t4s(end:-1:1,end:-1:1,:);

s.bands.C.ths =  s.bands.C.ths(end:-1:1,end:-1:1,:);
s.bands.C.tvs =  s.bands.C.tvs(end:-1:1,end:-1:1,:);
s.bands.C.t3s =  s.bands.C.t3s(end:-1:1,end:-1:1,:);
s.bands.C.t4s =  s.bands.C.t4s(end:-1:1,end:-1:1,:);

s.bands.X.ths =  s.bands.X.ths(end:-1:1,end:-1:1,:);
s.bands.X.tvs =  s.bands.X.tvs(end:-1:1,end:-1:1,:);
s.bands.X.t3s =  s.bands.X.t3s(end:-1:1,end:-1:1,:);
s.bands.X.t4s =  s.bands.X.t4s(end:-1:1,end:-1:1,:);

s.bands.Ku.ths =  s.bands.Ku.ths(end:-1:1,end:-1:1,:);
s.bands.Ku.tvs =  s.bands.Ku.tvs(end:-1:1,end:-1:1,:);
s.bands.Ku.t3s =  s.bands.Ku.t3s(end:-1:1,end:-1:1,:);
s.bands.Ku.t4s =  s.bands.Ku.t4s(end:-1:1,end:-1:1,:);

s.bands.Ka.ths =  s.bands.Ka.ths(end:-1:1,end:-1:1,:);
s.bands.Ka.tvs =  s.bands.Ka.tvs(end:-1:1,end:-1:1,:);
s.bands.Ka.t3s =  s.bands.Ka.t3s(end:-1:1,end:-1:1,:);
s.bands.Ka.t4s =  s.bands.Ka.t4s(end:-1:1,end:-1:1,:);

%s.bands.L.t3s(find(s.bands.L.t3s < 0.0)) = 0.0;
%s.bands.L.t4s(find(s.bands.L.t4s < 0.0)) = 0.0;

%s.bands.C.t3s(find(s.bands.C.t3s < 0.0)) = 0.0;
%s.bands.C.t4s(find(s.bands.C.t4s < 0.0)) = 0.0;

%s.bands.X.t3s(find(s.bands.X.t3s < 0.0)) = 0.0;
%s.bands.X.t4s(find(s.bands.X.t4s < 0.0)) = 0.0;

%s.bands.Ku.t3s(find(s.bands.Ku.t3s < 0.0)) = 0.0;
%s.bands.Ku.t4s(find(s.bands.Ku.t4s < 0.0)) = 0.0;

%s.bands.Ka.t3s(find(s.bands.Ka.t3s < 0.0)) = 0.0;
%s.bands.Ka.t4s(find(s.bands.Ka.t4s < 0.0)) = 0.0;

%s.bands.L.t3s(:) = 0.0;
%s.bands.L.t4s(:) = 0.0;
%s.bands.C.t3s(:) = 0.0;
%s.bands.C.t4s(:) = 0.0;
%s.bands.X.t3s(:) = 0.0;
%s.bands.X.t4s(:) = 0.0;
%s.bands.Ku.t3s(:) = 0.0;
%s.bands.Ku.t4s(:) = 0.0;
%s.bands.Ka.t3s(:) = 0.0;
%s.bands.Ka.t4s(:) = 0.0;


s.bands.L.ths_default = mean(s.bands.L.ths(:));
s.bands.L.tvs_default = mean(s.bands.L.tvs(:));
s.bands.L.t3s_default = 0.0;
s.bands.L.t4s_default = 0.0;

s.bands.C.ths_default = mean(s.bands.C.ths(:));
s.bands.C.tvs_default = mean(s.bands.C.tvs(:));
s.bands.C.t3s_default = 0.0;
s.bands.C.t4s_default = 0.0;

s.bands.X.ths_default = mean(s.bands.X.ths(:));
s.bands.X.tvs_default = mean(s.bands.X.tvs(:));
s.bands.X.t3s_default = 0.0;
s.bands.X.t4s_default = 0.0;

s.bands.Ku.ths_default = mean(s.bands.Ku.ths(:));
s.bands.Ku.tvs_default = mean(s.bands.Ku.tvs(:));
s.bands.Ku.t3s_default = 0.0;
s.bands.Ku.t4s_default = 0.0;

s.bands.Ka.ths_default = mean(s.bands.Ka.ths(:));
s.bands.Ka.tvs_default = mean(s.bands.Ka.tvs(:));
s.bands.Ka.t3s_default = 0.0;
s.bands.Ka.t4s_default = 0.0;

clear gg;
gg.scene = s;

gg = adjust_card(gg);

eval(['save -V7.3 -nocompression ' ofn ' -struct gg']);

return





%========================================================================================


function [s] = adjust_card(s)


vb = {'L','C','X','Ku','Ka'};
vv = {'ths','tvs','t3s','t4s'};

%%%

is = 5;
for it=1:(is-1)
  s.scene.lsm(it,:)       = s.scene.lsm(is,:);
  s.scene.lsm(it,:)       = s.scene.lsm(is,:);
end

is = size(s.scene.lsm,2)-11;
for it = (is+1):size(s.scene.lsm,2)
  s.scene.lsm(it,:) = s.scene.lsm(is,:);
  s.scene.lsm(it,:) = s.scene.lsm(is,:);
end



is = 5;
for it=1:(is-1)
  s.scene.lsm(:,it)       = s.scene.lsm(:,is);
  s.scene.lsm(:,it)       = s.scene.lsm(:,is);
end

is = size(s.scene.lsm,2)-11;
for it = (is+1):size(s.scene.lsm,2)
  s.scene.lsm(:,it) = s.scene.lsm(:,is);
  s.scene.lsm(:,it) = s.scene.lsm(:,is);
end


%==

for iv = 1:4;

  tvar = vv{iv};

  for i=1:5;

    bn = vb{i};

%s.scene.bands.(bn).(tvar)(1,:,:)        = s.scene.bands.(bn).(tvar)(2,:,:);
%s.scene.bands.(bn).(tvar)(1,:,:)        = s.scene.bands.(bn).(tvar)(2,:,:);
%s.scene.bands.(bn).(tvar)(end,:,:)      = s.scene.bands.(bn).(tvar)(end-1,:,:);
%s.scene.bands.(bn).(tvar)(end,:,:)      = s.scene.bands.(bn).(tvar)(end-1,:,:);

    is = 12;

    for it=1:(is-1)
      s.scene.bands.(bn).(tvar)(it,:,:)       = s.scene.bands.(bn).(tvar)(is,:,:);
      s.scene.bands.(bn).(tvar)(it,:,:)       = s.scene.bands.(bn).(tvar)(is,:,:);
    end

    is = size(s.scene.bands.(bn).(tvar),2)-4;

    for it = (is+1):size(s.scene.bands.(bn).(tvar),2)
      s.scene.bands.(bn).(tvar)(it,:,:) = s.scene.bands.(bn).(tvar)(is,:,:);
      s.scene.bands.(bn).(tvar)(it,:,:) = s.scene.bands.(bn).(tvar)(is,:,:);
    end

    is = 12;

    for it=1:(is-1)
      s.scene.bands.(bn).(tvar)(:,it,:)       = s.scene.bands.(bn).(tvar)(:,is,:);
      s.scene.bands.(bn).(tvar)(:,it,:)       = s.scene.bands.(bn).(tvar)(:,is,:);
    end

    is = size(s.scene.bands.(bn).(tvar),2)-4;


    for it = (is+1):size(s.scene.bands.(bn).(tvar),2)
      s.scene.bands.(bn).(tvar)(:,it,:) = s.scene.bands.(bn).(tvar)(:,is,:);
      s.scene.bands.(bn).(tvar)(:,it,:) = s.scene.bands.(bn).(tvar)(:,is,:);
    end

  end

end



for i=1:5;

  bn = vb{i};

  s.scene.bands.(bn).t3s(find(s.scene.bands.(bn).t3s < -200)) = 0.0;
  s.scene.bands.(bn).t4s(find(s.scene.bands.(bn).t4s < -200)) = 0.0;

  %s.scene.bands.(bn).t3s(:) = 0.0;
  %s.scene.bands.(bn).t4s(:) = 0.0;

end


return





%========================================================================================


function [] = create_cimr_sample_times_bigstep(date_center, time_length_hours, ofn)

%==

bv(1).bn = 'L';
bv(2).bn = 'C';
bv(3).bn = 'X';
bv(4).bn = 'Ku';
bv(5).bn = 'Ka';
bv(1).nh = 1;
bv(2).nh = 4;
bv(3).nh = 4;
bv(4).nh = 8;
bv(5).nh = 8;

bv(1).pointings = [ 0.048,0.016 ];
bv(2).pointings = [ -0.0072,0.0192; -0.0016,-0.0208; 0.004,0.0352; 0.0096,-0.0368];
bv(3).pointings = [ -0.0072222,0.019444; -0.0016667,-0.021667; 0.0044444,0.036111; 0.01,-0.037778 ];
bv(4).pointings = [ -0.0081416,0; -0.0053097,-0.0074336; -0.0024779,0.0063717; 0.00035398,-0.0010619; 0.0028319,-0.0084956; 0.0060177,0.0099115; 0.0084956,0.0028319; 0.011327,-0.0046018];
bv(5).pointings = [ -0.0082927,0; -0.0053659,-0.0078049; -0.0026829,0.0065854; 0.0002439,-0.00097561; 0.0029268,-0.0085366; 0.0060976,0.010244; 0.0087805,0.0026829; 0.011463,-0.004878];

for ib=1:5
txi   = -bv(ib).pointings(:,2);
teta  = -bv(ib).pointings(:,1);
tzeta = -sqrt(1 - txi.*txi - teta.*teta);
[bv(ib).uo, bv(ib).vo] = xe_to_uv(txi, teta, tzeta, 46.886);
end

%==

cfg.L.num_horns                              = bv(1).nh;
cfg.C.num_horns                              = bv(2).nh;
cfg.X.num_horns                              = bv(3).nh;
cfg.Ku.num_horns                             = bv(4).nh;
cfg.Ka.num_horns                             = bv(5).nh;

% Big time steps.
cfg.L.big_time_step_seconds                  = 0.0556;
cfg.C.big_time_step_seconds                  = 0.014;
cfg.X.big_time_step_seconds                  = 0.0137;
cfg.Ku.big_time_step_seconds                 = 0.005;
cfg.Ka.big_time_step_seconds                 = 0.0037;

cfg.L.num_small_time_steps                   = int32(1);
cfg.C.num_small_time_steps                   = int32(1);
cfg.X.num_small_time_steps                   = int32(1);
cfg.Ku.num_small_time_steps                  = int32(1);
cfg.Ka.num_small_time_steps                  = int32(1);

cfg.L.small_time_step_seconds                = cfg.L.big_time_step_seconds/double(cfg.L.num_small_time_steps);
cfg.C.small_time_step_seconds                = cfg.C.big_time_step_seconds/double(cfg.C.num_small_time_steps);
cfg.X.small_time_step_seconds                = cfg.X.big_time_step_seconds/double(cfg.X.num_small_time_steps);
cfg.Ku.small_time_step_seconds               = cfg.Ku.big_time_step_seconds/double(cfg.Ku.num_small_time_steps);
cfg.Ka.small_time_step_seconds               = cfg.Ka.big_time_step_seconds/double(cfg.Ka.num_small_time_steps);

cfg.reflector_off_boresight_angle = 46.886;
cfg.reflector_rpm                 = -7.8;

%==

clear o;

cfg.date_start                    = date_center-time_length_hours/2/24;
cfg.date_end                      = date_center+time_length_hours/2/24;
cfg.reflector_time_az0            = cfg.date_start;
cfg.reflector_az0_deg             = 90.0;

cfg.ddate = cfg.date_end - cfg.date_start;

% Time per scanline in seconds.
cfg.dtscan_sec = 1/abs(cfg.reflector_rpm)*60;

% Loop over scanlines.
cfg.nscans = floor(cfg.ddate*86400/cfg.dtscan_sec);

% Store config information.
o.config = cfg;

% Loop over the bands.
for ib = 1:length(bv);

% Number of horns.
nh    = bv(ib).nh;

% Band name.
bname = bv(ib).bn;

% Compute the number of samples per scanline for this band. In general there is a small
% excess jump in sample time between the last sample of a scan and the first in the next scan.
nsamples = floor(cfg.dtscan_sec / cfg.(bname).small_time_step_seconds);

% Compute the array of sample times for all horns in this band.
dtscan_sec       = cfg.dtscan_sec;
dtsample_sec     = o.config.(bname).small_time_step_seconds;
[gsample, gscan] = ndgrid(1:nsamples, 1:cfg.nscans);
dtime_msec       = 1000*(gsample-1)*dtsample_sec + 1000*(gscan-1)*dtscan_sec;
  
for ih = 1:nh;

hname = sprintf('%s_Horn_%02d', bv(ib).bn, ih);

dims(1) = nh;
dims(2) = nsamples;
dims(3) = cfg.nscans;

o.bands.(bname).uo = bv(ib).uo;
o.bands.(bname).vo = bv(ib).vo;

o.bands.(bname).time_matlab            = cfg.date_start + dtime_msec/1000/86400;
o.bands.(bname).time_from_start_msec   = dtime_msec;

% The reflector azimuth is simply a linear function of time.
o.bands.(bname).razd        = wrapTo360(cfg.reflector_az0_deg + o.bands.(bname).time_from_start_msec/1000/60*360*cfg.reflector_rpm);

% Convert reflector azimuth to L1B standard scan angle (0-180=for; 180-360=aft).
o.bands.(bname).sang        = wrapTo360(90 - o.bands.(bname).razd);

o.bands.(bname).is   = int32(zeros(dims(2:3)));
gind_aft = find(o.bands.(bname).sang > 180);
o.bands.(bname).is(gind_aft) = 1;

% Create a prototype array to be filled in by the integrator.
o.bands.(bname).lat  = nan(dims);

%o.bands.(bname).lon  = nan(dims);
%o.bands.(bname).eia  = nan(dims);
%o.bands.(bname).eaa  = nan(dims);
%o.bands.(bname).pra  = nan(dims);

%o.bands.(bname).ta_x = nan(dims);
%o.bands.(bname).ta_y = nan(dims);
%o.bands.(bname).ta_3 = nan(dims);
%o.bands.(bname).ta_4 = nan(dims);

fprintf('%2d %2d %2s %10s\n', ib, ih, bname, hname);

end

end

eval(['save -V6 ' ofn ' -struct o']);

return






%========================================================================================

function [] = create_cimr_sample_times_smallstep(date_center, time_length_hours, ofn)

%%%%%

bv(1).bn = 'L';
bv(2).bn = 'C';
bv(3).bn = 'X';
bv(4).bn = 'Ku';
bv(5).bn = 'Ka';
bv(1).nh = 1;
bv(2).nh = 4;
bv(3).nh = 4;
bv(4).nh = 8;
bv(5).nh = 8;

bv(1).pointings = [ 0.048,0.016 ];
bv(2).pointings = [ -0.0072,0.0192; -0.0016,-0.0208; 0.004,0.0352; 0.0096,-0.0368];
bv(3).pointings = [ -0.0072222,0.019444; -0.0016667,-0.021667; 0.0044444,0.036111; 0.01,-0.037778 ];
bv(4).pointings = [ -0.0081416,0; -0.0053097,-0.0074336; -0.0024779,0.0063717; 0.00035398,-0.0010619; 0.0028319,-0.0084956; 0.0060177,0.0099115; 0.0084956,0.0028319; 0.011327,-0.0046018];
bv(5).pointings = [ -0.0082927,0; -0.0053659,-0.0078049; -0.0026829,0.0065854; 0.0002439,-0.00097561; 0.0029268,-0.0085366; 0.0060976,0.010244; 0.0087805,0.0026829; 0.011463,-0.004878];

for ib=1:5
txi   = -bv(ib).pointings(:,2);
teta  = -bv(ib).pointings(:,1);
tzeta = -sqrt(1 - txi.*txi - teta.*teta);
[bv(ib).uo, bv(ib).vo] = xe_to_uv(txi, teta, tzeta, 46.886);
end

%%%%

cfg.L.num_horns                              = bv(1).nh;
cfg.C.num_horns                              = bv(2).nh;
cfg.X.num_horns                              = bv(3).nh;
cfg.Ku.num_horns                             = bv(4).nh;
cfg.Ka.num_horns                             = bv(5).nh;

% Big time steps.
cfg.L.big_time_step_seconds                  = 0.0556;
cfg.C.big_time_step_seconds                  = 0.014;
cfg.X.big_time_step_seconds                  = 0.0137;
cfg.Ku.big_time_step_seconds                 = 0.005;
cfg.Ka.big_time_step_seconds                 = 0.0037;

cfg.L.num_small_time_steps                   = int32(5);
cfg.C.num_small_time_steps                   = int32(5);
cfg.X.num_small_time_steps                   = int32(5);
cfg.Ku.num_small_time_steps                  = int32(5);
cfg.Ka.num_small_time_steps                  = int32(5);

cfg.L.small_time_step_seconds                = cfg.L.big_time_step_seconds/double(cfg.L.num_small_time_steps);
cfg.C.small_time_step_seconds                = cfg.C.big_time_step_seconds/double(cfg.C.num_small_time_steps);
cfg.X.small_time_step_seconds                = cfg.X.big_time_step_seconds/double(cfg.X.num_small_time_steps);
cfg.Ku.small_time_step_seconds               = cfg.Ku.big_time_step_seconds/double(cfg.Ku.num_small_time_steps);
cfg.Ka.small_time_step_seconds               = cfg.Ka.big_time_step_seconds/double(cfg.Ka.num_small_time_steps);

cfg.reflector_off_boresight_angle = 46.886;
cfg.reflector_rpm                 = -7.8;

%%%%%

clear o;

cfg.date_start                    = date_center-time_length_hours/2/24;
cfg.date_end                      = date_center+time_length_hours/2/24;
cfg.reflector_time_az0            = cfg.date_start;
cfg.reflector_az0_deg             = 90.0;

cfg.ddate = cfg.date_end - cfg.date_start;

% Time per scanline in seconds.
cfg.dtscan_sec = 1/abs(cfg.reflector_rpm)*60;

% Loop over scanlines.
cfg.nscans = floor(cfg.ddate*86400/cfg.dtscan_sec);

% Store config information.
o.config = cfg;

% Loop over the bands.
for ib = 1:length(bv);

% Number of horns.
nh    = bv(ib).nh;

% Band name.
bname = bv(ib).bn;

% Compute the number of samples per scanline for this band. In general there is a small
% excess jump in sample time between the last sample of a scan and the first in the next scan.
nsamples = floor(cfg.dtscan_sec / cfg.(bname).small_time_step_seconds);

% Compute the array of sample times for all horns in this band.
dtscan_sec       = cfg.dtscan_sec;
dtsample_sec     = o.config.(bname).small_time_step_seconds;
[gsample, gscan] = ndgrid(1:nsamples, 1:cfg.nscans);
dtime_msec       = 1000*(gsample-1)*dtsample_sec + 1000*(gscan-1)*dtscan_sec;
  
for ih = 1:nh;

hname = sprintf('%s_Horn_%02d', bv(ib).bn, ih);

dims(1) = nh;
dims(2) = nsamples;
dims(3) = cfg.nscans;

o.bands.(bname).uo = bv(ib).uo;
o.bands.(bname).vo = bv(ib).vo;

o.bands.(bname).time_matlab            = cfg.date_start + dtime_msec/1000/86400;
o.bands.(bname).time_from_start_msec   = dtime_msec;

% The reflector azimuth is simply a linear function of time.
o.bands.(bname).razd        = wrapTo360(cfg.reflector_az0_deg + o.bands.(bname).time_from_start_msec/1000/60*360*cfg.reflector_rpm);

% Convert reflector azimuth to L1B standard scan angle (0-180=for; 180-360=aft).
o.bands.(bname).sang        = wrapTo360(90 - o.bands.(bname).razd);

o.bands.(bname).is   = int32(zeros(dims(2:3)));
gind_aft = find(o.bands.(bname).sang > 180);
o.bands.(bname).is(gind_aft) = 1;

% Create a prototype array to be filled in by the integrator.
o.bands.(bname).lat  = nan(dims);

%o.bands.(bname).lon  = nan(dims);
%o.bands.(bname).eia  = nan(dims);
%o.bands.(bname).eaa  = nan(dims);
%o.bands.(bname).pra  = nan(dims);

%o.bands.(bname).ta_x = nan(dims);
%o.bands.(bname).ta_y = nan(dims);
%o.bands.(bname).ta_3 = nan(dims);
%o.bands.(bname).ta_4 = nan(dims);

fprintf('%2d %2d %2s %10s\n', ib, ih, bname, hname);

end

end

eval(['save -V6 ' ofn ' -struct o']);

return



%========================================================================================


% Transform (xi,eta) coordinates to (u,v,z) coordinates for a given reflector tilt angle (off-nadir angle).
%
%  eta  =  u
%  xi   = -v
%  zeta =  z
%
function [u, v, z] = xe_to_uv(xi, eta, zeta, tilt_angle_deg)

if (zeta < -10)
  zeta = -sqrt(1 - xi.^2 - eta.^2);
end
  
if (zeta > 10)
  zeta = sqrt(1 - xi.^2 - eta.^2);
end
  
alt_uv = 90 + tilt_angle_deg;
alt_xe = 90;

[u, v, z] = rot_uv2(eta, -xi, zeta, 0.0, alt_xe, 0.0, alt_uv);

u =  real(u);
v =  real(v);
z =  real(z);

return;


%========================================================================================


% Transform (u,v) coordinates to (xi,eta,zeta) coordinates for a given reflector tilt angle (off-nadir angle).
%
%  eta  =  u
%  xi   = -v
%  zeta =  z
%
function [xi, eta, zeta] = uv_to_xe(u, v, z, tilt_angle_deg)

if (z < -10)
  z = -sqrt(1 - u.^2 - v.^2);
end
  
if (z > 10)
  z = sqrt(1 - u.^2 - v.^2);
end
  
alt_uv = 90 + tilt_angle_deg;
alt_xe = 90;

[eta, xi, zeta] = rot_uv2(u, v, z, 0.0, alt_uv, 0.0, alt_xe);

xi   = -real(xi);
eta  =  real(eta);
zeta =  real(zeta);

return;




%========================================================================================

% Rotate b (director cosine coordinates) from (az0,alt0) to (az1,alt1). Output in director cosine coordinates.
function [ur, vr, zr] = rot_uv2(u, v, z, az0, alt0, az1, alt1)

% Steps:
%  Rotate to 0 azimuth.
%  Rotate about the y-axis to desired altitude.
%  Rotate about z-axis from 0 azimuth to the desired azimuth.

dco = [u(:)'; v(:)'; z(:)'];

dcr = rot_uv(dco, az0, alt0, az1, alt1);

ur = reshape(dcr(1,:), size(u));
vr = reshape(dcr(2,:), size(u));
zr = reshape(dcr(3,:), size(u));

return;


%========================================================================================


% Rotate b (director cosine coordinates) from (az0,alt0) to (az1,alt1). Output in director cosine coordinates.
function [br] = rot_uv(b, az0, alt0, az1, alt1)

% Steps:
%  Rotate to 0 azimuth.
%  Rotate about the y-axis to desired altitude.
%  Rotate about z-axis from 0 azimuth to the desired azimuth.

altrot = 0.0;
azrot  = 0.0;

% Compute cartesian coordinates of the fixed alt-azimuth grid.
%b(1,:) = cos(deg2rad(galt(:))).*cos(deg2rad(gaz(:)));
%b(2,:) = cos(deg2rad(galt(:))).*sin(deg2rad(gaz(:)));
%b(3,:) = sin(deg2rad(galt(:)));

% Rotate about z-axis to 0 azimuth.
rot0(1,1) =  cos(deg2rad(-az0)); 
rot0(1,2) = -sin(deg2rad(-az0)); 
rot0(1,3) =  0;
rot0(2,1) =  sin(deg2rad(-az0));
rot0(2,2) =  cos(deg2rad(-az0));
rot0(2,3) =  0;
rot0(3,1) =  0;
rot0(3,2) =  0;
rot0(3,3) =  1;

% Rotate about the y-axis to the new altitude alt1.
dalt = alt1 - alt0;
rot1(1,1) =  cos(deg2rad(dalt));
rot1(1,2) =  0;
rot1(1,3) = -sin(deg2rad(dalt));
rot1(2,1) =  0;
rot1(2,2) =  1;
rot1(2,3) =  0;
rot1(3,1) =  sin(deg2rad(dalt));
rot1(3,2) =  0;
rot1(3,3) =  cos(deg2rad(dalt));

% Rotate about the z-axis to the new azimuth az1.
daz = az1;
rot2(1,1) =  cos(deg2rad(daz)); 
rot2(1,2) = -sin(deg2rad(daz)); 
rot2(1,3) =  0;
rot2(2,1) =  sin(deg2rad(daz));
rot2(2,2) =  cos(deg2rad(daz));
rot2(2,3) =  0;
rot2(3,1) =  0;
rot2(3,2) =  0;
rot2(3,3) =  1;

% Create rotated cartesian vectors.
br = rot2*rot1*rot0*b;

% Convert back to spherical coordinates.
%raz  = wrapTo360(reshape(rad2deg(atan2(br(2,:),br(1,:))),size(gaz)));
%ralt = reshape(rad2deg(asin(br(3,:))),size(galt));

return
