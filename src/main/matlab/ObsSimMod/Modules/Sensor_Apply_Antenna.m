% dummy script to test things

function Sensor_Apply_Antenna( configurationParameters, inputs, outputs)


global E2E_HOME


idfunction = 'Sensor_Apply_Antenna';

% making global the name of the module to parse
%	 parse the output files names in other
%	 modules
global DUMMY

global SENSOR_SIMULATION




%=== Handling module inputs


%= get log class saved as global variable

global LOG


%= initialize command line parsing class

clp = CLP (configurationParameters, inputs, outputs);


%= Get inputs, outputs and configuration files using

conf1 = clp.getConfFile(1);
conf2 = clp.getConfFile(2);

dirin  = clp.getInputFile(1);
dirout = clp.getOutputFile(1);

LOG.info([ idfunction, ' ** Input folder: ', dirin ])
LOG.info([ idfunction, ' ** Output folder: ', dirout ])



%= Creating folder input and outputs if not existing already


if ~exist( dirin, 'dir' )
  LOG.info( [ idfunction, ' ** Creating folder ', dirin ]);
  mkdir( dirin );
end  



if ~exist( dirout, 'dir' )
  LOG.info( [ idfunction, ' ** Creating folder ', dirout ]);
  mkdir( dirout );
end  

SENSOR_SIMULATION = [ dirout, '/', idfunction ];



%= Parse configuration files 

cfm1 = ConFM(conf1);
cfm2 = ConFM(conf2);



%= Read parameters

LOG.info( [ idfunction, ' ** Reading parameters from global configuration file']);

geodata_version = cfm1.getParameter('geodata_version').getValue;
LOG.info( [ idfunction, ' ** geodata_version has value ', geodata_version ]);

software_version = cfm1.getParameter('software_version').getValue;
LOG.info( [ idfunction, ' ** software_version has value ', software_version ]);

DUMMY = cfm2.getParameter('dummy').getValue;
LOG.info( [ idfunction, ' ** Reading ', DUMMY ]);


%=== Some checks

if ~strcmp( geodata_version,'v1.3') 
  osfi_error('The geodata version given cannot be treated by the software');
end

if ~strcmp( software_version,'v1.3') 
  osfi_error('The software version given do not correspond to this software package');
end


%=== dummy action

novalue = -999;

filesave = [ dirout, '/', idfunction, '_Output.asc' ];  
LOG.info( [ idfunction, ' ** Saving ascii file ', filesave ] );
save(filesave,'novalue','-ascii','-double');



%=== to finish coding

if 0

% Set the center date of simulation.
date_center = cdate_asc;

% length of simulation in hours.
time_length_hours = 0.32;

% Create the scene Tb input file.
ifn = 'cimr_sceps_toa_card_devalgo_polarscene_1_20161217_v2p0_aa_000.nc';
ofn = 'cimr_sceps_toa_card_devalgo_polarscene_1_20161217_v2p0_aa_000.card.mat';
create_cimr_test_scene_polarscene_1_20161217_v2p0_aa_000(ifn, ofn);

% Create the skeletion files for small and big time step runs.
create_cimr_sample_times_bigstep(date_center,   time_length_hours, 'cimr_sceps_toa_card_devalgo_polarscene_1_20161217_v2p0_aa_000.bigstep.skeleton.mat');
create_cimr_sample_times_smallstep(date_center, time_length_hours, 'cimr_sceps_toa_card_devalgo_polarscene_1_20161217_v2p0_aa_000.smallstep.skeleton.mat');

% Now run one of:
%
% ./run_cimr_compute_l1b.csh cimr_l1b_carlos4_bigstep_tot.cfg
% ./run_cimr_compute_l1b.csh cimr_l1b_carlos4_smallstep_tot.cfg
%
% Then (re)convert to netCDF:
%
% ./run_cimr_l1b_mat_to_netcdf_minvars_nom_nedt_apc_tot.csh cimr_l1b_devalgo_polarscene_1.bigstep.mat     cimr_l1b_devalgo_polarscene_1.bigstep.nc
% ./run_cimr_l1b_mat_to_netcdf_minvars_nom_nedt_apc_tot.csh cimr_l1b_devalgo_polarscene_1.smallstep.mat   cimr_l1b_devalgo_polarscene_1.smallstep.nc
%



%= lat lon  and data thinning with lat lon filters

datafile = [ idatafile, '_atmosphere', '.nc' ]; 


lat        = single(ncread( datafile, 'latitude'));
nla        = size(lat,2);
lon        = single(ncread( datafile, 'longitude'));
nlo        = size(lon,1);
ioa        = find( lat >= min(latitude_filter) & lat < max(latitude_filter) & lon >= min(longitude_filter) & lon < max(longitude_filter));

end

return

