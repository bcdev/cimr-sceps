%-------------------------------------------------------------------------------
%
% MODULE  GeoInputs_Extract
%
%    This module extracts from the CIMR netcdf database the geophysical
%    fields needed for a global or polar scene forward model simulation
%    and store them in the adequate format to interact with the forward 
%    models.
%
% FORMAT   GeoInputs_Extract( configurationParameters, inputs, outputs)
%        
% IN    configurationParameters   char		Typically 'Global_Configuration.xml,
%  				  		XXX_Local_Configuration.xml' with XXX
%						the name of module. Used to pass
%						the name of the global and local
%						configuration xml files. Respect naming
%						to be compliant with ESA E2E ICD.
%					
% 	inputs			  char		Name of folder wit input data files,
%				  		typically 'XXX_input' to be compliant
%						with ESA E2E ICD.
%
%	outputs			  char		Name of folder wit output data files,
%				  		typically 'XXX_output' to be compliant
%						with ESA E2E ICD.
%
%-------------------------------------------------------------------------------
% Project:	  CIMR Algorithm Performance Evaluation
% Package:	  CIMR Scientific Work Bench
% Developer:	  Estellus 
% Contact:	  carlos.jimenez@estellus.fr 
% Initiated:	  2019-01-08
%-------------------------------------------------------------------------------

function GeoInputs_Extract( configurationParameters, inputs, outputs)

global E2E_HOME


idfunction = 'GeoInputs_Extract';

% making global the name of the module to parse
%	 parse the output files names in other
%	 modules
global GEOINPUT_SIMULATION

% Making scene_date as a global variable for this simulation
global SCENE_DATE

% Making scene_type as global variable
global SCENE_TYPE


%=== Handling module inputs


%= get log class saved as global variable

global LOG


%= initialize command line parsing class

clp = CLP (configurationParameters, inputs, outputs);


%= Get inputs, outputs and configuration files using

conf1 = clp.getConfFile(1);
conf2 = clp.getConfFile(2);

dirout = clp.getOutputFile(1);

LOG.info([ idfunction, ' ** Output folder: ', dirout ])


%= Creating folder for outputs if not existing already


if ~exist( dirout, 'dir' )
  LOG.info( [ idfunction, ' ** Creating folder ', dirout ]);
  mkdir( dirout );
end  

GEOINPUT_SIMULATION = [ dirout, '/', idfunction ];


%= Parse configuration files 

cfm1 = ConFM(conf1);
cfm2 = ConFM(conf2);


%= Read parameters

LOG.info( [ idfunction, ' ** Reading parameters from global configuration file']);

geodata_version = cfm1.getParameter('geodata_version').getValue;
LOG.info( [ idfunction, ' ** geodata_version has value ', geodata_version ]);

software_version = cfm1.getParameter('software_version').getValue;
LOG.info( [ idfunction, ' ** software_version has value ', software_version ]);

SCENE_TYPE = cfm2.getParameter('scene_type').getValue;
LOG.info( [ idfunction, ' ** Extracting fields for ', SCENE_TYPE, ' scene' ]);

SCENE_DATE = cfm2.getParameter('scene_date').getValue;
LOG.info( [ idfunction, ' ** Extracting fields for day ', SCENE_DATE ]);

%=data_folder = cfm2.getParameter('data_folder').getValue;
%=data_folder = strrep( data_folder, 'E2E_HOME', E2E_HOME );
%=data_file_name = cfm2.getParameter('data_file_name').getValue;
%=datafile = [ data_folder, '/', data_file_name ];

datafile = clp.getInputFile(1);
LOG.info( [ idfunction, ' ** Extracting fields stored at ', datafile ]);
idatafile = [ datafile ];

latitude_filter = cfm2.getParameter('latitude_filter').getValue;
LOG.info( [ idfunction, ' ** Extracting fields for latitude ', num2str(latitude_filter)]);

longitude_filter = cfm2.getParameter('longitude_filter').getValue;
LOG.info( [ idfunction, ' ** Extracting fields for longitude ', num2str(longitude_filter)]);



%=== Some checks

if ~strcmp( geodata_version,'v1.3') 
  osfi_error('The geodata version given cannot be treated by the software');
end

if ~strcmp( software_version,'v1.3') 
  osfi_error('The software version given do not correspond to this software package');
end



%=== Atmosphere ===============================================================
%    profiles			    


%= lat lon  and data thinning with lat lon filters

%=datafile = [ idatafile, '_atmosphere', '.nc' ]; 
datafile = clp.getInputFile(1);
LOG.info( [ idfunction, ' ** Read atmosphere file']);


lat        = single(ncread( datafile, 'latitude'));
nla        = size(lat,2);
lon        = single(ncread( datafile, 'longitude'));
nlo        = size(lon,1);
ioa        = find( lat >= min(latitude_filter) & lat < max(latitude_filter) & lon >= min(longitude_filter) & lon < max(longitude_filter));


clon       = lon(:,1);
io         = find(  clon >= min(longitude_filter) & clon < max(longitude_filter));
lio        = length(io);

clat       = lat(1,:);
ia        = find( clat >= min(latitude_filter) & clat < max(latitude_filter));
lia        = length(ia);

if isempty( io ) | isempty( ia )
  osfi_error('There is no data left after applying the latitude and longitude filters');
end



lat        = single(ncread( datafile, 'latitude',[io(1) ia(1) 1],[lio lia 1]));
lon        = single(ncread( datafile, 'longitude',[io(1) ia(1) 1],[lio lia 1]));
ioa        = 1:(size(lat,1)*size(lat,2));
nlo        = size(lon,1);
nla        = size(lat,2);
clear li lia clat clon


% lons are in 1E-180E-179W-1W, changing to 1-360

lon( lon < 0 ) = lon( lon < 0 ) + 360;


% original lats and lons
oLAT = lat;
oLON = lon;



%= P
data    = single(ncread( datafile, 'pressure_level'));
% canging form ECMWF convention to make from ground to the top
ind_pre = [ length(data):-1:1];
data    = data(ind_pre);  




sparam  = 'P';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' field for day ', SCENE_DATE ]);
save( outfile, 'data' );
P       = data;
npr     = length(P);

if P(1) < P(2)
  osfi_error('Atmospheric pressures are expected as a vector from the ground to the top ');
end

%= T

aux    = single(zeros(nlo,nla,npr));
for p = 1:npr
    ip         = ind_pre(p);
    adata      = single((ncread( datafile, 'temperature_profile',[io(1) ia(1) ip 1],[nlo nla 1 1])));
    aux(:,:,p) = adata;
end



data    = single(zeros(npr,nlo*nla));
for p=1:npr
 auc       = squeeze(aux(:,:,p));
 data(p,:) = auc(:);
end
data = data';
sparam  = 'T';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
data = single(data);
save( outfile, 'data','-v7.3' );
T       = data;
clear aux auc data


%= LW
%  Kg/Kg to g/m3

aux    = single(zeros(nlo,nla,npr));
for p = 1:npr
  ip         = ind_pre(p);
  adata      = single((ncread( datafile, 'tclw_profile',[io(1) ia(1) ip 1],[nlo nla 1 1])));
  aux(:,:,p) = adata;
end


data    = single(zeros(npr,nlo*nla));
for p=1:npr
 auc        = squeeze(aux(:,:,p));
 data(p,:) = auc(:);
end
data = data';


iP   = 1e2 * repmat( P', length(ioa), 1); % to Pascal
data = 1e3 * mixing_ratio_to_density(data, T, iP); % to g/m3

sparam  = 'LW';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
data = single(data);
save( outfile, 'data','-v7.3' );
clear aux auc data



%= IW g/m3

data     = single(zeros(npr,nlo*nla));
data = data(:,ioa);
data = data';
sparam  = 'IW';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
data = single(data);
save( outfile, 'data','-v7.3' );
clear data




%= RHO
%  Kg/Kg to g/m3


aux    = single(zeros(nlo,nla,npr));
for p = 1:npr
    ip         = ind_pre(p);
    adata      = single((ncread( datafile, 'tcwv_profile',[io(1) ia(1) ip 1],[nlo nla 1 1])));
    aux(:,:,p) = adata;
end


data    = single(zeros(npr,nlo*nla));
for p=1:npr
 auc        = squeeze(aux(:,:,p));
 data(p,:) = auc(:);
end
data = data';
iP   = 1e2 * repmat( P', length(ioa), 1); % to Pascal
data = 1e3 * mixing_ratio_to_density(data, T, iP); % to g/m3

sparam  = 'RHO';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
data = single(data);
save( outfile, 'data','-v7.3' );
clear aux auc data



%=== Atmosphere ===============================================================
%    total column			    


%=datafile = [ idatafile, '_surface', '.nc' ]; 
datafile = clp.getInputFile(3);


%= RHO total column
%  Kg/m2 


data    = single(ncread( datafile, 'total_column_water_vapour',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:);
data    = data(ioa);
sparam  = 'CWVC';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
data = single(data);
save( outfile, 'data' );
clear data





%= LW total column
%  Kg/m2 


data    = single(ncread( datafile, 'total_column_liquid_water',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:);
data    = data(ioa);
sparam  = 'CLWC';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
data = single(data);
save( outfile, 'data' );
clear data





%=== Extraction surface ===============================================================
%    Land-sea-ice-coast masks


%=datafile = [ idatafile, '_surface', '.nc' ]; 
datafile = clp.getInputFile(3);

%= SIC 

iSIC    = single(ncread( datafile, 'sea_ice_concentration',[io(1) ia(1) 1],[nlo nla 1]));
iSIC    = iSIC/100;
iSIC    = iSIC(:);


%= mask conventions
%  1=open-sea, 2=land, 5=open-lake, 
%  9=open-sea with ice in the grid, 
%  13=open-lake with ice in the grid"

lsi    = single(ncread( datafile, 'land_sea_ice_mask',[io(1) ia(1) 1],[nlo nla 1]));
ind_ice = find(lsi == 9 | lsi == 13);
ind_sea = find(lsi == 1);
ind_lan = find(lsi == 2 | lsi == 5);


mask          = zeros(nlo*nla,1);
mask(ind_sea) = 1;
mask(ind_ice) = 2;
mask(ind_lan) = 3;
mask = single(mask);

% saving mask


data    = mask(:);
data    = data(ioa);
sparam  = 'MASK';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );
clear data aux_lan aux_sea aux_ice lat lon lsi %mask




%= indexes to selected lon-lat

aux     = zeros(nlo,nla);
aux(ind_sea) = 1;
ind_sea = aux(ioa);
ind_sea = find( ind_sea == 1);

aux     = zeros(nlo,nla);
aux(ind_lan) = 1;
ind_lan = aux(ioa);
ind_lan = find( ind_lan == 1);

aux     = zeros(nlo,nla);
aux(ind_ice) = 1;
ind_ice = aux(ioa);
ind_ice = find( ind_ice == 1);

aSIC    = iSIC;
iSIC    = iSIC(ioa);



%=== MIXED SCENES
% we need sea fields for mixed scenes of
% sea and ice, i.e., when 0 < SIC < 100 and
% coastal pixels

ind_sea = union( ind_sea, ind_ice );

% coastal pixels need to be as land as well

non_sea    = setdiff( 1:length(ioa), ind_sea);
non_lan    = setdiff( 1:length(ioa), ind_lan);
non_ice    = setdiff( 1:length(ioa), ind_ice);



%=== time stamp
%    Using amsr_scan time covreted to matlab datenum
%    i.e. number of days from  01-01-0000 midnight

%=datafile = [ idatafile, '_observations', '.nc' ]; 
datafile = clp.getInputFile(2);

data    = single(ncread( datafile, 'amsr2_scan_time',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:);
data    = data(ioa);
data    = data/60/60/24 + datenum('1993-01-01'); 

sparam  = 'UTC';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );



%=== Sea ===============================================================
%			    SST	[K]		Surface Temperature
%			    SSS [g/kg]		Sea Surface Salinity
%			    OWS	[m/s]		Wind velocity
%			    SLP [mbar]		Mean Sea Level Pressure


%=datafile = [ idatafile, '_surface', '.nc' ]; 
datafile = clp.getInputFile(3);

data    = single(ncread( datafile, 'sea_surface_salinity',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:);
data    = data(ioa);

data(non_sea) = -999;    

sparam  = 'SSS';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );
SSS = data;


%= SST K

data    = single(ncread( datafile, 'sea_surface_temperature',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:);
data    = data(ioa);



sparam  = 'SST';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );




%= OWS m/s


auxu    = single(ncread( datafile, 'surface_wind_u',[io(1) ia(1) 1],[nlo nla 1]));
auxv    = single(ncread( datafile, 'surface_wind_v',[io(1) ia(1) 1],[nlo nla 1]));

data    = sqrt( auxu.^2 + auxv.^2 );
data    = data(:);
data    = data(ioa);
data(non_sea) = -999;    
sparam  = 'OWS';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );


data    = auxu;
data    = data(:);
data    = data(ioa);
data(non_sea) = -999;    
sparam  = 'UWS';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );


data    = auxv;
data    = data(:);
data    = data(ioa);
data(non_sea) = -999;    
sparam  = 'VWS';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );



%= SLP Pa to hPa=mbar

data    = single(ncread( datafile, 'mean_sea_level_pressure',[io(1) ia(1) 1],[nlo nla 1]));
data    = data/1e2;
SLP     = data;
SLP     = SLP(:);
SLP     = SLP(ioa);
data    = SLP;
data(non_sea) = -999;    
sparam  = 'SLP';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );




%=== Ice
%			    IST	[K]		Ice Surface Temperature 
%			    T2M [K]		2-meter Air Temperature
%			    SIC [-]		Sea Ice Concentration 0-1
%			    SIT [m]		Sea ice Thickness
%			    SND [m]		Snow Depth 
%			    LAT [degrees]	Latitude, 90S to 90N
%			    ILP [mbar]          Land Level Pressure
%

%=datafile = [ idatafile, '_surface', '.nc' ]; 
datafile = clp.getInputFile(3);

%= LAT

data    = oLAT(:); 
LAT     = data;
sparam  = 'LAT';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );

data    = oLAT; 
sparam  = 'nonvectorizedLAT';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );


%= IST 


data    = single(ncread( datafile, 'sea_ice_temperature',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:); 
data    = data(ioa);
data(non_ice) = -999;

sparam  = 'IST';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );


%= T2M

data    = single(ncread( datafile, 'temperature_2m',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:);

data    = data(ioa);
sparam  = 'T2M';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );


%= SIC

data          = iSIC;
data(non_ice) = -999;    
sparam  = 'SIC';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );
clear iSIC


%= SNOW DEPTH

data    = single(ncread( datafile, 'snow_thickness',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:);
data    = data(ioa);
data(non_ice) = -999;

sparam  = 'SND';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );



%= FREEZING DAYS NUMBER

data    = single(ncread( datafile, 'freezing_days_number',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:);
data    = data(ioa);
data(non_ice) = -999;    
sparam  = 'FDN';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );




%= SIT m

data    = single(ncread( datafile, 'sea_ice_thickness',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:);
data    = data(ioa);

data(non_ice) = -999;    
sparam  = 'SIT';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );



%= Ice type 0 FYI / 1 MYI

data    = single(ncread( datafile, 'sea_ice_type',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:);
data    = data(ioa);

data(non_ice) = -999;    
sparam  = 'XYI';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );



%= ILP using slp
data    = SLP;
data(non_ice) = -999;    



sparam  = 'ILP';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );



%=== Land
%			    LST	[K]		Land Surface Temperature
%			    LAT [degrees]	Latitude, 90S to 90N
%			    LON [degrees]	Longitude, 0 to 360
%			    MONTH []		Month of the year
%			    LLP [mbar]          Land Level Pressure
%			    LSP [0-1]		Land Surface Percentage,
%						optional, required for a
%						mixed land and sea scene
%			    LIP [0-1]		Land Surface Percentage,
%						optional, required for a
%						mixed land and ice scene



%= LST  
data    = single(ncread( datafile, 'skin_temperature',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:);
data    = data(ioa);

data(non_lan) = -999;    
sparam  = 'LST';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );




%= MONTH
imonth  = SCENE_DATE(5:6);
data    = str2num(imonth)*ones(size(data));
sparam  = 'MONTH';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );



%= LON  to 0-360
data    = oLON(:);
sparam  = 'LON';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );

data    = oLON;
sparam  = 'nonvectorizedLON';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );



%= LLP 
data    = SLP;
data(non_lan) = -999;    
sparam  = 'LLP';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );
clear SLP



%= hard coded switch to output the observations
%  for test-debugging purposes

do_obs = 1;

if do_obs

datafile = [ idatafile, '_observations', '.nc' ]; 
datafile = clp.getInputFile(2);


%= TB 1.4V

data    = single(ncread( datafile, 'smap_tb_v',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:);
data    = data(ioa); 
sparam  = 'SMAP_1V';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );


%= TB 1.4H

data    = single(ncread( datafile, 'smap_tb_h',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:);
data    = data(ioa); 
sparam  = 'SMAP_1H';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );

%= TB 6V

data    = single(ncread( datafile, 'amsr2_tb_res06_6v',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:);
data    = data(ioa); 
sparam  = 'AMSR2_6V';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );



%= TB 6H

data    = single(ncread( datafile, 'amsr2_tb_res06_6h',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:);
data    = data(ioa); 
sparam  = 'AMSR2_6H';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );



%= TB 10V

data    = single(ncread( datafile, 'amsr2_tb_res10_10v',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:);
data    = data(ioa); 
sparam  = 'AMSR2_10V';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );


%= TB 10H

data    = single(ncread( datafile, 'amsr2_tb_res10_10h',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:);
data    = data(ioa); 
sparam  = 'AMSR2_10H';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );


%= TB 18V

data    = single(ncread( datafile, 'amsr2_tb_res23_18v',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:);
data    = data(ioa); 
sparam  = 'AMSR2_18V';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );



%= TB 18H

data    = single(ncread( datafile, 'amsr2_tb_res23_18h',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:);
data    = data(ioa); 
sparam  = 'AMSR2_18H';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );


%= TB 23V

data    = single(ncread( datafile, 'amsr2_tb_res23_23v',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:);
data    = data(ioa); 
sparam  = 'AMSR2_23V';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );


%= TB 23H

data    = single(ncread( datafile, 'amsr2_tb_res23_23h',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:);
data    = data(ioa); 
sparam  = 'AMSR2_23H';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );



%= TB 36V

data    = single(ncread( datafile, 'amsr2_tb_res36_36v',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:);
data    = data(ioa); 
sparam  = 'AMSR2_36V';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );




%= TB 36H

data    = single(ncread( datafile, 'amsr2_tb_res36_36h',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:);
data    = data(ioa); 
sparam  = 'AMSR2_36H';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );



%= TB 89V

data    = single(ncread( datafile, 'amsr2_tb_orig_89av',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:);
data    = data(ioa); 
sparam  = 'AMSR2_89V';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );


%= TB 89H

data    = single(ncread( datafile, 'amsr2_tb_orig_89ah',[io(1) ia(1) 1],[nlo nla 1]));
data    = data(:);
data    = data(ioa); 
sparam  = 'AMSR2_89H';
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save( outfile, 'data' );

end

return
