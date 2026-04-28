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
% Project:	  CIMR SCEPS 
% Package:	  SCEPScodes 
% Developer:	  Estellus 
% Contact:	  carlos.jimenez@estellus.fr 
% Initiated:	  2019-01-08
% Updated:	  2025-04-29
%-------------------------------------------------------------------------------


function GeoInputs_Extract( configurationParameters, inputs, outputs)

warning off


idfunction = 'GeoInputs_Extract';


%=== Handling module inputs


%= get log class saved as global variable
%  global LOG - original, but for runtime 
%  better to call the Logger

LOG = Logger();


%= initialize command line parsing class

clp = CLP (configurationParameters, inputs, outputs);


%= Get inputs, outputs and configuration files using

conf1 = clp.getConfFile(1);
conf2 = clp.getConfFile(2);

dirout = clp.getOutputFile(1);

LOG.info([ idfunction, ' ** Output folder: ', dirout ])


%= Parse configuration files 

cfm1 = ConFM(conf1);
cfm2 = ConFM(conf2);


BYPASS_READING = cfm2.getParameter('bypass_reading').getValue;


if BYPASS_READING == 1  

  sparam  = 'SNTH';
  outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
  if exist( outfile, 'file' ) 
    LOG.info( [ idfunction, ' **  Bypassing the reading of fields, there exist already' ]);
  else
    LOG.info( [ idfunction, ' **  Bypassing the reading of fields demanded, but no fields, so they are preapred' ]);
    BYPASS_READING == 0  
  end

end

if ~BYPASS_READING


%= Creating folder for outputs if not existing already

if exist( dirout, 'dir' )

  % removing folder in case left from previous run
   eval(['!rm -r ', dirout ])
  
end

% creating folder
LOG.info( [ idfunction, ' ** Creating folder ', dirout ]);
mkdir( dirout );
  



%= Read parameters

LOG.info( [ idfunction, ' ** Reading parameters from global configuration file']);

geodata_version = cfm1.getParameter('geodata_version').getValue;
LOG.info( [ idfunction, ' ** geodata_version has value ', geodata_version ]);

software_version = cfm1.getParameter('software_version').getValue;
LOG.info( [ idfunction, ' ** software_version has value ', software_version ]);

SCENE_TYPE = cfm2.getParameter('scene_type').getValue;
LOG.info( [ idfunction, ' ** Extracting fields for ', SCENE_TYPE, ' scene' ]);



%= with more flexibility in case observation file do not exsit


datafile = clp.getInputFile(1);
if ~isempty(strfind(datafile, '_atmosphere') )
  datafile_atmosphere = datafile;
elseif ~isempty( strfind(datafile, '_surface'))
  datafile_surface = datafile;
end


datafile = clp.getInputFile(2);
if ~isempty(strfind(datafile, '_atmosphere') )
  datafile_atmosphere = datafile;
elseif ~isempty( strfind(datafile, '_surface'))
  datafile_surface = datafile;
end


if ~exist( 'datafile_atmosphere', 'var' ) | ~exist( 'datafile_surface', 'var' )
  osfi_error('There are missing input data files');
else
  LOG.info( [ idfunction, ' ** Extracting fields stored at ', datafile_atmosphere ]);
  LOG.info( [ idfunction, ' ** Extracting fields stored at ', datafile_surface ]);
end


if ~exist( 'datafile_obs', 'var' )
  do_obs = 0;
else
  do_obs = 1;
  LOG.info( [ idfunction, ' ** Extracting fields stored at ', datafile_obs ]);
end



latitude_filter = cfm2.getParameter('latitude_filter').getValue;
LOG.info( [ idfunction, ' ** Extracting fields for latitude ', num2str(latitude_filter)]);

longitude_filter = cfm2.getParameter('longitude_filter').getValue;
LOG.info( [ idfunction, ' ** Extracting fields for longitude ', num2str(longitude_filter)]);



%=== Some checks

if ~strcmp( geodata_version,'v2.1') 
  osfi_error('The geodata version given cannot be treated by the software');
end

if ~strcmp( software_version,'v2.1') 
  osfi_error('The software version given do not correspond to this software package');
end





%= lat lon  and data thinning with lat lon filters
%  if grids are regular

datafile = datafile_atmosphere;

if 0

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

  clon       = lon(1,:);
  io         = find(  clon >= min(longitude_filter) & clon < max(longitude_filter));
  lio        = length(io);

  clat       = lat(:,1);
  ia        = find( clat >= min(latitude_filter) & clat < max(latitude_filter));
  lia        = length(ia);
  
end 


if isempty( io ) | isempty( ia )
  osfi_error('There is no data left after applying the latitude and longitude filters');
end

else

  LOG.info( [ idfunction, ' ** Lat-lon filtering not activated']);  

  lat        = single(ncread( datafile, 'latitude'));
  nla        = size(lat,2);
  lon        = single(ncread( datafile, 'longitude'));
  nlo        = size(lon,1);
  io         = 1:nlo;
  lio        = nlo;
  ia         = 1:nla;
  lia        = nla;
  ioa        = 1:(nlo*nla);

end


% lons are in 1E-180E-179W-1W, changing to 1-360 
% lon( lon < 0 ) = lon( lon < 0 ) + 360;


% original lats and lons
oLAT = lat;
oLON = lon;


if min(oLON(:)) < min(longitude_filter) | max(oLON(:)) > max(longitude_filter) |min(oLAT(:)) < min(latitude_filter) | max(oLAT(:)) > max(latitude_filter)
  LOG.info( [ idfunction, ' ** Lat-lon filtering not exact as the lat lon are not regular']);
end

%= saving variable with datatfile name


data = datafile;
aux  = strfind( data, '/');
data = data( aux(end)+1:end );
aux  = strfind( data, '_');
data = data( 1:aux(end)-1);


sparam  = 'DATAFILENAME';
outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' field ' ]);
save( outfile, 'data' );


%=== Ancillary ===============================================================
%    


%= LAT

data    = oLAT(:); 
LAT     = data;
sparam  = 'LAT';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data

data    = oLAT; 
sparam  = 'nonvectorizedLAT';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data


%= LON -180 180
data    = oLON(:);
sparam  = 'LON';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data


data    = oLON;
sparam  = 'nonvectorizedLON';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data
		    


%= time stamp
%  Days from 2000/1/1


data    = single(ncread( datafile, 'time',[io(1) ia(1)],[nlo nla]));
data    = data(:);
data    = data(ioa);

sparam  = 'TIME';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );


%= local time
%  derived from local time

cf = 60*60*24;
[kyear, kmonth, kday, khour]  = unixsecs2date( cf * data + date2unixsecs(2000,1,1) );
ltime = khour + (12 * oLON(:) /180);
kyear  = unique( kyear );
kmonth = unique( kmonth );
kday   = unique( kday );
clear data




%= LOCAL TIME

data    = ltime(ioa);
sparam  = 'LTIME';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data


%= YEAR
data    = kyear;
sparam  = 'YEAR';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data


%= MONTH
data    = kmonth;
sparam  = 'MONTH';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data


%= DAY
data    = kday;
sparam  = 'DAY';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data


clear ltime kyear kmonth kday khour





%=== Atmosphere ===============================================================
%    profiles	
%			   T	[K]		Temperature
%			   P	[Pa]		Pressure
%			   SPRES [Pa]		Surface Pressure
%			   Q    [Kg/Kg]		Specific Humidity
%			   CLWC [Kg/Kg]		Specific Cloud Liquid Water Content 
%			   CRWC [Kg/Kg]		Specific Cloud Rain Water Content 		    




datafile = datafile_atmosphere;


%= P from mbar/hPa to Pa
data    = single(ncread( datafile, 'pressure_level'));
% changing form ECMWF convention to make from ground to the top
ind_pre = [ length(data):-1:1];
data    = data(ind_pre); 

sparam  = 'P';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' field ' ]);
save( outfile, 'data' );
P       = data;
npr     = length(P);
clear data


if P(1) < P(2)
  osfi_error('Atmospheric pressures are expected as a vector from the ground to the top ');
end


%= T

aux    = single(zeros(nlo,nla,npr));
for p = 1:npr
    ip         = ind_pre(p);
    adata      = single((ncread( datafile, 'temperature_profile',[io(1) ia(1) ip],[nlo nla 1])));
    aux(:,:,p) = adata;
end


data    = single(zeros(npr,nlo*nla));
for p=1:npr
 auc       = squeeze(aux(:,:,p));
 data(p,:) = auc(:);
end
data = data';
sparam  = 'T';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
data = single(data);
save( outfile, 'data','-v7.3' );
T       = data;
clear aux auc data adata


%= LW
%  Kg/Kg 

aux    = single(zeros(nlo,nla,npr));
for p = 1:npr
  ip         = ind_pre(p);
  adata      = single((ncread( datafile, 'tclw_profile',[io(1) ia(1) ip],[nlo nla 1])));
  aux(:,:,p) = adata;
end


data    = single(zeros(npr,nlo*nla));
for p=1:npr
 auc        = squeeze(aux(:,:,p));
 data(p,:) = auc(:);
end
data = data';


%iP   = 1e2 * repmat( P', length(ioa), 1); % to Pascal
%data = 1e3 * mixing_ratio_to_density(data, T, iP); % to g/m3

sparam  = 'CLWC';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
data = single(data);
save( outfile, 'data','-v7.3' );
clear aux auc data adata


%= RW
%  Kg/Kg

aux    = single(zeros(nlo,nla,npr));
for p = 1:npr
  ip         = ind_pre(p);
  adata      = single((ncread( datafile, 'tcrw_profile',[io(1) ia(1) ip],[nlo nla 1])));
  aux(:,:,p) = adata;
end


data    = single(zeros(npr,nlo*nla));
for p=1:npr
 auc        = squeeze(aux(:,:,p));
 data(p,:) = auc(:);
end
data = data';


%iP   = 1e2 * repmat( P', length(ioa), 1); % to Pascal
%data = 1e3 * mixing_ratio_to_density(data, T, iP); % to g/m3

sparam  = 'CRWC';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
data = single(data);
save( outfile, 'data','-v7.3' );
clear aux auc data adata




%= Q
%  Kg/Kg


aux    = single(zeros(nlo,nla,npr));
for p = 1:npr
    ip         = ind_pre(p);
    adata      = single((ncread( datafile, 'tcwv_profile',[io(1) ia(1) ip],[nlo nla 1])));
    aux(:,:,p) = adata;
end


data    = single(zeros(npr,nlo*nla));
for p=1:npr
 auc        = squeeze(aux(:,:,p));
 data(p,:) = auc(:);
end
data = data';

sparam  = 'Q';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
data = single(data);
save( outfile, 'data','-v7.3' );
clear aux auc data adata





%=== Atmosphere ===============================================================
%    total column
%			   TCWV	[kg/m2]		Total Column Water Vapour
%			   TCLW	[kg/m2]		Total Column Liquid Water
%			   TCRW	[kg/m2]		Total Column Rain Water			    


datafile = datafile_surface;


%= WV total column
%  Kg/m2 


data    = single(ncread( datafile, 'total_column_water_vapour',[io(1) ia(1)],[nlo nla]));
data    = data(:);
data    = data(ioa);
sparam  = 'TCWV';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
data = single(data);
datr = data;
save( outfile, 'data' );
clear data


%= daily WV total column, this is an optional variable, as it was noted that 
%  for some simulations at different times of the same day the daily variability
%  impacted the simulated emis. If exists it will be used, otherwise the existing
%  WV total colum will be copied here so the RTM gets a value. So it's up to the 
%  user, if a daily value is added to the geo file, it is used, if not the RTM 
%  runs on the WV total column variable.
%
%  Kg/m2 

sparam  = 'DTCWV';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);

try

  %= daily available
  data    = single(ncread( datafile, 'daily_total_column_water_vapour',[io(1) ia(1)],[nlo nla]));
  data    = data(:);
  data    = data(ioa);
  data = single(data);

catch

  %= daily si replaced with existing instant value
  data = datr;

end

save( outfile, 'data' );
clear data datr


%= the latest iteration of RTMs do not need these 
%  fields as it uses the atmospheric profiles
%  readjusted with surface pressure, and an update
%  of the column integrations is used. Removing
%  but leaving the code.

if 0


%= LW total column
%  Kg/m2 


data    = single(ncread( datafile, 'total_column_liquid_water',[io(1) ia(1)],[nlo nla]));
data    = data(:);
data    = data(ioa);
sparam  = 'TCLW';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
data = single(data);
save( outfile, 'data' );
clear data



%= RW total column
%  Kg/m2 


data    = single(ncread( datafile, 'total_column_rain_water',[io(1) ia(1)],[nlo nla]));
data    = data(:);
data    = data(ioa);
sparam  = 'TCRW';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
data = single(data);
save( outfile, 'data' );
clear data

end





%=== Extraction surface ===============================================================
%    Land-sea-ice-coast masks


datafile = datafile_surface;


%= SIC 

iSIC    = single(ncread( datafile, 'sea_ice_concentration',[io(1) ia(1)],[nlo nla]));
iSIC    = iSIC(:);


%= mask conventions
%  1=open-sea, 2=land, 3=seaice
%  9=open-sea with ice in the grid, 
%  13=open-lake with ice in the grid"

lsi    = single(ncread( datafile, 'land_sea_ice_mask',[io(1) ia(1)],[nlo nla]));


ind_sea = find(lsi == 1);
ind_lan = find(lsi == 2);
ind_ice = find(lsi == 3 | lsi == 9);


mask          = zeros(nlo*nla,1);
mask(ind_sea) = 1;
mask(ind_lan) = 2;
mask(ind_ice) = 3;
mask = single(mask);



% saving mask


data    = mask(:);
data    = data(ioa);
sparam  = 'MASK';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data aux_lan aux_sea aux_ice lat lon lsi mask




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



%=== Common surface  ===============================================================
%
%			    SP [Pa]          Surface Air Pressure

%= SP Pa 

data    = single(ncread( datafile, 'surface_air_pressure',[io(1) ia(1)],[nlo nla]));
data    = data;
data    = data(:);
data    = data(ioa);
sparam  = 'SP';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data



%=== Sea surface
%			    SST	[K]		Surface Temperature
%			    SSS [g/kg]		Sea Surface Salinity
%			    OWS	[m/s]		Wind velocity
%			    UWS	[m/s]		U wind component
%			    VWS	[m/s]		V wind component



datafile = datafile_surface;


%=== adding to non_sea when we do not have sss or sst

data    = single(ncread( datafile, 'sea_surface_salinity',[io(1) ia(1)],[nlo nla]));
data    = data(:);
data    = data(ioa);
data(non_sea) = -999;    


sparam  = 'SSS';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data


%= SST K

data    = single(ncread( datafile, 'sea_surface_temperature',[io(1) ia(1)],[nlo nla]));
data    = data(:);
data    = data(ioa);

data(non_sea) = -999;    

sparam  = 'SST';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data


%= OWS m/s

auxu    = single(ncread( datafile, 'surface_wind_u',[io(1) ia(1)],[nlo nla]));
auxv    = single(ncread( datafile, 'surface_wind_v',[io(1) ia(1)],[nlo nla]));

data    = sqrt( auxu.^2 + auxv.^2 );
data    = data(:);
data    = data(ioa);
data(non_sea) = -999;    
sparam  = 'OWS';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );


data    = auxu;
data    = data(:);
data    = data(ioa);
data(non_sea) = -999;    
sparam  = 'UWS';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );


data    = auxv;
data    = data(:);
data    = data(ioa);
data(non_sea) = -999;    
sparam  = 'VWS';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data auxv auxu






%=== Sea ice surface
%			    SIT [K]		Sea ice Surface Temperature
%			    SIC [-]		Sea Ice Concentration 0-1
%			    SIC [-]		Sea Ice Concentration 0-1
%			    SITH [m]		Sea Ice Thickness
%			    SIA [years]		Sea Ice Age
%			    SIR [m]		Sea Ice Roughness
%			    SISNTH [m]		Sea Ice Snow Thickness
%			    LAT [degrees]	Latitude, 90S to 90N
%			    LAT [degrees]	Latitude, 180W to 180E
%			    TCLW [kg/m2]		Total Column Liquid Water



%=datafile = [ idatafile, '_surface', '.nc' ]; 
datafile = datafile_surface;


%= SIC

data          = iSIC;


data(non_ice) = -999;    
sparam  = 'SIC';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear iSIC



%= SITH m

data    = single(ncread( datafile, 'sea_ice_thickness',[io(1) ia(1)],[nlo nla]));
data    = data(:);
data    = data(ioa);

data(non_ice) = -999;    
sparam  = 'SITH';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data


%= SIR m

data    = single(ncread( datafile, 'sea_ice_roughness',[io(1) ia(1)],[nlo nla]));
data    = data(:);
data    = data(ioa);

data(non_ice) = -999;    
sparam  = 'SIR';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data


%= SIA m

data    = single(ncread( datafile, 'sea_ice_age',[io(1) ia(1)],[nlo nla]));
data    = data(:);
data    = data(ioa);

data(non_ice) = -999;    
sparam  = 'SIA';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data


%= SIT 
%  m

data    = single(ncread( datafile, 'sea_ice_temperature',[io(1) ia(1)],[nlo nla]));
data    = data(:); 
data    = data(ioa);
data(non_ice) = -999;
datr    = data;

sparam  = 'SIT';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data


%= daily SIT 
%  Same considerations for daily values as for daily WV total column


sparam  = 'DSIT';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);

try

  data    = single(ncread( datafile, 'daily_sea_ice_temperature',[io(1) ia(1)],[nlo nla]));
  data    = data(:); 
  data    = data(ioa);
  data(non_ice) = -999;

catch

  data = datr;

end

save( outfile, 'data' );
clear data datr



%= SISNT m

data    = single(ncread( datafile, 'sea_ice_snow_thickness',[io(1) ia(1)],[nlo nla]));
data    = data(:);
data    = data(ioa);

data(non_ice) = -999;    
sparam  = 'SISNTH';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data




%=== Land climato
%			    LST	[K]		Land Skin Surface Temperature
%			    LAT [degrees]	Latitude, 90S to 90N
%			    LON [degrees]	Longitude, 0 to 360
%			    MONTH []		Month of the year



%= LST  
data    = single(ncread( datafile, 'skin_temperature',[io(1) ia(1)],[nlo nla]));
data    = data(:);
data    = data(ioa);

data(non_lan) = -999;    
sparam  = 'LST';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data




%=== Snowed land
%			    DAY[]		day of month
%			    T2M [K]		2-meter Air Temperature
%			    LTIME [hours]	Decimal local time
%			    LC[]		Lake cover
%			    SDOR[m/km?]		Standard deviation of Orography
%			    SNA[]		Snow Albedo
%			    SND[]		Snow Density
%			    SNTH[m]		Snow thickness
%			    SNTL1[m]		Soil temperature Layer 1
%			    SNT[K]		Snow temperature
%			    SNMLT[m]		Snow melt




%= T2M

data    = single(ncread( datafile, 'temperature_2m',[io(1) ia(1)],[nlo nla]));
data    = data(:);
data    = data(ioa);

sparam  = 'T2M';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data



%= LAKE COVER

data    = single(ncread( datafile, 'lake_cover',[io(1) ia(1)],[nlo nla]));
data    = data(:);
data    = data(ioa);
sparam  = 'LC';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data


%= Std of orography

data    = single(ncread( datafile, 'std_orography',[io(1) ia(1)],[nlo nla]));
data    = data(:);
data    = data(ioa);
sparam  = 'SDOR';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data


%= Snow albedo

data    = single(ncread( datafile, 'snow_albedo',[io(1) ia(1)],[nlo nla]));
data    = data(:);
data    = data(ioa);
data(non_lan) = -999;    
datr    = data;

sparam  = 'SNA';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data


%= Daily Snow albedo
%  Same considerations for daily values as for daily WV total column

sparam  = 'DSNA';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);

try

  data    = single(ncread( datafile, 'daily_snow_albedo',[io(1) ia(1)],[nlo nla]));
  data    = data(:);
  data    = data(ioa);
  data(non_lan) = -999;    

catch
 
  data = datr;

end

save( outfile, 'data' );
clear data datr



%= Snow Density
data    = single(ncread( datafile, 'snow_density',[io(1) ia(1)],[nlo nla]));
data    = data(:);
data    = data(ioa);
data(non_lan) = -999;    

sparam  = 'SND';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data


%= SNOW thickness

data    = single(ncread( datafile, 'snow_thickness',[io(1) ia(1)],[nlo nla]));
data    = data(:);
data    = data(ioa);
data(non_lan) = -999;

sparam  = 'SNTH';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data


%= SNOW temperature

data    = single(ncread( datafile, 'snow_temperature',[io(1) ia(1)],[nlo nla]));
data    = data(:);
data    = data(ioa);
data(non_lan) = -999;

sparam  = 'SNT';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data


%= SNOW temperature 

data    = single(ncread( datafile, 'snow_temperature',[io(1) ia(1)],[nlo nla]));
data    = data(:);
data    = data(ioa);
data(non_lan) = -999;

sparam  = 'SNT';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data


%= Soil layer 1 temperature

data    = single(ncread( datafile, 'soil_temperature_l1',[io(1) ia(1)],[nlo nla]));
data    = data(:);
data    = data(ioa);
data(non_lan) = -999;
datr = data;

sparam  = 'SNTL1';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data


%= Daily soil layer 1 temperature
%  Same considerations for daily values

sparam  = 'DSNTL1';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);

try

  data    = single(ncread( datafile, 'daily_soil_temperature_l1',[io(1) ia(1)],[nlo nla]));
  data    = data(:);
  data    = data(ioa);
  data(non_lan) = -999;

catch

  data = datr;

end

save( outfile, 'data' );
clear data datr


%= Snow melt


data    = single(ncread( datafile, 'snow_melt',[io(1) ia(1)],[nlo nla]));
data    = data(:);
data    = data(ioa);
data(non_lan) = -999;

sparam  = 'SNMLT';
outfile = [ dirout, '/', idfunction, '_Output_',  sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save( outfile, 'data' );
clear data

end

return



