%-------------------------------------------------------------------------------
%
% MODULE  Orbit_Geolocation_Extract
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
% Initiated:	  2024-03-25
%-------------------------------------------------------------------------------

function Orbit_Geolocation_Extract( configurationParameters, inputs, outputs)


idfunction = 'Orbit_Geolocation_Extract';


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


%= Parse configuration files 

cfm1 = ConFM(conf1);
cfm2 = ConFM(conf2);

%= Read parameters

LOG.info( [ idfunction, ' ** Reading parameters from global configuration file']);

geodata_version = cfm1.getParameter('geodata_version').getValue;
LOG.info( [ idfunction, ' ** geodata_version has value ', geodata_version ]);

software_version = cfm1.getParameter('software_version').getValue;
LOG.info( [ idfunction, ' ** software_version has value ', software_version ]);

SCENE_FILE_NAME = clp.getInputFile(1);
LOG.info( [ idfunction, ' ** Scene TOA-TBs extracted from file ', SCENE_FILE_NAME, ' scene' ]);

ORBIT_FILE_NAME = clp.getInputFile(2);
LOG.info( [ idfunction, ' ** Orbit information from internal file ', ORBIT_FILE_NAME ]);


%=== Some checks

if ~strcmp( geodata_version,'v1.3') 
  osfi_error('The geodata version given cannot be treated by the software');
end

if ~strcmp( software_version,'v1.3') 
  osfi_error('The software version given do not correspond to this software package');
end


LOG.info( [ idfunction, ' ** Orbit information from internal file -2 ', ORBIT_FILE_NAME ]);

%=== Reading scene geolocation 

ifile = SCENE_FILE_NAME;

if  ~contains( ifile, '.nc' ) 
  ifile = [ ifile, '.nc' ];
end

lon  = ncread( ifile, 'longitude' );
lat  = ncread( ifile, 'latitude' );

%=== Internal file with orbit info

ofile = ORBIT_FILE_NAME;
if  ~contains( ofile, '.mat' ) 
  ofile = [ ofile, '.mat' ];
end

LOG.info( [ idfunction, ' ** Orbit information from internal file -4 ', ORBIT_FILE_NAME ]);

%=== Locate orbit section corresponding to given scene


clon = mean( lon(:) );
clat = mean( lat(:) );


% Find closest asc and desc pass center times in MATLAB datenum() format.
[date_asc_orbit, date_des_orbit] = cimr_find_closest_time_latlon(clat, clon, ofile);

sprintf('asc = %s %s :: desc = %s %s\n', date_asc_orbit, datestr(date_des_orbit) );
LOG.info( [ idfunction, ' ** Finding simulation dates in file orbit' ] );


%=== Saving orbit center dates of simulation

filesave = [ dirout, '/', idfunction, '_Output_date_asc_orbit.asc' ];  
LOG.info( [ idfunction, ' ** Saving ascii file ', filesave ] );
save(filesave,'date_asc_orbit','-ascii','-double');

filesave = [ dirout, '/', idfunction, '_Output_date_des_orbit.asc' ];  
LOG.info( [ idfunction, ' ** Saving ascii file ', filesave ] );
save(filesave,'date_des_orbit','-ascii','-double');


return


%==========================================================================
%  Finding closes lon-lat-time in orbit file

function [ta, td, ta0, ta1, td0, td1] = cimr_find_closest_time_latlon(clat, clon, orbit_file)

o = load( orbit_file);

gia = find(o.output.sat_vel_north > 0);                                          
gid = find(o.output.sat_vel_north < 0);

%=keyboard

if 0

  %= using matlab map toolbox

  tda = distance(clat, clon, o.output.sp_lat_geod(gia), o.output.sp_lon_geod(gia));
  tdd = distance(clat, clon, o.output.sp_lat_geod(gid), o.output.sp_lon_geod(gid));

else

   %= using free toolbox m_map
   
   ng           = length(gia);

   lon          = nan(1,2*ng);
   lon(1:2:end) = clon;
   lon(2:2:end) = o.output.sp_lon_geod(gia);

   lat          = nan(1,2*ng);
   lat(1:2:end) = clat;
   lat(2:2:end) = o.output.sp_lat_geod(gia);
  
   tda = m_lldist(lon,lat);
   tda = tda(1:2:end);

   ng           = length(gid);

   lon          = nan(1,2*ng);
   lon(1:2:end) = clon;
   lon(2:2:end) = o.output.sp_lon_geod(gid);

   lat          = nan(1,2*ng);
   lat(1:2:end) = clat;
   lat(2:2:end) = o.output.sp_lat_geod(gid);
  
   tdd = m_lldist(lon,lat);
   tdd = tdd(1:2:end);

end


dtwin = 10*60/86400;

[td, mia] = min(tda);
[td, mid] = min(tdd);
mia = gia(mia);
mid = gid(mid);
ta0 = o.output.time(mia)-dtwin/2;
ta1 = o.output.time(mia)+dtwin/2;
td0 = o.output.time(mid)-dtwin/2;
td1 = o.output.time(mid)+dtwin/2;

ta = o.output.time(mia);
td = o.output.time(mid);

gia = find((o.output.time >= ta0) & (o.output.time <= ta1));
gid = find((o.output.time >= td0) & (o.output.time <= td1));


return


%==========================================================================
function [dist,lons,lats] = m_lldist(long,lat,N)
% M_LLDIST Spherical earth distance between points in long/lat coordinates. 
%   RANGE=M_LLDIST(LONG,LAT) gives the distance in kilometers between
%   successive points in the vectors LONG and LAT, computed
%   using the Haversine formula on a spherical earth of radius
%   6378.137km. Distances are probably good to better than 1% of the
%   "true" distance on the ellipsoidal earth
%
%   [RANGE,LONGS,LATS]=M_LLDIST(LONG,LAT,N) computes the N-point geodesics
%   between successive points. Each geodesic is returned on its
%   own row of length N+1.
%
%   See also M_XYDIST

% Rich Pawlowicz (rich@ocgy.ubc.ca) 6/Nov/00
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.
%
% 30/Dec/2005 - added n-point geodesic computations, based on an algorithm
%               coded by Jeff Barton at Johns Hopkins APL in an m-file
%               I looked at at mathworks.com.


pi180=pi/180;
earth_radius=6378.137;

m=length(long)-1;

long1=reshape(long(1:end-1),m,1)*pi180;
long2=reshape(long(2:end)  ,m,1)*pi180;
lat1= reshape(lat(1:end-1) ,m,1)*pi180;
lat2= reshape(lat(2:end)   ,m,1)*pi180;

dlon = long2 - long1; 
dlat = lat2 - lat1; 
a = (sin(dlat/2)).^2 + cos(lat1) .* cos(lat2) .* (sin(dlon/2)).^2;
angles = 2 * atan2( sqrt(a), sqrt(1-a) );
dist = earth_radius * angles;


if nargin==3 && nargout>1   % Compute geodesics.

  % Cartesian unit vectors in rows of v1,v2
  v1=[cos(long1).*cos(lat1)   sin(long1).*cos(lat1)   sin(lat1) ];
  v2=[cos(long2).*cos(lat2)   sin(long2).*cos(lat2)   sin(lat2) ];

  % We want to get a unit vector tangent to the great circle.
  n1=cross(v1,v2,2); 
  t1=cross(n1,v1,2);
  t1=t1./repmat(sqrt(sum(t1.^2,2)),1,3);

  lons=zeros(m,N+1);
  lats=zeros(m,N+1);
  for k=1:m

   % Radials for all points
   p1=v1(k,:)'*cos(angles(k)*(0:N)/N) + t1(k,:)'*sin(angles(k)*(0:N)/N);

   lons(k,:)=atan2(p1(2,:),p1(1,:))/pi180;
   lats(k,:)=asin(p1(3,:))/pi180;

  end

end


