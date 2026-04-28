%-------------------------------------------------------------------------------
%
% MODULE   Forward_Model
%
%    This module does the forward modeling for Earth surfaces.
%    The surface and atmosphere inputs need to be previoulsy extracted 
%    from the netcdf repository by the module GeoInputs_Extract. 
%    
%
% FORMAT   Forward_Model( configurationParameters, inputs, outputs)
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
% Updated:	  2024-07-01
%-------------------------------------------------------------------------------


function Forward_Model( configurationParameters, inputs, outputs)

% Reading from system 
E2E_HOME = getenv( 'E2E_HOME' );
 

%= hard switch to test netcdf saving
%  only for debugging


LOG = Logger();
LOG.info( [' ** E2E_HOME is ', E2E_HOME ])


idfunction = 'Forward_Model';



%= initialize command line parsing class

clp = CLP (configurationParameters, inputs, outputs);


%= Get inputs, outputs and configuration files using

conf1 = clp.getConfFile(1);
conf2 = clp.getConfFile(2);

dirin  = clp.getInputFile(1);
dirout = clp.getOutputFile(1);


LOG.info([ idfunction, ' ** Input folder: ', dirin ])
LOG.info([ idfunction, ' ** Output folder: ', dirout ])


FM_SIMULATION = [ dirout, '/', idfunction ];


% creating folder
LOG.info( [ idfunction, ' ** Creating folder ', dirout ]);
mkdir( dirout );
  



%= Parse configuration files 

cfm1 = ConFM(conf1);
cfm2 = ConFM(conf2);



%= Read parameters

LOG.info( [ idfunction, ' ** Reading parameters from global configuration file']);

geodata_version = cfm1.getParameter('geodata_version').getValue;
LOG.info( [ idfunction, ' ** geodata_version has value ', geodata_version ]);

software_version = cfm1.getParameter('software_version').getValue;
LOG.info( [ idfunction, ' ** software_version has value ', software_version ]);

workers_number = cfm1.getParameter('workers_number').getValue;
LOG.info( [ idfunction, ' ** workers_number has value ', num2str(workers_number) ]);


LOG.info([ idfunction, ' ** Reading parameters from local configuration file']);

frequencies	       = cfm2.getParameter('frequencies').getValue;
bands    	       = cfm2.getParameter('bands').getValue;
und                    = strfind( bands, ' ');

nfr       = length(frequencies);

for a = 1:nfr

  if a == 1
    ind        = 1:und(a);
  elseif a == nfr
    ind = und(a-1)+1:length(bands);
  else
    ind = und(a-1)+1:und(a)-1;
  end
  bnames(a).txt = deblank(bands(ind));  
end



v_polarization	       = cfm2.getParameter('v_polarization').getValue;
h_polarization	       = cfm2.getParameter('h_polarization').getValue;
atmosphere_adding      = cfm2.getParameter('atmosphere_adding').getValue;
zenith_angle	       = cfm2.getParameter('zenith_angle').getValue;
azimuth_angle	       = cfm2.getParameter('azimuth_angle').getValue;
sensor_height	       = cfm2.getParameter('sensor_height').getValue;
data_thinning	       = cfm2.getParameter('data_thinning').getValue;
only_netcdf	       = cfm2.getParameter('only_netcdf').getValue;
only_land_climato      = cfm2.getParameter('only_land_climato').getValue;
only_seaice_climato    = cfm2.getParameter('only_seaice_climato').getValue;



%= checking versions

if ~strcmp( geodata_version,'v2.1') 
  osfi_error('The geodata version given cannot be treated by the software');
end

if ~strcmp( software_version,'v2.1') 
  osfi_error('The software version given do not correspond to this software package');
end


%= checking zenith angle range

iza = unique( zenith_angle(:) );

za_ok = 1;

for f = 1:length(iza)
  if ~(iza(f) >= 0 & iza(f) <= 89)
    za_ok = 0;
    break
  end
end

if ~za_ok
  osfi_error('Zenith angles outside the permitted range');    
end



%= checking frequency range

fre_ok = 1;
afreq  = [1.410 6.925 10.650 18.700 36.500];
agap   = 0.1 * afreq;

for f = 1:nfr

  a =  nearest_in_vector( afreq, frequencies(f) );

  if ~(frequencies(f) >= afreq(a)-agap(a) & frequencies(f) <= afreq(a)+agap(a))
    fre_ok = 0;
    break
  end

end

if ~fre_ok
  osfi_error('Frequencies and outside the permitted range');    
end 

% sizes for storage matrices

nza       = size(zenith_angle,2);
naa       = size(azimuth_angle,2);
nb        = length(bnames);


%= to recover filename

GEOINPUT_SIMULATION = [ dirin, '/GeoInputs_Extract' ];


% datafile name for saving netcdf

outfile = [ GEOINPUT_SIMULATION, '_Output_DATAFILENAME.mat' ];
LOG.info( [ idfunction, ' ** Loading DATAFILENAME' ]);
datafilename = load( outfile );
datafilename = datafilename.data;



if ~only_netcdf





%=== Reading forward model inputs
%    and organizing in structures
%    as required by the FM
%    Reminder:
%    MASK 1 ocean
%	  2 land
%	  3 sea-ice


sparam{1}	 = 'MASK';
sparam{2}	 = 'LAT';
sparam{3}	 = 'LON';
na		 = length( sparam );



for a = 1:na

  outfile = [ GEOINPUT_SIMULATION, '_Output_', sparam{a}, '.mat' ];
  LOG.info( [ idfunction, ' ** Loading ', sparam{a}, ' fields' ]);
  eval( [ 'S.', sparam{a}, ' = load_data_single( outfile );' ]);

end



%=== Checking original lat lons to
%    validate data thinning
%    NOTE: Only for global scenes that have a 
%          regular lon-lat. For polar 
%	   scenes the Sensor_Apply_Antenna
%	   needs to interpolate differently



%= simplifying this for current implementation
%  only option is data thinning same in both
%  lat lon coordinates



if data_thinning >= 1

    outfile = [ GEOINPUT_SIMULATION, '_Output_nonvectorizedLAT', '.mat' ];
    LOG.info( [ idfunction, ' ** Loading nonvectorizedLAT fields for day']);
    nvLAT = load_data_single( outfile );
    na = size(nvLAT,2);
    no = size(nvLAT,1);
    clear nvLAT

    try  
 
      pro = reshape(1:length(S.LAT), no, na);
      pro = pro(1:data_thinning:end,1:data_thinning:end);
      ro  = size( pro,1);
      ra  = size( pro,2);
      pro = pro(:);
      

    catch

      osfi_error('The selected data_thinning does not result in regular lat-lon matrices and cannot be applied');

    end

end



%=== Processing
%    overwriting pro for some debugging tests


if data_thinning == -2

  % for tests, random sampling of
  % 4 surface types and SIC == 1

  pro = [];
  nc = 100;

  if 0 % for SP tests
 
    p = 1;
    sparam{p} = 'SP';
    outfile = [ GEOINPUT_SIMULATION, '_Output_', sparam{p}, '.mat' ];
    LOG.info( [ idfunction, ' ** Loading ', sparam{p}, ' fields for day']);
    aux = load_data_single( outfile );
    aur = reshape( aux, 1500, 1400 );

    p = 1;
    sparam{p} = 'P';
    outfile = [ GEOINPUT_SIMULATION, '_Output_', sparam{p}, '.mat' ];
    LOG.info( [ idfunction, ' ** Loading ', sparam{p}, ' fields for day']);
    aup = load_data_single( outfile );

    pro =  find (aux > 95500 & aux < 96000 ); 

  end
  
  if 1
  % sea
  ind = find( S.MASK == 1);
  if ~isempty(ind)
    ind = ind( randperm(length(ind)));
    pro = [ pro; ind(1:nc )];
  end
  end

  % lan
  if 1
  ind = find( S.MASK == 2);
  if ~isempty(ind)
    ind = ind( randperm(length(ind)));
    pro = [ pro; ind(1:nc )];
  end
  end


  % ice
  if 1
  ind = find( S.MASK == 3);
  if ~isempty(ind)
    ind = ind( randperm(length(ind)));
    pro = [ pro; ind(1:nc )];
  end
  end

  % greenland
  if 0
  ind = find( S.MASK == 2 & S.LAT > 60 & S.LAT < 61 & S.LON < 91 & S.LON > 90);
  whos ind
  if ~isempty(ind)
    pro = [ pro; ind];
  end
  end

  ro  = size( pro,1);
  ra  = size( pro,2);


end


npr       = length(pro);


%=== reducing to the selected cases

S.MASK = S.MASK( pro, :);
S.LON  = S.LON( pro, :);
S.LAT  = S.LAT( pro, :);



%=== reading with data thinning

%= ancillary

a=1;
sparam{a}	 = 'MONTH';
a=a+1;
sparam{a}	 = 'DAY';

np		 = length( sparam );
for p = 1:np

  outfile = [ GEOINPUT_SIMULATION, '_Output_', sparam{p}, '.mat' ];
  LOG.info( [ idfunction, ' ** Loading ', sparam{p}, ' fields for day']);
  aux = load_data_single( outfile );
  if ~(strcmp( sparam{p}, 'DAY') | strcmp( sparam{p}, 'MONTH') | strcmp( sparam{p}, 'P'))
    aux = aux( pro, :);
  end
  eval( [ 'S.', sparam{p}, ' = aux;' ]);

end

clear sparam
a=1;
sparam{a}	 = 'TIME';
p = 1;
outfile = [ GEOINPUT_SIMULATION, '_Output_', sparam{p}, '.mat' ];
LOG.info( [ idfunction, ' ** Loading ', sparam{p}, ' fields for day']);
aux = load_data_single( outfile );
aux = aux( pro, :);
eval( [ 'S.', sparam{p}, ' = aux;' ]);



%=== common Q definitions

%= Sensor

Q.sensor_input.F	  = frequencies;
Q.sensor_input.VPOL	  = v_polarization;
Q.sensor_input.HPOL	  = h_polarization;
Q.sensor_input.H	  = sensor_height;
Q.sensor_input.TOA        = atmosphere_adding;


%= Inputs

Q.dir_input             = [ E2E_HOME, '/InputData/ForwardModelData']; 



%==========================================================================
% IN    atmos_input	structure array		Structure containing the
%						atmospheric profiles:
%
%			   T	[K]		Temperature
%			   P	[Pa]		Pressure
%			   SP   [Pa]		Surface Pressure
%			   Q    [Kg/Kg]		Specific Humidity
%			   CLWC [Kg/Kg]		Specific Cloud Liquid Water Content 
%			   CRWC [Kg/Kg]		Specific Cloud Rain Water Content 	
%
%						and the OPTIONAL altitude-
%						integrated contents
%
%			   TCWV	[kg/m2]		Total Column Water Vapour
%			   TCLW	[kg/m2]		Total Column Liquid Water
%			   TCRW	[kg/m2]		Total Column Rain Water

clear sparam
a=1;
sparam{a}	 = 'P';
a=a+1;
sparam{a}	 = 'SP';
a=a+1;
sparam{a}	 = 'T';
a=a+1;
sparam{a}	 = 'Q';
a=a+1;
sparam{a}	 = 'CLWC';
a=a+1;
sparam{a}	 = 'CRWC';


% NOTE: in the forward modeling of scenes we do not
%       use the integrated contents, as in the forward
%       model code the profiles are adjusted to the
%       column values. That's useful when using the 
%       forward model for the column variable retrievals.
%       But as we modify the atmos profile to be consistent
%       with the surface pressure the column content is not
%       valid any more, and here we do not check in the code
%       that consistency any more. Therefore not imported: 
%        
%	     'TCWV';
%	     'TCLW';
%	     'TCRW';

%


np		 = length( sparam );
for p = 1:np

  outfile = [ GEOINPUT_SIMULATION, '_Output_', sparam{p}, '.mat' ];
  LOG.info( [ idfunction, ' ** Loading ', sparam{p}, ' fields for day']);
  aux = load_data_single( outfile );
  if ~strcmp( sparam{p}, 'P')
    aux = aux( pro, :);
  end
  eval( [ 'A.', sparam{p}, ' = aux;' ]);

end



%= combining CLWW and CRWC as CLRWC
%  as the RT uses a common quantity

A.CLRWC = A.CLWC + A.CRWC;
A = rmfield( A, 'CLWC' );
A = rmfield( A, 'CRWC' );


%= using the full line calculations

A.LOWFREQ = 0;



%==========================================================================
% Sea surface
%			    SST	[K]		Surface Temperature
%			    SSS [g/kg]		Sea Surface Salinity
%			    OWS	[m/s]		Wind velocity
%			    UWS	[m/s]		U wind component
%			    VWS	[m/s]		V wind component


%= notice that sea is also needed for sea-ice

ind = find( S.MASK == 1 | S.MASK == 3);

if ~isempty( ind )

 clear sparam
 a=1;
 sparam{a}	 = 'SST';
 a=a+1;
 sparam{a}	 = 'SSS';
 a=a+1;
 sparam{a}	 = 'OWS';
 a=a+1;
 sparam{a}	 = 'UWS';
 a=a+1;
 sparam{a}	 = 'VWS';

 np		 = length( sparam );
 for p = 1:np

   outfile = [ GEOINPUT_SIMULATION, '_Output_', sparam{p}, '.mat' ];
   LOG.info( [ idfunction, ' ** Loading ', sparam{p}, ' fields for day']);
   aux = load_data_single( outfile );
   aux = aux( pro, :);
   eval( [ 'S.', sparam{p}, ' = aux;' ]);

 end


 %==  adding wind direction 
 %    CONVENTION: angle between the direction toward which the wind is blowing 
 %    and the east direction, counterclockwise positive from due east

 S.OWD = mod(atan2d(S.VWS,S.UWS),360);

end




%==========================================================================
%  Sea ice surface
%			    SIT [K]		Sea ice Surface Temperature
%			    SIC [-]		Sea Ice Concentration 0-1
%			    SITH [m]		Sea Ice Thickness
%			    SIA [years]		Sea Ice Age
%			    SIR [m]		Sea Ice Roughness
%			    SISNTH [m]		Sea Ice Snow Thickness
%			    LAT [degrees]	Latitude, 90S to 90N
%			    LON [degrees]	Latitude, 180W to 180E
%			    TCLW [kg/m2]	Total Column Liquid Water


ind = find( S.MASK == 3 );

if ~isempty( ind )
 
 %=== option -1 can be used to reread emis from a previous
 %    sim for daily simulations. With the current convention 
 %    p00 is the first one and an emis is not existing yet, 
 %    so we overwrite the -1 option if p00. This avoids having separate
 %    control files for p00 and the following pXX

 aux     = strfind( dirout, '/');
 aux     = aux(end)-1;
 pXX     = dirout( aux-2:aux);

 if strcmp( 'p00', pXX ) & only_seaice_climato == -1
   only_seaice_climato = 0;
 end


 if only_seaice_climato == 0

  %= Using SURFEM-SeaIce to derive dynamical emis values

  clear sparam
  a=1;
  sparam{a}	 = 'SIT';
  a=a+1;
  sparam{a}	 = 'DSIT';
  a=a+1;
  sparam{a}	 = 'SIC';
  a=a+1;
  sparam{a}	 = 'SITH';
  a=a+1;
  sparam{a}	 = 'SIA';
  a=a+1;
  sparam{a}	 = 'SIR';
  a=a+1;
  sparam{a}	 = 'SISNTH';
  a=a+1;
  sparam{a}	 = 'DTCWV';

  np		 = length( sparam );
  for p = 1:np

    outfile = [ GEOINPUT_SIMULATION, '_Output_', sparam{p}, '.mat' ];
    LOG.info( [ idfunction, ' ** Loading ', sparam{p}, ' fields for day']);
    aux = load_data_single( outfile );
    aux = aux( pro, :);
    eval( [ 'S.', sparam{p}, ' = aux;' ]);

  end


 elseif only_seaice_climato == 1

   %= climato sea-ice from land climato
        %  used for the moment for regions like the Bering
        %  Sea where there is no emis from the previous one
        %  because nextSIM did not simulate sea-ice
        %  like the Bering Sea. In this case we set pc_sea to zero
        %  as the TELSEM climato emis already includes ocean 
        %  contribution. Notice that the L, C, an X is borrowed 
        %  from the Ka, so the reuslts are likely to be not 
        %  very realistic. Also, the ST we are multiplying
        % may be an issue as well.    

  clear sparam
  a=1;
  sparam{a}	 = 'SIC';
  a=a+1;
  sparam{a}	 = 'SIT';

  np		 = length( sparam );
  for p = 1:np
    outfile = [ GEOINPUT_SIMULATION, '_Output_', sparam{p}, '.mat' ];
    LOG.info( [ idfunction, ' ** Loading ', sparam{p}, ' fields for day']);
    aux = load_data_single( outfile );
    aux = aux( pro, :);
    eval( [ 'S.', sparam{p}, ' = aux;' ]);
  end

  clear sparam

  S.SIMONTH = S.MONTH;


 elseif only_seaice_climato == -1

  %= special case to reuse the sea ice emis from a previous p00
  %  simulation of the same scene

  clear sparam
  a=1;
  sparam{a}	 = 'SIC';
  a=a+1;
  sparam{a}	 = 'SIT';

  np		 = length( sparam );
  for p = 1:np
    outfile = [ GEOINPUT_SIMULATION, '_Output_', sparam{p}, '.mat' ];
    LOG.info( [ idfunction, ' ** Loading ', sparam{p}, ' fields for day']);
    aux = load_data_single( outfile );
    aux = aux( pro, :);
    eval( [ 'S.', sparam{p}, ' = aux;' ]);
  end

  clear sparam

  if strcmp( 'p00', pXX )
    osfi_error( [ 'Reading of p00 simulation sea ice emis demanded but not compatible with current ', pXX, ' simulation' ]);
  end 

  for b = 1:nb
 
    sparam  = [ 'EMIS_ICE_', bnames(b).txt];
    outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
    outfile = strrep( outfile, pXX, 'p00');
    load( [ outfile ]);

    if size( data,1) ~= npr | size(data,3) ~= nza
      osfi_error( [ 'Reading of p00 simulation sea ice emis demanded but not compatible with current ', pXX, ' simulation' ]);
    end 

    if b == 1
      LOG.info( [ idfunction, ' ** Loading EMIS_ICE netcdf fields from p00 simulation' ]);
      S.SIEMIS = single(nan(2,npr,nfr,nza));
    end

    %= saving
    for z = 1:nza
      S.SIEMIS(:,:,b,z) = squeeze(data(:,:,z))';
    end
    clear data

  end

  
 end

end


%==========================================================================
% Land 


clear sparam
a=1;
sparam{a}	 = 'LST';

p = 1;
outfile = [ GEOINPUT_SIMULATION, '_Output_', sparam{p}, '.mat' ];
LOG.info( [ idfunction, ' ** Loading ', sparam{p}, ' fields for day']);
aux = load_data_single( outfile );
aux = aux( pro, :);
eval( [ 'S.', sparam{p}, ' = aux;' ]);




%=== Snowed land
%			    DAY[]		day of month
%			    LTIME [hours]	Decimal local time
%			    T2M [K]		2-meter Air Temperature
%			    LC[]		Lake cover
%			    SDOR[m/km?]		Standard deviation of Orography
%			    SNA[]		Snow Albedo
%			    SND[Kg/m3]		Snow Density
%			    SNTH[m]		Snow thickness
%			    SNTL1[m]		Soil temperature Layer 1
%			    SNT[K]		Snow temperature
%			    SNMLT[m]		Snow melt


ind = find( S.MASK == 2 );

if ~isempty( ind )

 if only_land_climato ==  0

  %=== dynamical emissivity from SURFEM-SnowLand
  %    with climatological values 

  clear sparam
  a=1;
  sparam{a}	 = 'LTIME';
  a=a+1;
  sparam{a}	 = 'T2M';
  a=a+1;
  sparam{a}	 = 'LC';
  a=a+1;
  sparam{a}	 = 'SDOR';
  a=a+1;
  sparam{a}	 = 'SNA';
  a=a+1;
  sparam{a}	 = 'DSNA';
  a=a+1;
  sparam{a}	 = 'SND';
  a=a+1;
  sparam{a}	 = 'SNTH';
  a=a+1;
  sparam{a}	 = 'SNTL1';
  a=a+1;
  sparam{a}	 = 'DSNTL1';
  a=a+1;
  sparam{a}	 = 'SNT';
  a=a+1;
  sparam{a}	 = 'SNMLT';
  a=a+1;
  sparam{a}	 = 'LST';

  np		 = length( sparam );
  for p = 1:np

    outfile = [ GEOINPUT_SIMULATION, '_Output_', sparam{p}, '.mat' ];
    LOG.info( [ idfunction, ' ** Loading ', sparam{p}, ' fields for day']);
    aux = load_data_single( outfile );
    aux = aux( pro, :);
    eval( [ 'S.', sparam{p}, ' = aux;' ]);

  end

 elseif only_land_climato ==  1

  %== special mode using SURFEM-Snowland to only output
  %   the climato emis eia-interpolated. For that SURFEM-SnowLand
  %   expects a field SNT = nan, so we indicate that nd only store
  %   that field

  S.SNT = single(nan( npr, 1 ));

 elseif only_land_climato ==  2

  %===  is a temporary special mode that uses the older land climato
  %     before snowland was developed, and can still be called for
  %     non-snow land bedore a dynamical non-snowed land emis
  %     module is developed, no snow parameters needed 


 end




 %=== extracting climato emis and agb for land emis calculation
 %    if needed to speed caluclations, instead of loading
 %    them internally for each dynamical emis calculation.
 %    Reminder of only_land_climato:
 %
 %  0  is dynamical emissivity from SURFEM-Snowland,   
 %     climato is an input as well so we pre-load here
 %
 %  1  is climatological emis from SURFEM-Snowland, so
 %     we  pre-loaded here.
 %
 %  2  is a temporary special mode that uses the older land climato
 %     before snowland was developed, and can still be called for
 %     non-snow land bedore a dynamical non-snowed land emis
 %     module is developed 


 if only_land_climato == 0 | only_land_climato == 1

    path_rsrc = [ Q.dir_input, '/EmisSnow' ];

    ifreq = Q.sensor_input.F;
    ifreq( ifreq == 1.4135 ) = 1.41;
    nf    = length( ifreq );     
    no    = length( S.LON );     
    ind = find( S.MASK == 2 );

    if ~isempty( ind )

      LOG.info( [ idfunction, ' ** Preparing CLEMIS and CLAGB fields for day']);
 
      ilon = S.LON(ind);
      ilat = S.LAT(ind);
      ilti = S.LTIME(ind);


      %= for robustness we calculate a climato as an average
      %  of the ascending and  descending, if this needs to
      %  be changed used the real ilti, to signal this to
      %  SURFEM-Snowland ilti = []

      alti = [];

      [ cemis, cagb ] = fom_snowland_emis_surfem( ifreq, [], ilat, ilon, S.MONTH, S.DAY, alti, path_rsrc  );

      S.CLEMIS = single(-999*ones(2,no,nf));
      S.CLEMIS( :, ind, : ) = cemis; 

      S.CLAGB  = single(-999*ones(no,1));
      S.CLAGB( ind ) = cagb;

      clear cemis cagb

    end   

 end

end






%= adapting pro
%  jdn=pro left as jnd different from pro
%  could eventually be used in the parallel loop
%  to only do a fraction of pro
spro = pro;
pro  = 1:length(pro);
jnd  = pro;





%= storage matrices


bts          = single(nan( npr, naa, 4, nfr, nza ));

bts_harm_iso = single(nan( npr, 4, nfr, nza ));
bts_harm_cs1 = single(nan( npr, 4, nfr, nza ));
bts_harm_cs2 = single(nan( npr, 4, nfr, nza ));

emis_sea     = single(nan( npr, 2, nfr, nza ));
emis_ice     = single(nan( npr, 2, nfr, nza ));
emis_lan     = single(nan( npr, 2, nfr, nza ));

tra          = single(nan( npr, nfr, nza ));
tbup         = single(nan( npr, nfr, nza ));
tbdown       = single(nan( npr, nfr, nza ));

tcwv         = single(nan( npr, 1));
tclrw        = single(nan( npr, 1 ));




LOG.info( [ idfunction, ' ** Processing ', num2str(length(pro)), ' fields out of ', num2str(npr) ]);


per_pro = round(1:npr);



%=== parallelizing code with 4 cluster

% dummy sentences so parpool can identify
% them as variables


A = A;
S = S;
Q           = Q;
LOG         = LOG;


if workers_number > 1 


  %=== parallel code

  % setting number of workers 

  myCluster = parcluster('local');
  myCluster.NumWorkers = workers_number;  % 'Modified' property now TRUE
  saveProfile(myCluster);    % 'local' profile now updated,
                           % 'Modified' property now FALSE

  % creating pool object

  poolobj = parpool('local', workers_number);

  % parallelized loop


  ico = 1:npr;
  tco = ico(1:round(npr/100):end);


  tic

  parfor co = ico

    if ~isempty(find(per_pro==co)) & ~isempty( find( co==tco ) )
      LOG.info( sprintf(' PROCESSED cell %0.2f making %0.2f per cent \n', co, 100*co/npr ));
    end

    [ y, y_harm_iso, y_harm_cs1, y_harm_cs2, esea, eice, elan, itra, itbup, itbdown, itcwv, itclrw ]  = module_parallel( co, pro, jnd, zenith_angle, azimuth_angle, Q, A, S, LOG );
    
    bts(co,:,:,:,:)         = y;
    bts_harm_iso(co,:,:,:) = y_harm_iso;
    bts_harm_cs1(co,:,:,:) = y_harm_cs1;
    bts_harm_cs2(co,:,:,:) = y_harm_cs2;
    emis_sea(co,:,:,:)     = esea; 
    emis_ice(co,:,:,:)     = eice; 
    emis_lan(co,:,:,:)     = elan;
    tra(co,:,:)            = itra; 
    tbup(co,:,:)           = itbup; 
    tbdown(co,:,:)         = itbdown; 
    tcwv(co)               = itcwv; 
    tclrw(co)              = itclrw; 


  end


  % closing parallel pool
  delete( poolobj );


else

  %=== regular code

  ico = 1:npr;
  tco = ico(1:round(npr/100):end);

  tic

  for co = ico

    if ~isempty(find(per_pro==co)) & ~isempty( find( co==tco ) )
      LOG.info( sprintf(' PROCESSED cell %0.2f making %0.2f per cent \n', co, 100*co/npr ));
    end
     
    [ y, y_harm_iso, y_harm_cs1, y_harm_cs2, esea, eice, elan, itra, itbup, itbdown, itcwv, itclrw ]  = module_parallel( co, pro, jnd, zenith_angle, azimuth_angle, Q, A, S, LOG );
 

    bts(co,:,:,:,:)        = y;
    bts_harm_iso(co,:,:,:) = y_harm_iso;
    bts_harm_cs1(co,:,:,:) = y_harm_cs1;
    bts_harm_cs2(co,:,:,:) = y_harm_cs2;
    emis_sea(co,:,:,:)     = esea; 
    emis_ice(co,:,:,:)     = eice; 
    emis_lan(co,:,:,:)     = elan;
    tra(co,:,:)            = itra; 
    tbup(co,:,:)           = itbup; 
    tbdown(co,:,:)         = itbdown; 
    tcwv(co)               = itcwv; 
    tclrw(co)              = itclrw; 

  end

end


aux  = num2str(toc/60);
LOG.info( [ idfunction, ' ** Processing finished in ', aux , ' minutes']);



%=== Saving the BTs
%    one file per frequency
%    in matlab internal format
%    and also the atmos contributions
%    and 

for f = 1:nb


  data    = bts(:,:,:,f,:);
  sparam  = [ 'BTS_', bnames(f).txt];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
  save_data_single( outfile, data );

  data    = bts_harm_iso(:,:,f,:);
  sparam  = [ 'BTS_HARM_ISO_', bnames(f).txt];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
  save_data_single( outfile, data );

  data    = bts_harm_cs1(:,:,f,:);
  sparam  = [ 'BTS_HARM_CS1_', bnames(f).txt];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
  save_data_single( outfile, data );

  data    = bts_harm_cs2(:,:,f,:);
  sparam  = [ 'BTS_HARM_CS2_', bnames(f).txt];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
  save_data_single( outfile, data );

  data    = emis_sea(:,:,f,:);
  sparam  = [ 'EMIS_SEA_', bnames(f).txt];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
  save_data_single( outfile, data );

  data    = emis_ice(:,:,f,:);
  sparam  = [ 'EMIS_ICE_', bnames(f).txt ];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
  save_data_single( outfile, data );

  data    = emis_lan(:,:,f,:);
  sparam  = [ 'EMIS_LAN_', bnames(f).txt ];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
  save_data_single( outfile, data  );

  data    = tra(:,f,:);
  sparam  = [ 'TRA_', bnames(f).txt ];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
  save_data_single( outfile, data  );

  data    = tbup(:,f,:);
  sparam  = [ 'TBUP_', bnames(f).txt ];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
  save_data_single( outfile, data  );

  data    = tbdown(:,f,:);
  sparam  = [ 'TBDOWN_', bnames(f).txt ];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
  save_data_single( outfile, data  );

end


%=== Saving column integrated H2O and LW

data    = tcwv;
sparam  = [ 'TCWV' ];
outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save_data_single( outfile, data  );

data    = tclrw;
sparam  = [ 'TCLRW' ];
outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save_data_single( outfile, data  );





%=== Saving Surface fields after data thinning
%    with index in PRO


%= PRO index
data    = int32(spro);
sparam  = [ 'IND'];
outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
save_data_single( outfile, data );


gf  = fieldnames( S );
np = length( gf );


for p = 1:np

  sparam  = gf{p};
  data    = getfield( S, sparam );
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields ' ]);
  save_data_single( outfile, data );

end


end













%=== saving BTS in netcdf file for public distribution
%    and input to OSS


%=== netcdf name

ifilesave  = [ dirout, '/', datafilename, '_toa_tbs' ];  
LOG.info( [ idfunction, ' ** Saving ', ifilesave]);

do_gf = 0;

comp(1).name = 'Vpo';
comp(2).name = 'Hpo';
comp(3).name = '3rd';
comp(4).name = '4th';

coml(1).name = 'V-polarized';
coml(2).name = 'H-polarized';
coml(3).name = '3rd Stokes';
coml(4).name = '4th Stokes';



sta.date     = 'to be added';




%= dimensions for reshaping

outfile = [ GEOINPUT_SIMULATION, '_Output_nonvectorizedLAT', '.mat' ];
LOG.info( [ idfunction, ' ** Loading nonvectorizedLAT fields ' ]);
nvLAT = load_data_single( outfile );
na = size(nvLAT,2);
no = size(nvLAT,1);
clear nvLAT



if data_thinning >= 1

  ro  = length( 1:data_thinning:no );
  ra  = length( 1:data_thinning:na );

else

  ro = no;
  ra = na;

end 




%=== filling structure to write netcdf file
%    NOTE: ORDER of non-TB fields is important  


%=== LAT

  sparam  = [ 'LAT'];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Saving ', sparam, ' netcdf fields ' ] );
  load( outfile );
  data = single(reshape(data,ro,ra));

  a = 1;
  str(a).value         = single(data);
  str(a).units         = 'degrees_north';
  str(a).netcdf_name   = 'lat';
  str(a).long_name     = 'latitude';
  str(a).valid_min     = single(-90);
  str(a).valid_max     = single(90);
  str(a).actual_range  = [ min(str(a).value(:)) max(str(a).value(:))];
  str(a).comment       = 'geographic latitude in decimal degrees';

%=== LON

  sparam  = [ 'LON'];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Saving ', sparam, ' netcdf fields ' ]);
  load( outfile );
  data = single(reshape(data,ro,ra));

  data( data > 180 ) = data(data > 180 ) - 360;

  a = a+1;
  str(a).value         = single(data);
  str(a).units         = 'degrees_east';
  str(a).netcdf_name   = 'lon';
  str(a).long_name     = 'longitude' ;
  str(a).valid_min     = single(-180);
  str(a).valid_max     = single(+180);
  str(a).actual_range  = [ min(str(a).value(:)) max(str(a).value(:))];
  str(a).comment       = 'geographic longitude in decimal degrees';
  clear data

%=== incidence angle

  a = a+1;
  str(a).value = single(squeeze(zenith_angle(1,:)));
  str(a).netcdf_name = 'eia';
  str(a).units       = 'degrees';
  str(a).long_name     = 'earth incidence angle';
  str(a).valid_min     = single(0);
  str(a).valid_max     = single(90);
  str(a).actual_range  = [ min(str(a).value(:)) max(str(a).value(:))];
  str(a).comment       = 'angle between the ground normal and the emission direction vector';
  clear data


%==== Geo data

a = 0;


%=== time
  

  sparam  = [ 'TIME'];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Loading ', sparam, ' netcdf fields ' ]);
  load( outfile );
  dato = single(reshape(data,ro,ra));

  a = a+1;
  stx(a).value        = dato;
  stx(a).netcdf_name  = [ 'time' ];
  stx(a).long_name    = [ 'time stamp' ];
  stx(a).units	      = 'days since 2000-01-01 00:00:00 UTC';
  stx(a).valid_min    = single(0);
  stx(a).valid_max    = single(18250);
  stx(a).coordinates  = 'lat lon';
  stx(a).fillvalue    = single(-32768);
  str(a).actual_range  = single([ min(str(a).value(:)) max(str(a).value(:))]);
  stx(a).comment       = 'UTC time stamp for each lat-lon expressed as the fractional number of days since since 2000-01-01 00:00:00';
  clear dato


%=== surface type
%  Original:
%    sea  1
%    land 2
%    ice  3

  sparam  = [ 'MASK'];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Loading ', sparam, ' netcdf fields ' ]);
  load( outfile );

  value = int8( zeros(ro,ra) ) ;
  value( data == 1 ) = int8(1);
  value( data == 2 ) = int8(2);
  value( data == 3 ) = int8(4);

  lsm = data;

  a = a +1;
  stx(a).value         = value;
  stx(a).netcdf_name   = 'surface_type';
  stx(a).long_name     = 'surface type' ;
  stx(a).flag_masks    = int8( [1 2 4]);
  stx(a).flag_meanings = 'water sea-ice land';
  stx(a).valid_min     = int8(1);
  stx(a).valid_max     = int8(4);
  stx(a).fillvalue     = int8(-32768);
  stx(a).actual_range  = [ min(stx(a).value(:)) max(stx(a).value(:))];
  stx(a).coordinates   = 'lat lon';
  stx(a).comment       = 'defined by the original surface variables used to characterize the surface conditions';



%=== Wind direction and wind speed
%    CONVENTION: angle between the direction toward which the wind is blowing 
%    and the east direction, counterclockwise positive from due east


  sparam  = [ 'OWD'];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Loading ', sparam, ' netcdf fields ' ]);
  load( outfile );
  owd = single(reshape(data,ro,ra));

  sparam  = [ 'OWS'];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Loading ', sparam, ' netcdf fields ' ]);
  load( outfile );
  ows = single(reshape(data,ro,ra));



%= wind direction


  ind  = find( ows >= -998 );
  dato = -999 * single( ones( size(ows) ) ); 
  dato(ind) = ows(ind);
  dato( isnan( dato ) == 1) = -999;

  ond       = find( dato < -998 );
  und       = find( dato >= -998 );
  lvalue    = single(nanmin(dato(und)));
  hvalue    = single(nanmax(dato(und)));

  mivalue   = 0;
  mavalue   = 360;

  ran       = hvalue - lvalue;
  ofs       = single(lvalue + ran/2);
  scf       = single(2 * 32000 / ran);
  dato      = int16(scf * (dato - ofs));
  dato(ond) = int16(-32768);    


  a = a+1;
  stx(a).value        = dato;
  stx(a).netcdf_name  = [ 'wind_speed' ];
  stx(a).long_name    = [ 'wind sped' ];
  stx(a).units	      = 'm/s';
  %stx(a).valid_min    = int16(scf*(mivalue-ofs));
  %stx(a).valid_max    = int16(scf*(mavalue-ofs));
  stx(a).coordinates  = 'lat lon';
  stx(a).fillvalue    = int16(-32768);
  stx(a).scale_factor = 1/scf;
  stx(a).add_offset   = ofs;
  stx(a).actual_range  = int16( scf * ([ lvalue hvalue ]-ofs));
  clear dato




%= wind speed

  ind  = find( ows >= -998 );
  dato = -999 * single( ones( size(ows) ) ); 
  dato(ind) = owd(ind);
  dato( isnan( dato ) == 1) = -999;

  ond       = find( dato < -998 );
  und       = find( dato >= -998 );
  lvalue    = single(nanmin(dato(und)));
  hvalue    = single(nanmax(dato(und)));

  mivalue   = 0;
  mavalue   = 360;

  ran       = hvalue - lvalue;
  ofs       = single(lvalue + ran/2);
  scf       = single(2 * 32000 / ran);
  dato      = int16(scf * (dato - ofs));
  dato(ond) = int16(-32768);    


  a = a+1;
  stx(a).value        = dato;
  stx(a).netcdf_name  = [ 'wind_direction' ];
  stx(a).long_name    = [ 'wind direction' ];
  stx(a).units	      = 'degrees';
  %stx(a).valid_min    = int16(scf*(mivalue-ofs));
  %stx(a).valid_max    = int16(scf*(mavalue-ofs));
  stx(a).coordinates  = 'lat lon';
  stx(a).fillvalue    = int16(-32768);
  stx(a).scale_factor = 1/scf;
  stx(a).add_offset   = ofs;
  stx(a).actual_range  = int16( scf * ([ lvalue hvalue ]-ofs));
  stx(a).comment       = 'The wind direction is the defined here as the angle between the direction toward which the wind is blowing and the east direction, counterclockwise positive from due east. The wind direction relative to the earth azimuth angle will be alculated within the instrument simulators as the difference of the wind direction and the earth azimuth agle, where the earth azimuth angle is defined with the same convention, i.e., angle between the ground projection of the emission direction vector and the east direction, counterclockwise positive';
  clear dato


%=== sst

  sparam  = [ 'SST'];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Loading ', sparam, ' netcdf fields ' ]);
  load( outfile );
  data = single(reshape(data,ro,ra));
  data( isnan( data ) == 1) = -999;

  ond       = find( data < -998 );
  und       = find( data >= -998 );
  lvalue    = single(nanmin(data(und)));
  hvalue    = single(nanmax(data(und)));

  mivalue   = 270;
  mavalue   = 310;

  ran       = hvalue - lvalue;
  ofs       = single(lvalue + ran/2);
  scf       = single(2 * 32000 / ran);
  dato      = int16(scf * (data - ofs));
  dato(ond) = int16(-32768);    


  a = a+1;
  stx(a).value        = dato;
  stx(a).netcdf_name  = [ 'sea_surface_temperature' ];
  stx(a).long_name    = [ 'sea surface temperature' ];
  stx(a).units	      = 'K';
  %stx(a).valid_min    = int16(scf*(mivalue-ofs));
  %stx(a).valid_max    = int16(scf*(mavalue-ofs));
  stx(a).coordinates  = 'lat lon';
  stx(a).fillvalue    = int16(-32768);
  stx(a).scale_factor = 1/scf;
  stx(a).add_offset   = ofs;
  stx(a).actual_range  = int16( scf * ([ lvalue hvalue ]-ofs));

  clear data



%=== sss

  sparam  = [ 'SSS'];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Loading ', sparam, ' netcdf fields ' ]);
  load( outfile );
  data = single(reshape(data,ro,ra));
  data( isnan( data ) == 1) = -999;

  ond       = find( data < -998 );
  und       = find( data >= -998 );
  lvalue    = single(nanmin(data(und)));
  hvalue    = single(nanmax(data(und)));

  mivalue   = 0;
  mavalue   = 40;

  ran       = hvalue - lvalue;
  ofs       = single(lvalue + ran/2);
  scf       = single(2 * 32000 / ran);
  dato      = int16(scf * (data - ofs));
  dato(ond) = int16(-32768);    


  a = a+1;
  stx(a).value        = dato;
  stx(a).netcdf_name  = [ 'sea_surface_salinity' ];
  stx(a).long_name    = [ 'sea surface salinity' ];
  stx(a).units	      = 'PSU';
  %stx(a).valid_min    = int16(scf*(mivalue-ofs));
  %stx(a).valid_max    = int16(scf*(mavalue-ofs));
  stx(a).coordinates  = 'lat lon';
  stx(a).fillvalue    = int16(-32768);
  stx(a).scale_factor = 1/scf;
  stx(a).add_offset   = ofs;
  stx(a).actual_range  = int16( scf * ([ lvalue hvalue ]-ofs));

  clear data




%=== tcwv

  sparam  = [ 'TCWV'];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Loading ', sparam, ' netcdf fields ' ]);
  load( outfile );
  data = single(reshape(data,ro,ra));
  data( isnan( data ) == 1) = -999;

  ond       = find( data < -998 );
  und       = find( data >= -998 );
  lvalue    = single(nanmin(data(und)));
  hvalue    = single(nanmax(data(und)));

  mivalue   = 0;
  mavalue   = 100;

  ran       = hvalue - lvalue;
  ofs       = single(lvalue + ran/2);
  scf       = single(2 * 32000 / ran);
  dato      = int16(scf * (data - ofs));
  dato(ond) = int16(-32768);    


  a = a+1;
  stx(a).value        = dato;
  stx(a).netcdf_name  = [ 'total_water_vapour' ];
  stx(a).long_name    = [ 'total colum water vapour' ];
  stx(a).units	      = 'Kg/m2';
  %stx(a).valid_min    = int16(scf*(mivalue-ofs));
  %stx(a).valid_max    = int16(scf*(mavalue-ofs));
  stx(a).coordinates  = 'lat lon';
  stx(a).fillvalue    = int16(-32768);
  stx(a).scale_factor = 1/scf;
  stx(a).add_offset   = ofs;
  stx(a).actual_range  = int16( scf * ([ lvalue hvalue ]-ofs));

  clear data

%=== tclrw

  sparam  = [ 'TCLRW'];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Loading ', sparam, ' netcdf fields ' ]);
  load( outfile );
  data = single(reshape(data,ro,ra));
  data( isnan( data ) == 1) = -999;

  ond       = find( data < -998 );
  und       = find( data >= -998 );
  lvalue    = single(nanmin(data(und)));
  hvalue    = single(nanmax(data(und)));

  mivalue   = 0;
  mavalue   = 1;

  ran       = hvalue - lvalue;
  ofs       = single(lvalue + ran/2);
  scf       = single(2 * 32000 / ran);
  dato      = int16(scf * (data - ofs));
  dato(ond) = int16(-32768);    


  a = a+1;
  stx(a).value        = dato;
  stx(a).netcdf_name  = [ 'total_water_liquid' ];
  stx(a).long_name    = [ 'total colum cloud liquid and rain water' ];
  stx(a).units	      = 'Kg/m2';
  %stx(a).valid_min    = int16(scf*(mivalue-ofs));
  %stx(a).valid_max    = int16(scf*(mavalue-ofs));
  stx(a).coordinates  = 'lat lon';
  stx(a).fillvalue    = int16(-32768);
  stx(a).scale_factor = 1/scf;
  stx(a).add_offset   = ofs;
  stx(a).actual_range  = int16( scf * ([ lvalue hvalue ]-ofs));
  stx(a).comment       = 'This is the sum of the toctal colum cloud liquid water and rain water';
  clear data



  

%=== TB fields in groups per band

pname(1).txt = 'v';
pname(2).txt = 'h';
pname(3).txt = 't3';
pname(4).txt = 't4';


%= TBs ISO at a given zenith angle 
%  only for V and H


ns = 0;

for b = 1:nb
 
  sparam  = [ 'BTS_HARM_ISO_', bnames(b).txt];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  if b == 1
    LOG.info( [ idfunction, ' ** Loading BTS_HARM_ISO netcdf fields ' ]);
  end
  load( [ outfile ]);

  n2 = size( data, 2 );  %eia
  n4 = size( data, 3 );  %pol

  date = single( nan(ro, ra, n2, n4 ) );

  for t2 = 1:n2
      for t4 = 1:n4

        var = reshape(squeeze(data(:,t2,t4)),ro, ra);
        date(:,:,t2,t4) = var;

      end
  end


  for t2 = 1:n2

    dato      = squeeze( date(:,:,t2,:) );
    if do_gf
      dato = gapfillclosest( dato );
    end
    dato( isnan( dato ) == 1) = -999;

    ond       = find( dato < -998 );
    und       = find( dato >= -998 );
    lvalue    = single(nanmin(dato(und)));
    hvalue    = single(nanmax(dato(und)));

    ran       = hvalue - lvalue;
    ofs       = single(lvalue + ran/2);
    scf       = single(2 * 32760 / ran);
    dato      = int16(scf * (dato - ofs));
    dato(ond) = int16(-32768);    

    mivalue   = 0;
    mavalue   = 400;

    stg(b,ns+t2).value        = dato;
    stg(b,ns+t2).netcdf_name  = [ 'brightness_temperature_', pname(t2).txt, '_iso' ];
    stg(b,ns+t2).long_name    = [ bnames(b).txt, '-band ', pname(t2).txt, ' isotropic component of brightness temperature' ];
    stg(b,ns+t2).units	  = 'kelvin';
    %stg(b,ns+t2).valid_min    = int16(scf*(mivalue-ofs));
    %stg(b,ns+t2).valid_max    = int16(scf*(mavalue-ofs));
    stg(b,ns+t2).coordinates  = 'eia lat lon';
    stg(b,ns+t2).fillvalue    = int16(-32768);
    stg(b,ns+t2).scale_factor = 1/scf;
    stg(b,ns+t2).add_offset   = ofs;
    stg(b,ns+t2).actual_range  = int16( scf * ([ lvalue hvalue ]-ofs));

  end  

end



%= TBs 1st harmonic 

ns = 4;

for b = 1:nb
 
  sparam  = [ 'BTS_HARM_CS1_', bnames(b).txt];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  if b == 1
    LOG.info( [ idfunction, ' ** Loading BTS_HARM_CS1 netcdf fields ' ]);
  end
  load( [ outfile ]);

  n2 = size( data, 2 );
  n4 = size( data, 3 );

  date = single( nan(ro, ra, n2, n4 ) );

  for t2 = 1:n2
      for t4 = 1:n4

        var = reshape(squeeze(data(:,t2,t4)),ro, ra);
        date(:,:,t2,t4) = var;

      end
  end

  for t2 = 1:n2

    dato      = squeeze( date(:,:,t2,:) );
    if do_gf
      dato = gapfillclosest( dato );
    end

    dato( isnan( dato ) == 1) = -999;
    if b == 1  
      % L band at 50deg
      dato      = squeeze( dato(:,:,4) );
    else
      % rest band at 55deg
      dato      = squeeze( dato(:,:,5) );
    end

    dato( isnan( dato ) == 1) = -999;
    ond       = find( dato < -998 );
    und       = find( dato >= -998 );
    lvalue    = single(nanmin(dato(und)));
    hvalue    = single(nanmax(dato(und)));

    mivalue   = -5;
    mavalue   = 5;

    ran       = hvalue - lvalue;
    ofs       = single(lvalue + ran/2);
    scf       = single(2 * 32760 / ran);
    dato      = int16(scf * (dato - ofs));
    dato(ond) = int16(-32768);    

    stg(b,ns+t2).value        = dato;
    if t2 <= 2
      stg(b,ns+t2).netcdf_name  = [ 'brightness_temperature_', pname(t2).txt, '_c1' ];
    else
      stg(b,ns+t2).netcdf_name  = [ 'brightness_temperature_', pname(t2).txt, '_s1' ];
    end
    stg(b,ns+t2).long_name    = [ bnames(b).txt, '-band ', pname(t2).txt, ' 1st harmonic component of brightness temperature' ];
    stg(b,ns+t2).units	  = 'kelvin';
    %stg(b,ns+t2).valid_min    = int16(scf*(mivalue-ofs));
    %stg(b,ns+t2).valid_max    = int16(scf*(mavalue-ofs));
    stg(b,ns+t2).coordinates  = 'lat lon';
    stg(b,ns+t2).fillvalue    = int16(-32768);
    stg(b,ns+t2).scale_factor = 1/scf;
    stg(b,ns+t2).add_offset   = ofs;
    stg(b,ns+t2).actual_range = int16( scf * ([ lvalue hvalue ]-ofs));
    if b == 1
      stg(b,ns+t2).comment    = 'These brightness temperatures have a small eia dependence and are only givne for an eia of 50 deg';
    else 
      stg(b,ns+t2).comment    = 'These brightness temperatures have a small eia dependence and are only givne for an eia of 55 deg';
    end

  end  

end





%= TBs 2nd harmonic at a given zenith angle 

ns = 8;

for b = 1:nb
 
  sparam  = [ 'BTS_HARM_CS2_', bnames(b).txt];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  if b == 1
    LOG.info( [ idfunction, ' ** Loading BTS_HARM_CS2 netcdf fields ' ]);
  end
  load( [ outfile ]);

  n2 = size( data, 2 );
  n4 = size( data, 3 );

  date = single( nan(ro, ra, n2, n4 ) );

  for t2 = 1:n2
      for t4 = 1:n4

        var = reshape(squeeze(data(:,t2,t4)),ro, ra);
        date(:,:,t2,t4) = var;

      end
  end

  for t2 = 1:n2

    dato      = squeeze( date(:,:,t2,:) );
    if do_gf
      dato = gapfillclosest( dato );
    end
    if b == 1  
      % L band at 50deg
      dato      = squeeze( dato(:,:,4));
    else
      % rest band at 55deg
      dato      = squeeze( dato(:,:,5));
    end

    dato( isnan( dato ) == 1) = -999;
    ond       = find( dato < -998 );
    und       = find( dato >= -998 );
    lvalue    = single(nanmin(dato(und)));
    hvalue    = single(nanmax(dato(und)));

    mivalue   = -5;
    mavalue   = 5;

    ran       = hvalue - lvalue;
    ofs       = single(lvalue + ran/2);
    scf       = single(2 * 32760 / ran);
    dato      = int16(scf * (dato - ofs));
    dato(ond) = int16(-32768);    


    stg(b,ns+t2).value        = dato;
    if t2 <=2
      stg(b,ns+t2).netcdf_name  = [ 'brightness_temperature_', pname(t2).txt, '_c2' ];
    else
      stg(b,ns+t2).netcdf_name  = [ 'brightness_temperature_', pname(t2).txt, '_s2' ];
    end
    stg(b,ns+t2).long_name    = [ bnames(b).txt, '-band ', pname(t2).txt, ' 2nd harmonic component of brightness temperature' ];
    stg(b,ns+t2).units	  = 'kelvin';
    %stg(b,ns+t2).valid_min    = int16(scf*(mivalue-ofs));
    %stg(b,ns+t2).valid_max    = int16(scf*(mavalue-ofs));
    stg(b,ns+t2).coordinates  = 'lat lon';
    stg(b,ns+t2).fillvalue    = int16(-32768);
    stg(b,ns+t2).scale_factor = 1/scf;
    stg(b,ns+t2).add_offset   = ofs;
    stg(b,ns+t2).actual_range  = int16( scf * ([ lvalue hvalue ]-ofs));
    if b == 1
      stg(b,ns+t2).comment    = 'These brightness temperatures have a small eia dependence and are only givne for an eia of 50 deg';
    else 
      stg(b,ns+t2).comment    = 'These brightness temperatures have a small eia dependence and are only givne for an eia of 55 deg';
    end

  end  

end





%= atmospheric transmittance

ns = 13;

for b = 1:nb
 
  sparam  = [ 'TRA_', bnames(b).txt ];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  if b == 1
    LOG.info( [ idfunction, ' ** Loading TRA netcdf fields ' ]);
  end
  load( [ outfile ]);

  n2 = size( data, 2 );

  dato = single( nan(ro, ra, n2) );

  for t2 = 1:n2

        var = reshape(squeeze(data(:,t2)),ro, ra);
        dato(:,:,t2) = var;

  end

    if do_gf
      dato = gapfillclosest( dato );
    end

    dato( isnan( dato ) == 1) = -999;
    ond       = find( dato < -998 );
    und       = find( dato >= -998 );
    lvalue    = single(nanmin(dato(und)));
    hvalue    = single(nanmax(dato(und)));

    mivalue   = 0;
    mavalue   = 1;

    ran       = hvalue - lvalue;
    ofs       = single(lvalue + ran/2);
    scf       = single(2 * 32760 / ran);
    dato      = int16(scf * (dato - ofs));
    dato(ond) = int16(-32768);    


    stg(b,ns).value        = dato;
    stg(b,ns).netcdf_name  = [ 'transmittance_atmosphere' ];
    stg(b,ns).long_name    = [ bnames(b).txt, '-band atmospheric transmittance' ];
    stg(b,ns).units	  = ' ';
    stg(b,ns).coordinates  = 'lat lon';
    stg(b,ns).fillvalue    = int16(-32768);
    stg(b,ns).scale_factor = 1/scf;
    stg(b,ns).add_offset   = ofs;
    stg(b,ns).actual_range  = int16( scf * ([ lvalue hvalue ]-ofs));


end



%= atmospheric tbup

ns = 14;

for b = 1:nb
 
  sparam  = [ 'TBUP_', bnames(b).txt ];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  if b == 1
    LOG.info( [ idfunction, ' ** Loading TBUP netcdf fields ' ]);
  end
  load( [ outfile ]);

  n2 = size( data, 2 );

  dato = single( nan(ro, ra, n2) );

  for t2 = 1:n2

        var = reshape(squeeze(data(:,t2)),ro, ra);
        dato(:,:,t2) = var;

  end

    if do_gf
      dato = gapfillclosest( dato );
    end
    dato( isnan( dato ) == 1) = -999;
    ond       = find( dato < -998 );
    und       = find( dato >= -998 );
    lvalue    = single(nanmin(dato(und)));
    hvalue    = single(nanmax(dato(und)));

    mivalue   = 0;
    mavalue   = 1;

    ran       = hvalue - lvalue;
    ofs       = single(lvalue + ran/2);
    scf       = single(2 * 32760 / ran);
    dato      = int16(scf * (dato - ofs));
    dato(ond) = int16(-32768);    


    stg(b,ns).value        = dato;
    stg(b,ns).netcdf_name  = [ 'brightness_temperature_atmosphere_up' ];
    stg(b,ns).long_name    = [ bnames(b).txt, '-band atmospheric upwelling brightness_temperature' ];
    stg(b,ns).units	  = 'K';
    stg(b,ns).coordinates  = 'lat lon';
    stg(b,ns).fillvalue    = int16(-32768);
    stg(b,ns).scale_factor = 1/scf;
    stg(b,ns).add_offset   = ofs;
    stg(b,ns).actual_range  = int16( scf * ([ lvalue hvalue ]-ofs));

end



%= atmospheric tbdwon

ns = 15;

for b = 1:nb
 
  sparam  = [ 'TBDOWN_', bnames(b).txt ];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  if b == 1
    LOG.info( [ idfunction, ' ** Loading TBDOWN netcdf fields ' ]);
  end
  load( [ outfile ]);

  n2 = size( data, 2 );

  dato = single( nan(ro, ra, n2) );

  for t2 = 1:n2

        var = reshape(squeeze(data(:,t2)),ro, ra);
        dato(:,:,t2) = var;

  end

    if do_gf
      dato = gapfillclosest( dato );
    end
    dato( isnan( dato ) == 1) = -999;
    ond       = find( dato < -998 );
    und       = find( dato >= -998 );
    lvalue    = single(nanmin(dato(und)));
    hvalue    = single(nanmax(dato(und)));

    mivalue   = 0;
    mavalue   = 1;

    ran       = hvalue - lvalue;
    ofs       = single(lvalue + ran/2);
    scf       = single(2 * 32760 / ran);
    dato      = int16(scf * (dato - ofs));
    dato(ond) = int16(-32768);    


    stg(b,ns).value        = dato;
    stg(b,ns).netcdf_name  = [ 'brightness_temperature_atmosphere_down'];
    stg(b,ns).long_name    = [ bnames(b).txt, '-band atmospheric downwelling brightness_temperature' ];
    stg(b,ns).units	  = 'K';
    stg(b,ns).coordinates  = 'lat lon';
    stg(b,ns).fillvalue    = int16(-32768);
    stg(b,ns).scale_factor = 1/scf;
    stg(b,ns).add_offset   = ofs;
    stg(b,ns).actual_range  = int16( scf * ([ lvalue hvalue ]-ofs));

end



%= Emissivities V and H ocean

ns = 15;

igf = find( lsm == 3 ); % land to nan


for b = 1:nb
 
  sparam  = [ 'EMIS_SEA_', bnames(b).txt];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  if b == 1
    LOG.info( [ idfunction, ' ** Loading EMIS_SEA netcdf fields ' ]);
  end
  load( [ outfile ]);

  n2 = size( data, 2 );
  n4 = size( data, 3 );

  date = single( nan(ro, ra, n2, n4 ) );

  for t2 = 1:n2
      for t4 = 1:n4

        var = reshape(squeeze(data(:,t2,t4)),ro, ra);
        date(:,:,t2,t4) = var;

      end
  end

  for t2 = 1:n2

    dato      = squeeze( date(:,:,t2,:) );
    if do_gf
      dato = gapfillclosest( dato, igf );
    end
    dato( isnan( dato ) == 1) = -999;

    ond       = find( dato < -998 );
    und       = find( dato >= -998 );
    lvalue    = single(nanmin(dato(und)));
    hvalue    = single(nanmax(dato(und)));

    mivalue   = 0;
    mavalue   = 1.5;

    ran       = hvalue - lvalue;
    ofs       = single(lvalue + ran/2);
    scf       = single(2 * 32760 / ran);
    dato      = int16(scf * (dato - ofs));
    dato(ond) = int16(-32768);    

    stg(b,ns+t2).value        = dato;
    stg(b,ns+t2).netcdf_name  = [ 'emissivity_sea_', pname(t2).txt ];
    stg(b,ns+t2).long_name    = [ bnames(b).txt, '-band ', pname(t2).txt, ' sea emissivity' ];
    stg(b,ns+t2).units	  = 'kelvin';
    %stg(b,ns+t2).valid_min    = int16(scf*(mivalue-ofs));
    %stg(b,ns+t2).valid_max    = int16(scf*(mavalue-ofs));
    stg(b,ns+t2).coordinates  = 'lat lon';
    stg(b,ns+t2).fillvalue    = int16(-32768);
    stg(b,ns+t2).scale_factor = 1/scf;
    stg(b,ns+t2).add_offset   = ofs;
    stg(b,ns+t2).actual_range  = int16( scf * ([ lvalue hvalue ]-ofs));

  end  

end




%= Emissivities V and H ice

ns = 17;

igf = find( lsm ~= 2  ); % land and sea

for b = 1:nb
 
  sparam  = [ 'EMIS_ICE_', bnames(b).txt];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  if b == 1
    LOG.info( [ idfunction, ' ** Loading EMIS_ICE netcdf fields ' ]);
  end
  load( [ outfile ]);

  n2 = size( data, 2 );
  n4 = size( data, 3 );

  date = single( nan(ro, ra, n2, n4 ) );

  for t2 = 1:n2
      for t4 = 1:n4

        var = reshape(squeeze(data(:,t2,t4)),ro, ra);
        date(:,:,t2,t4) = var;

      end
  end


  for t2 = 1:n2

    dato      = squeeze( date(:,:,t2,:) );

    if do_gf
      dato = gapfillclosest( dato, igf );
    end

    % to flag the nans in the conversion
    dato( isnan(dato) == 1) = -999;

    ond       = find( dato < -998 );
    und       = find( dato >= -998 );

    if ~isempty( und )

      lvalue    = single(nanmin(dato(und)));
      hvalue    = single(nanmax(dato(und)));

      mivalue   = 0;
      mavalue   = 1.5;

      ran       = hvalue - lvalue;
      ofs       = single(lvalue + ran/2);
      scf       = single(2 * 32760 / ran);
      dato      = int16(scf * (dato - ofs));
      dato(ond) = int16(-32768);    

    else

      mivalue   = 0;
      mavalue   = 1.5;

      ofs       = 1;
      scf       = 0;
      dato      = int16(zeros(size(dato)));
      dato(ond) = int16(-32768);    
 
    end

    stg(b,ns+t2).value        = dato;
    stg(b,ns+t2).netcdf_name  = [ 'emissivity_ice_', pname(t2).txt ];
    stg(b,ns+t2).long_name    = [ bnames(b).txt, '-band ', pname(t2).txt, ' sea emissivity' ];
    stg(b,ns+t2).units	  = 'kelvin';
    %stg(b,ns+t2).valid_min    = int16(scf*(mivalue-ofs));
    %stg(b,ns+t2).valid_max    = int16(scf*(mavalue-ofs));
    stg(b,ns+t2).coordinates  = 'lat lon';
    stg(b,ns+t2).fillvalue    = int16(-32768);
    stg(b,ns+t2).scale_factor = 1/scf;
    stg(b,ns+t2).add_offset   = ofs;
    stg(b,ns+t2).actual_range  = int16( scf * ([ lvalue hvalue ]-ofs));

  end  

end




%= Emissivities V and H lan

ns = 19;

igf = find( lsm  ~= 3  ); % sea-ice and sea

for b = 1:nb
 
  sparam  = [ 'EMIS_LAN_', bnames(b).txt];
  outfile = [ dirout, '/', idfunction, '_Output_', sparam, '.mat' ];  
  if b == 1
    LOG.info( [ idfunction, ' ** Loading EMIS_LAN netcdf fields ' ]);
  end
  load( [ outfile ]);

  n2 = size( data, 2 );
  n4 = size( data, 3 );

  date = single( nan(ro, ra, n2, n4 ) );

  for t2 = 1:n2
      for t4 = 1:n4

        var = reshape(squeeze(data(:,t2,t4)),ro, ra);
        date(:,:,t2,t4) = var;

      end
  end

  for t2 = 1:n2

    dato      = squeeze( date(:,:,t2,:) );

    if do_gf
      dato = gapfillclosest( dato, igf );
    end

    % to flag the nans in the conversion
    dato( isnan(dato) == 1) = -999;

    ond       = find( dato < -998 );
    und       = find( dato >= -998 );

    if ~isempty( und )

      lvalue    = single(nanmin(dato(und)));
      hvalue    = single(nanmax(dato(und)));

      mivalue   = 0;
      mavalue   = 1.5;

      ran       = hvalue - lvalue;
      ofs       = single(lvalue + ran/2);
      scf       = single(2 * 32760 / ran);
      dato      = int16(scf * (dato - ofs));
      dato(ond) = int16(-32768);    

    else

      mivalue   = 0;
      mavalue   = 1.5;

      ofs       = 1;
      scf       = 0;
      dato      = int16(zeros(size(dato)));
      dato(ond) = int16(-32768);    
 
    end

    stg(b,ns+t2).value        = dato;
    stg(b,ns+t2).netcdf_name  = [ 'emissivity_lan_', pname(t2).txt ];
    stg(b,ns+t2).long_name    = [ bnames(b).txt, '-band ', pname(t2).txt, ' sea emissivity' ];
    stg(b,ns+t2).units	  = 'kelvin';
    %stg(b,ns+t2).valid_min    = int16(scf*(mivalue-ofs));
    %stg(b,ns+t2).valid_max    = int16(scf*(mavalue-ofs));
    stg(b,ns+t2).coordinates  = 'lat lon';
    stg(b,ns+t2).fillvalue    = int16(-32768);
    stg(b,ns+t2).scale_factor = 1/scf;
    stg(b,ns+t2).add_offset   = ofs;
    stg(b,ns+t2).actual_range  = int16( scf * ([ lvalue hvalue ]-ofs));

  end  

end





forward_model_netcdf_writing( sta, stx, str, stg, [ ifilesave, '.nc' ] );


return



%==========================================================================


function [ bts, bts_harm_iso, bts_harm_cs1, bts_harm_cs2, esea, eice, elan, tra, tbup, tbdown, tcwv, tclrw ]  = module_parallel( co, pro, dopro, zenith_angle, azimuth_angle, Q, A, S, LOG );




%=== Frequencies
%

afreq	= [1.410 6.925 10.650 18.700 36.500];
nf      = length( Q.sensor_input.F );
nza     = size(zenith_angle,2);
nc      = 4;


p = pro(co);

nza    = length(zenith_angle);
naa    = length(azimuth_angle);
nch    = length(Q.sensor_input.F);
nfr    = length(Q.sensor_input.F);
npo    = 4;


%= initializing


%= zenith angle

Q.sensor_input.ZA = zenith_angle;


%= signaling no colums in the RT calculation

atmos_input.P     = A.P;
atmos_input.TCWV  = [];
atmos_input.TCLRW = [];


%= individual structures for RT script

atmos_input.LOWFREQ = A.LOWFREQ;
atmos_input.P       = A.P;

if isfield( S, 'MONTH' )
  surf_input.MONTH = S.MONTH;
  S = rmfield( S, 'MONTH' );
end

if isfield( S, 'DAY' )
  surf_input.DAY = S.DAY;
  S = rmfield( S, 'DAY' );
end

if isfield( S, 'SIMONTH' )
  surf_input.SIMONTH = S.SIMONTH;
  S = rmfield( S, 'SIMONTH' );
end





Afn = fieldnames( A );
naf = length( Afn );
Sfn = fieldnames( S );
nsf = length( Sfn );


if sum(dopro==p) >0 

  atmos_input.SP    = A.SP(p); 
  atmos_input.T     = A.T(p,:)'; 
  atmos_input.Q     = A.Q(p,:)'; 
  atmos_input.CLRWC = A.CLRWC(p,:)'; 

  for a = 1:nsf
    if ~strcmp( Sfn{a}, 'CLEMIS' ) & ~strcmp( Sfn{a}, 'SIEMIS' )
      eval( ['surf_input.', Sfn{a}, ' = S.', Sfn{a}, '(p,:);' ]);
    end
  end    

  % CLEMIS if exist , size is 2 x npro x freq
  if isfield( S, 'CLEMIS' )
    surf_input.CLEMIS = squeeze(S.CLEMIS(:,p,:));
  end

  % SIEMIS if exist , size is 2 x npro x freq x za
  if isfield( S, 'SIEMIS' )
    surf_input.SIEMIS = squeeze(S.SIEMIS(:,p,:,:));
  end


  %try

    %=== simulating harmonics so AA does not matter, set to zero
    Q.sensor_input.AA = 0;

    [~, bts_harm_iso, bts_harm_cs1, bts_harm_cs2, esea, eice, elan, tra, tbup, tbdown, tcwv, tclrw ]= fom_atmos_ocean_seaice_land_tbs_core( atmos_input, surf_input, Q.sensor_input, Q.dir_input );


    %=== simuating AA dependent TBs using the t3 and t4 harmonics
    %    only stored in internal mat format, in the netcdf file
    %    we only give away the components as there is no 
    %    azimuthal dependence

    bts          = nan( naa, 4, nfr, nza );

    for aa = 1:naa

      rowd = deg2rad(surf_input.OWD - azimuth_angle(aa));

      %=== including azimuthal dependence for v and pol
      %    but not used, see above 
    
      bts(aa,1,:,:) = bts_harm_iso(1,:,:) + bts_harm_cs1(1,:,:) *  cos(rowd) + bts_harm_cs2(1,:,:) *  cos(2*rowd); 
      bts(aa,2,:,:) = bts_harm_iso(2,:,:) + bts_harm_cs1(2,:,:) *  cos(rowd) + bts_harm_cs2(2,:,:) *  cos(2*rowd); 
      bts(aa,3,:,:) = bts_harm_cs1(3,:,:) *  sin(rowd) + bts_harm_cs2(3,:,:) *  sin(2*rowd); 
      bts(aa,4,:,:) = bts_harm_cs1(4,:,:) *  sin(rowd) + bts_harm_cs2(4,:,:) *  sin(2*rowd); 

    end


  %catch me


    if 0
    LOG.info( sprintf(' CELL %0.2f DID NOT PROCESS DUE TO ERROR \n', p ));
    LOG.info( [ 'ERROR message: ', me.message]);

    bts          = nan( 4, nfr, nza );
    bts_harm_iso = nan( 4, nfr, nza );
    bts_harm_cs1 = nan( 4, nfr, nza );
    bts_harm_cs2 = nan( 4, nfr, nza );
    esea         = nan( 2, nfr, nza );
    eice         = nan( 2, nfr, nza );
    elan         = nan( 2, nfr, nza );
    tra          = nan( nfr, nza );
    tbup         = nan( nfr, nza );
    tbdown       = nan( nfr, nza );
    tcwv         = nan( nfr, nza );
    tclrw        = nan( nfr, nza );
    end

  %end

else

  LOG.info( sprintf(' CELL %0.2f NOT TO BE PROCESSED \n', p ));

end



return




%=======================================================================

function forward_model_netcdf_writing( sta, stx, str, stg, savename )

eval( [ '!rm ', savename ]);

%=== setting no compression

dflevel  = 1; 


%=== removing file if exists

if exist( savename, 'file' )
  delete( savename );
end



%=== checking sizes


if strcmp(str(1).netcdf_name,'lat')
  ny = size(str(1).value,1);
  nx = size(str(1).value,2);
else
  error('Latitude should be the first cell in str');
end
if ~strcmp(str(2).netcdf_name,'lon')
  error('Longitude should be the second cell in str');
end
if strcmp(str(3).netcdf_name,'eia')
  nz = length( str(3).value);
else
  error('Earth incidence angle should be the 3rd cell in str');
end

%=== writing coordinate variables


for f = 1:length(str)

    itype = class( str(f).value );

    if f == 3

      nccreate( savename, str(f).netcdf_name, 'Dimensions', { 'eia' nz }, 'DataType', itype,  'DeflateLevel', dflevel, 'Format', 'netcdf4'    );

    else

      nccreate( savename, str(f).netcdf_name, 'Dimensions', { 'y' ny 'x' nx }, 'DataType', itype,  'DeflateLevel', dflevel, 'Format', 'netcdf4');

    end

    ncwrite( savename, str(f).netcdf_name, str(f).value);

    if ~isempty( str(f).actual_range )
      ncwriteatt( savename, str(f).netcdf_name, 'actual_range', str(f).actual_range );
    end
 
    ncwriteatt( savename, str(f).netcdf_name, 'long_name', str(f).long_name );

    ncwriteatt( savename, str(f).netcdf_name, 'standard_name', str(f).netcdf_name );

    ncwriteatt( savename, str(f).netcdf_name, 'units', str(f).units );

    ncwriteatt( savename, str(f).netcdf_name, 'comment', str(f).comment );


end



%=== writing ancillary variables 
%    in group



for f = 1:length(stx)

    itype = class( stx(f).value );

    varname = ['/AUXILIARY/', stx(f).netcdf_name ];

    nccreate( savename, varname, 'Dimensions', { 'y' ny 'x' nx }, 'DataType', itype,  'DeflateLevel', dflevel, 'Format', 'netcdf4', 'FillValue', stx(f).fillvalue);

    ncwrite( savename, varname, stx(f).value );

    if ~isempty( stx(f).add_offset )
      ncwriteatt( savename, varname, 'add_offset', stx(f).add_offset );
    end

    if ~isempty( stx(f).scale_factor )
      ncwriteatt( savename, varname, 'scale_factor', stx(f).scale_factor );
    end

    if 0
    if ~isempty( stx(f).valid_min )
      ncwriteatt( savename, varname, 'valid_min', stx(f).valid_min );
    end

    if ~isempty( stx(f).valid_max )
      ncwriteatt( savename, varname, 'valid_max', stx(f).valid_max );
    end
    end

    if ~isempty( stx(f).coordinates )
      ncwriteatt( savename, varname, 'coordinates', stx(f).coordinates );
    end

    if ~isempty( stx(f).actual_range )
      ncwriteatt( savename, varname, 'actual_range', stx(f).actual_range );
    end
 
    ncwriteatt( savename, varname, 'long_name', stx(f).long_name );

    ncwriteatt( savename, varname, 'standard_name', stx(f).netcdf_name );

    ncwriteatt( savename, varname, 'units', stx(f).units );


end




%=== writing TB variables
%    in groups


bname(1).txt = 'L_BAND';
bname(2).txt = 'C_BAND';
bname(3).txt = 'X_BAND';
bname(4).txt = 'Ku_BAND';
bname(5).txt = 'Ka_BAND';


nb = size(stg,1);
nv = size(stg,2);

%= index for not writing _iso for t3 and t4

indv  = [1 2   5 6 7 8   9 10 11 12 13 14 15 16 17 18 19 20 21];


for b = 1:nb

  for v = indv

    itype = class( stg(b,v).value );

    varname = ['/', bname(b).txt, '/', stg(b,v).netcdf_name ];

    if size(stg(b,v).value,3) > 1   

      nccreate( savename, varname, 'Dimensions', { 'y' ny 'x' nx 'eia' nz }, 'DataType', itype,  'DeflateLevel', dflevel, 'Format', 'netcdf4', 'FillValue',  stg(b,v).fillvalue    );

      eval( [ 'auf = ', itype, '( zeros( ny, nx, nz) );' ])
      for z = 1:nz
        auk = stg(b,v).value;
        auk = squeeze(auk(:,:,z));
        auf(:,:,z) = auk;
      end 

      ncwrite( savename, varname, auf );

    else

      nccreate( savename, varname, 'Dimensions', { 'y' ny 'x' nx }, 'DataType', itype,  'DeflateLevel', dflevel, 'Format', 'netcdf4', 'FillValue',  stg(b,v).fillvalue    );

      ncwrite( savename, varname, stg(b,v).value );

    end




    if ~isempty( stg(b,v).add_offset )
      ncwriteatt( savename, varname, 'add_offset', stg(b,v).add_offset );
    end

    if ~isempty( stg(b,v).scale_factor )
      ncwriteatt( savename, varname, 'scale_factor', stg(b,v).scale_factor );
    end

    if 0
    if ~isempty( stg(b,v).valid_min )
      ncwriteatt( savename, varname, 'valid_min', stg(b,v).valid_min );
    end

    if ~isempty( stg(b,v).valid_max )
      ncwriteatt( savename, varname, 'valid_max', stg(b,v).valid_max );
    end
    end

    if ~isempty( stg(b,v).coordinates )
      ncwriteatt( savename, varname, 'coordinates', stg(b,v).coordinates );
    end

    if ~isempty( stg(b,v).actual_range )
      ncwriteatt( savename, varname, 'actual_range', stg(b,v).actual_range );
    end
 
    ncwriteatt( savename, varname, 'long_name', stg(b,v).long_name );

    ncwriteatt( savename, varname, 'standard_name', stg(b,v).netcdf_name );

    ncwriteatt( savename, varname, 'units', stg(b,v).units );

  end

end





%=== writing global attributes

ncwriteatt( savename, '/', 'title', [ 'Top-of-atmosphere brigthness temperatures for geophysical scene' ] );
ncwriteatt( savename, '/', 'project', 'ESA CIMR-SCEPS-II' );
ncwriteatt( savename, '/', 'summary', 'Top-of-atmosphere brigthness temperatures produced with CIMR-SCEPS Work Package 2.1'  );
ncwriteatt( savename, '/', 'description', 'Top-of-atmosphere brigthness temperatures for geophysical scene'  );
ncwriteatt( savename, '/', 'work_package_manager_name', 'Carlos Jimenez' );
ncwriteatt( savename, '/', 'work_package_manager_email', 'carlos.jimenez@estellus.fr' );
ncwriteatt( savename, '/', 'source', 'CIMR Products Test Data Sets' );
ncwriteatt( savename, '/', 'grid_name', 'EASE2_N01KM' );
ncwriteatt( savename, '/', 'scene_type', 'harmonics' );
ncwriteatt( savename, '/', 'acknowledgement' , 'Data providers are acknowledged by making their data available' ); 
ncwriteatt( savename, '/', 'creator_name', 'Carlos Jimenez, Estellus' );
ncwriteatt( savename, '/', 'creator_email', 'carlos.jimenez@estellus.fr' );
ncwriteatt( savename, '/', 'creator_url', 'www.estellus.fr' );
ncwriteatt( savename, '/', 'history', [ string(datetime) ] );




return


