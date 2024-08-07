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

do_just_netcdf = 0;



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

for a = 1:length(frequencies)

  if a == 1
    ind        = 1:und(a);
  elseif a == length(frequencies)
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
no_land		       = cfm2.getParameter('no_land').getValue;
emis_land_external     = cfm2.getParameter('emis_land_external').getValue;


%= removed partameters from previous versions
%  left as reminder of previous capabilities
%  in case they need to be put back in
do_rem = 1;
if ~do_rem
  emis_ice_smoothing     = cfm2.getParameter('emis_ice_smoothing').getValue;
end


%= checking versions

if ~strcmp( geodata_version,'v1.3') 
  osfi_error('The geodata version given cannot be treated by the software');
end

if ~strcmp( software_version,'v1.3') 
  osfi_error('The software version given do not correspond to this software package');
end


%= checking ice calculation

if ~do_rem
  if emis_ice_smoothing ~= 1
    osfi_error('so far emis_ice_smoothing can only bet set to 1');    
  end
end

%= checking consistency frequencies and angles


if size( zenith_angle, 1) ~= length( frequencies )
  osfi_error('Frequencies and zenith angles are not consistent');    
end

if size( azimuth_angle, 1) ~= length( frequencies )
  osfi_error('Frequencies and zenith angles are not consistent');    
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


%= checking that zenith_angles are the same 
%  for all freqs, to be consistne with
%  current emis angle interpolation

if 0
if sum(sum(zenith_angle,1)/size(zenith_angle,1)) ~= sum(zenith_angle(1,:),2)
  osfi_error('Zenith angles need to be the same for all frequencies');
end 
end   



%= checking frequency range

fre_ok = 1;
afreq  = [1.410 6.925 10.650 18.700 36.500];
agap   = 0.1 * afreq;

for f = 1:length(frequencies)

  a =  nearest_in_vector( afreq, frequencies(f) );

  if ~(frequencies(f) >= afreq(a)-agap(a) & frequencies(f) <= afreq(a)+agap(a))
    fre_ok = 0;
    break
  end

end

if ~fre_ok
  osfi_error('Frequencies and outside the permitted range');    
end 

% angle an frequency sizes

nza       = size(zenith_angle,2);
naa       = size(azimuth_angle,2);
nb        = length(bnames);


%= recovering filename

GEOINPUT_SIMULATION = [ dirin, '/GeoInputs_Extract' ];
outfile    = [ GEOINPUT_SIMULATION, '_Output_*MASK.mat' ];
aux        = dir(outfile);
auc        = strfind( aux.name, '_MASK.mat');
SCENE_DATE = aux.name(auc-8:auc-1);




if ~do_just_netcdf



%=== Reading forward model inputs
%    and organizing in structures
%    as required by the FM


sparam{1}	 = 'MASK';
sparam{2}	 = 'LAT';
sparam{3}	 = 'LON';
np		 = length( sparam );



for p = 1:np

  outfile = [ GEOINPUT_SIMULATION, '_Output_', SCENE_DATE, '_', sparam{p}, '.mat' ];
  LOG.info( [ idfunction, ' ** Loading ', sparam{p}, ' fields for day ', SCENE_DATE ]);
  eval( [ 'I.', sparam{p}, ' = load_data_single( outfile );' ]);

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

    outfile = [ GEOINPUT_SIMULATION, '_Output_', SCENE_DATE, '_nonvectorizedLAT', '.mat' ];
    LOG.info( [ idfunction, ' ** Loading nonvectorizedLAT fields for day ', SCENE_DATE ]);
    nvLAT = load_data_single( outfile );
    na = size(nvLAT,2);
    no = size(nvLAT,1);
    clear nvLAT

    try  
 
      pro = reshape(1:length(I.LAT), no, na);
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
  nc = 10;



  if 1
  % sea
  ind = find( I.MASK == 1);
  if ~isempty(ind)
    %ind = ind( randperm(length(ind)));
    pro = [ pro; ind(1:nc )];
  end
  end

  % land
  if 1
  ind = find( I.MASK == 3);
  if ~isempty(ind)
    %ind = ind( randperm(length(ind)));
    pro = [ pro; ind(1:nc )];
  end
  end

  if 1
  % ice
  ind = find( I.MASK == 2);
  if ~isempty(ind)
    %ind = ind( randperm(length(ind)));
    pro = [ pro; ind(1:nc )];
  end
  end




  if 0
  % coast
  ind = find( I.MASK == 4);
  if ~isempty(ind)
    %ind = ind( randperm(length(ind)));
    if length(ind) > nc
      pro = [ pro; ind(1:nc )]';
    else
      pro = [ pro; ind]';
    end 
  end
  end

  ro  = size( pro,1);
  ra  = size( pro,2);


end



%=== reading with data thinning


sparam{1}	 = 'P';
sparam{2}	 = 'IST';
sparam{3}	 = 'IW';
sparam{4}	 = 'LAT';
sparam{5}	 = 'LLP';
sparam{6}	 = 'LON';
sparam{7}	 = 'LST';
sparam{8}	 = 'LW';
sparam{9}	 = 'MONTH';
sparam{10}	 = 'OWS';
sparam{11}	 = 'UWS';
sparam{12}	 = 'VWS';
sparam{13}	 = 'ILP';
sparam{14}	 = 'RHO';
sparam{15}	 = 'SIC';
sparam{16}	 = 'SND';
sparam{17}	 = 'SLP';
sparam{18}	 = 'SSS';
sparam{19}	 = 'SST';
sparam{20}	 = 'T';
sparam{21}	 = 'T2M';
sparam{22}	 = 'XYI';
sparam{23}	 = 'MASK';
sparam{24}	 = 'FDN';
sparam{25}	 = 'CWVC';
sparam{26}	 = 'SIT';

% NOTE: in the forward modeling of scenes we do not
%       use the integrated contents, as in the forward
%       model code the profiles are adjusted to the
%       column values. That's useful when using the 
%       forward model for the column variable retrievals.
%       But as we modify the atmos profile to be consistent
%       with the surface pressure the column content is not
%       valid any more, and here we do not check in the code
%       that consistency any more. Therefore change: 
%        
%	     sparam{23}	 = 'CWVC';
%	     sparam{24}	 = 'CLWC';

np		 = length( sparam );



for p = 1:np

  outfile = [ GEOINPUT_SIMULATION, '_Output_', SCENE_DATE, '_', sparam{p}, '.mat' ];
  LOG.info( [ idfunction, ' ** Loading ', sparam{p}, ' fields for day ', SCENE_DATE ]);
  aux = load_data_single( outfile );
  if p > 1
    aux = aux( pro, :);
  end
  eval( [ 'I.', sparam{p}, ' = aux;' ]);

end


% adapting pro
spro = pro;
pro  = 1:length(pro);


%= index to make land or not

if no_land
  % we only do sea, ice, and coast
  jnd = find(I.MASK ~= 2 );
else
  jnd = pro;
end



%= storage matrices

npr       = length(pro);
% v + h + the 3rd and 4th Stokes for each freq
nfr       = sum(v_polarization) + sum(h_polarization) + 2*length(frequencies);

bts          = nan( npr, nza, naa, nfr );
bts_harm_iso = nan( npr, nza, naa, nfr );
bts_harm_cs1 = nan( npr, nza, naa, nfr );
bts_harm_cs2 = nan( npr, nza, naa, nfr );
emis_sea     = nan( npr, nza, naa, nfr );
emis_ice     = nan( npr, nza, naa, nfr );
emis_lan     = nan( npr, nza, naa, nfr );



%=== common Q definitions

Q.atmos_input.P	  = I.P;


%= Sensor

Q.sensor_input.F	  = frequencies;
Q.sensor_input.VPOL	  = v_polarization;
Q.sensor_input.HPOL	  = h_polarization;
Q.sensor_input.H	  = sensor_height;
Q.sensor_input.TOA        = atmosphere_adding;


%= Inputs

Q.dir_input             = [ E2E_HOME, '/SCEPSdata/InputData/ForwardModelData']; 

%=== Indicating variable for jacobians

Q.J_DO		  = 0;


%=== Indicating if mean emis or with added variability for ice

if ~do_rem
  Q.surf_input.EMI_ICE_COV	  = ~emis_ice_smoothing;
end

%=== Inputing an external comon land emis for all pixels if given
%    [ V-emis x no.freqs  H-emis x no.freqs ]

if emis_land_external ~= -999
  Q.surf_input.EMI_LAN		  = emis_land_external;
end


%=== Loading parameterized emissivity tables
%    OBSOLETE: using SURFEM-Ocean

afreq	= [1.410 6.925 10.650 18.700 36.500];
nf      = length( Q.sensor_input.F );


%= Sea emis

if 0

aangle	= [0 20 40 60 80 85 89];;  


for f = 1:nf

 for z = 1:nza

    [ ~, rfreq ]= nearest_in_vector( afreq, Q.sensor_input.F(f) );
    rfreq = floor(rfreq); 
    [ ~, rangle ]= nearest_in_vector( aangle, zenith_angle(f,z));

    edir	= [ Q.dir_input, '/EmisSea' ];
    efile	= [ edir, '/sea_emis_surfem_ocean_', sprintf('%02.0f', rfreq), 'GHz_', sprintf('%0.0f',rangle), 'deg.mat'];
    eval(['global ', 'sea_emis_surfem_ocean_', sprintf('%02.0f', rfreq), 'GHz_', sprintf('%0.0f',rangle) ]);
    m       = matfile(efile); 
    m = m.data;
    eval( ['sea_emis_surfem_ocean_', sprintf('%02.0f', rfreq), 'GHz_', sprintf('%0.0f',rangle), '=m;']);
    clear m 


  end

end

end


%= Ice emis

%  originally loading in memory but with parfor does not
%  work, so reading each time from matfile as with land emis

if 0

for f = 1:nf


 %= angle dependence not added yet 
 if f == 1
   aangle = [40 55 ];  
 else
   aangle = [55 ];   
 end 

 for z = 1:nza

    [ ~, rfreq ]= nearest_in_vector( afreq, Q.sensor_input.F(f) );
    rfreq = floor(rfreq); 
    [ ~, rangle ]= nearest_in_vector( aangle, zenith_angle(f,z));

    edir	= [ Q.dir_input, '/EmisIce' ];
    efile	= [ edir, '/ice_emis_hd_', sprintf('%02.0f', rfreq), 'GHz_', sprintf('%0.0f',rangle), 'deg.mat']
    eval(['global ', 'ice_emis_hd_', sprintf('%02.0f', rfreq), 'GHz_', sprintf('%0.0f',rangle) ]);
    m       = matfile(efile); 
    m = m.data;
    eval( ['ice_emis_hd_', sprintf('%02.0f', rfreq), 'GHz_', sprintf('%0.0f',rangle), '=m;'])
    clear m 

  end

end



end



%=== loading sea scatter table
%    NOT IN USE


if 0

global scatterm
edir	= [ Q.dir_input, '/ScatterSea/sea_scatter.mat'];
load(edir);
scatterm = reshape(data,91,50,26,13,2);
clear data
  
end

LOG.info( [ idfunction, ' ** Processing ', num2str(length(pro)), ' fields out of ', num2str(npr), ' for day ', SCENE_DATE ]);

%tic

per_pro = round(1:npr);



%=== parallelizing code with 4 cluster

% dummy sentences so parpool can identify
% them as variables

I   = I;
Q   = Q;
LOG = LOG;




%= Parallel mode not finished, so disable. The problem is that the global
%  variables cannot be handle in the parfor loop, so if the FM needs to reload
%  the  emis table for each call any advantages of the multiple processing are
%  lost due to the loading time of these files. Needs to be investigated.



if workers_number > 0 


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
 
    [ y, y_harm_iso, y_harm_cs1, y_harm_cs2, esea, eice, elan ]  = module_parallel( co, pro, jnd, zenith_angle, azimuth_angle, Q, I, LOG );
    

    bts(co,:,:,:)          = y;
    bts_harm_iso(co,:,:,:) = y_harm_iso;
    bts_harm_cs1(co,:,:,:) = y_harm_cs1;
    bts_harm_cs2(co,:,:,:) = y_harm_cs2;
    emis_sea(co,:,:,:)     = esea; 
    emis_ice(co,:,:,:)     = eice; 
    emis_lan(co,:,:,:)     = elan;

  end
  toc


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
    
    [ y, y_harm_iso, y_harm_cs1, y_harm_cs2, esea, eice, elan ]  = module_parallel( co, pro, jnd, zenith_angle, azimuth_angle, Q, I, LOG );
 
    bts(co,:,:,:)          = y;
    bts_harm_iso(co,:,:,:) = y_harm_iso;
    bts_harm_cs1(co,:,:,:) = y_harm_cs1;
    bts_harm_cs2(co,:,:,:) = y_harm_cs2;
    emis_sea(co,:,:,:)     = esea; 
    emis_ice(co,:,:,:)     = eice; 
    emis_lan(co,:,:,:)     = elan; 

  end
  toc

end




aux  = num2str(toc/60/60);
LOG.info( [ idfunction, ' ** Processing finished in ', aux , ' hours for day ', SCENE_DATE ]);



%=== Saving the sensor-free BTs
%    one file per frequency
%    in matlab internal format


bts      = single(bts);
emis_sea = single(emis_sea);
emis_ice = single(emis_ice);
emis_lan = single(emis_lan);



for f = 1:nb

  ind     = (4*(f-1)+1):(4*f);

  data    = bts(:,:,:,ind);
  sparam  = [ 'BTS_', bnames(f).txt];
  outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
  save_data_single( outfile, data );

  data    = bts_harm_iso(:,:,:,ind);
  sparam  = [ 'BTS_HARM_ISO_', bnames(f).txt];
  outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
  save_data_single( outfile, data );

  data    = bts_harm_cs1(:,:,:,ind);
  sparam  = [ 'BTS_HARM_CS1_', bnames(f).txt];
  outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
  save_data_single( outfile, data );

  data    = bts_harm_cs2(:,:,:,ind);
  sparam  = [ 'BTS_HARM_CS2_', bnames(f).txt];
  outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
  save_data_single( outfile, data );

  data    = emis_sea(:,:,:,ind);
  sparam  = [ 'EMIS_SEA_', bnames(f).txt];
  outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
  save_data_single( outfile, data );

  data    = emis_ice(:,:,:,ind);
  sparam  = [ 'EMIS_ICE_', bnames(f).txt ];
  outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
  save_data_single( outfile, data );

  data    = emis_lan(:,:,:,ind);
  sparam  = [ 'EMIS_LAN_', bnames(f).txt ];
  outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
  save_data_single( outfile, data  );


end


%=== Also saving lat, lon and PRO index and wind after data thinning

%= latitude
data    = single(I.LAT(pro));

sparam  = [ 'LAT'];
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save_data_single( outfile, data );


%= longitude

data    = single(I.LON(pro));
sparam  = [ 'LON'];
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save_data_single( outfile, data );


%= PRO index
data    = int32(spro);
sparam  = [ 'IND'];
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save_data_single( outfile, data );

%= MASK

data    = single(I.MASK(pro));
sparam  = [ 'MASK'];
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save_data_single( outfile, data );


%= UWS

data    = single(I.UWS(pro));
sparam  = [ 'UWS'];
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save_data_single( outfile, data );


%= VWS

data    = single(I.VWS(pro));
sparam  = [ 'VWS'];
outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
LOG.info( [ idfunction, ' ** Saving ', sparam, ' fields for day ', SCENE_DATE ]);
save_data_single( outfile, data );

end


%=== saving BTS in netcdf file for public distribution
%    and input to OSS


filesave  = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_BTS' ];  
SGM_NETCDF_SAVE_FILE = filesave;

LOG.info( [ idfunction, ' ** Saving netcdf file ', filesave ] );

comp(1).name = 'Vpo';
comp(2).name = 'Hpo';
comp(3).name = '3rd';
comp(4).name = '4th';

sta.date     = SCENE_DATE;



%= dimensions for reshaping

  outfile = [ GEOINPUT_SIMULATION, '_Output_', SCENE_DATE, '_nonvectorizedLAT', '.mat' ];
  LOG.info( [ idfunction, ' ** Loading nonvectorizedLAT fields for day ', SCENE_DATE ]);
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



%= files per azimuth angle

for az = 1:naa

  if length(unique(azimuth_angle(:,1))) > 1
    osfi_error('Azimuth angles are not constant for all frequencies so netcdf files cannot be saved')
  end

end


for z = 1:naa

  aa = azimuth_angle(1,z);
  ifilesave = [ filesave, '_aa_', sprintf( '%03i', aa ) ]; 
  LOG.info( [ idfunction, ' ** Saving ', ifilesave]);

  sparam  = [ 'LAT'];
  outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Saving ', sparam, ' netcdf fields for day ', SCENE_DATE, ' and azimuth angle ', sprintf( '%03i', aa ) ]);
  load( outfile );

  data = single(reshape(data,ro,ra));

  a = 1;
  str(a).value = data;
  str(a).netcdf_name = 'latitude';
  str(a).long_name   = 'Latitude';
  str(a).units       = 'degrees_north';
  clear data

  sparam  = [ 'LON'];
  outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Saving ', sparam, ' netcdf fields for day ', SCENE_DATE ]);
  load( outfile );

  data = single(reshape(data,ro,ra));

  a = a+1;
  str(a).value = data;
  str(a).netcdf_name = 'longitude';
  str(a).long_name   = 'Longitude';
  str(a).units       = 'degrees_east';
  clear data

  sparam  = [ 'MASK'];
  outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Loading ', sparam, ' netcdf fields for day ', SCENE_DATE ]);
  load( outfile );

  data = single(reshape(data,ro,ra));

  a = a+1;
  str(a).value = data;
  str(a).netcdf_name = 'mask';
  str(a).long_name   = 'Sea(1)-land(2)-ice(3) mask';
  str(a).units       = ' ';
  clear data

  sparam  = [ 'UWS'];
  outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Loading ', sparam, ' netcdf fields for day ', SCENE_DATE ]);
  load( outfile );

  data = single(reshape(data,ro,ra));

  a = a+1;
  str(a).value = data;
  str(a).netcdf_name = 'uws';
  str(a).long_name   = 'u wind speed';
  str(a).units       = 'm/s';
  clear data



  sparam  = [ 'VWS'];
  outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
  LOG.info( [ idfunction, ' ** Loading ', sparam, ' netcdf fields for day ', SCENE_DATE ]);
  load( outfile );

  data = single(reshape(data,ro,ra));

  a = a+1;
  str(a).value = data;
  str(a).netcdf_name = 'vws';
  str(a).long_name   = 'v wind speed';
  str(a).units       = 'm/s';
  clear data

  
  a = a+1;
  str(a).value = squeeze(zenith_angle(1,:));
  str(a).netcdf_name = 'incidence_angle';
  str(a).long_name   = 'Earth incidence angle';
  str(a).units       = 'degrees';


  for b = 1:nb
 
    sparam  = [ 'BTS_', bnames(b).txt];
    outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
    if b == 1
      LOG.info( [ idfunction, ' ** Loading BTS netcdf fields for day ', SCENE_DATE ]);
    end
    load( [ outfile ]);

    ng = size( data, 2 );
    nz = size( data, 3 );
    nt = size( data, 4 );

    for t = 1:nt

      dato = single( nan(ro, ra, ng ) );
      
      for g = 1:ng
        dato(:,:,g) = reshape(squeeze(data(:,g,z,t)),ro, ra);
      end

      ran = max(dato(:)) - min(dato(:));
      ofs = min(dato(:)) + ran/2;
      scf = 2 * 32767 / ran;
      dato = int16(scf * (dato - ofs));
    
      a = a+1;
      str(a).netcdf_name  = [ 'toa_tbs_', bnames(b).txt, '_', comp(t).name ];
      str(a).scale_factor  = 1/scf;
      str(a).add_offset    = ofs;
      str(a).value  = dato;

    end

  end  

  for b = 1:nb

    sparam  = [ 'BTS_HARM_ISO_', bnames(b).txt];
    outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
    if b == 1
      LOG.info( [ idfunction, ' ** Loading BTS_HARM_ISO netcdf fields for day ', SCENE_DATE ]);
    end
    load( [ outfile ]);

    ng = size( data, 2 );
    nz = size( data, 3 );
    nt = size( data, 4 );

    for t = 1:nt

      dato = single( nan(ro, ra, ng ) );
      % to faciliate rescaling, values to be identified with mask later

      for g = 1:ng
        dato(:,:,g) = reshape(squeeze(data(:,g,z,t)),ro, ra);
      end

      dato( dato == -999 ) = 0;      

      ran = max(dato(:)) - min(dato(:));
      ofs = min(dato(:)) + ran/2;
      scf = 2 * 32767 / ran;
      dato = int16(scf * (dato - ofs));
    
      a = a+1;
      str(a).netcdf_name  = [ 'toa_tbs_harmonic_iso_', bnames(b).txt, '_', comp(t).name ];
      str(a).scale_factor  = 1/scf;
      str(a).add_offset    = ofs;
      str(a).value  = dato;

    end

  end

  for b = 1:nb
 
    sparam  = [ 'BTS_HARM_CS1_', bnames(b).txt];
    outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
    if b == 1
      LOG.info( [ idfunction, ' ** Loading BTS_HARM_CS1 netcdf fields for day ', SCENE_DATE ]);
    end
    load( [ outfile ]);

    ng = size( data, 2 );
    nz = size( data, 3 );
    nt = size( data, 4 );

    for t = 1:nt

      dato = single( nan(ro, ra, ng ) );

      for g = 1:ng
        dato(:,:,g) = reshape(squeeze(data(:,g,z,t)),ro, ra);
      end

      % to faciliate rescaling, values to be identified with mask later
      dato( dato == -999 ) = 0;      

      ran = max(dato(:)) - min(dato(:));
      ofs = min(dato(:)) + ran/2;
      scf = 2 * 32767 / ran;
      dato = int16(scf * (dato - ofs));
    
      a = a+1;
      str(a).netcdf_name  = [ 'toa_tbs_harmonic_1_', bnames(b).txt, '_', comp(t).name ];
      str(a).scale_factor  = 1/scf;
      str(a).add_offset    = ofs;
      str(a).value  = dato;

    end

  end


  for b = 1:nb
 
    sparam  = [ 'BTS_HARM_CS2_', bnames(b).txt];
    outfile = [ dirout, '/', idfunction, '_Output_', SCENE_DATE, '_', sparam, '.mat' ];  
    if b == 1
       LOG.info( [ idfunction, ' ** Loading BTS_HARM_CS2 netcdf fields for day ', SCENE_DATE ]);
    end
    load( [ outfile ]);

    ng = size( data, 2 );
    nz = size( data, 3 );
    nt = size( data, 4 );

    for t = 1:nt

      dato = single( nan(ro, ra, ng ) );

      for g = 1:ng
        dato(:,:,g) = reshape(squeeze(data(:,g,z,t)),ro, ra);
      end

      % to faciliate rescaling, values to be identified with mask later
      dato( dato == -999 ) = 0;      

      ran = max(dato(:)) - min(dato(:));
      ofs = min(dato(:)) + ran/2;
      scf = 2 * 32767 / ran;
      dato = int16(scf * (dato - ofs));
    
      a = a+1;
      str(a).netcdf_name  = [ 'toa_tbs_harmonic_2_', bnames(b).txt, '_', comp(t).name ];
      str(a).scale_factor  = 1/scf;
      str(a).add_offset    = ofs;
      str(a).value  = dato;

    end

  end  
  
  dat_aux_write_geo_card_tbs_netcdf( str, sta, [ ifilesave, '.nc' ] );

end



return


%-------------------------------------


function dat_aux_write_geo_card_tbs_netcdf( str, sta, savename )


%=== setting no compression

dflevel  = 1; 


%=== removing file if exists

if exist( savename, 'file' )
  delete( savename );
end



%=== checking sizes


if strcmp(str(1).netcdf_name,'latitude')
  nj = size(str(1).value,1);
  ni = size(str(1).value,2);
else
  error('Latitude should be the first cell in str');
end
if ~strcmp(str(2).netcdf_name,'longitude')
  error('Longitude should be the second cell in str');
end



%=== time dimension
nt  = 1;


%=== writing variables
nf = length(str);



%= incindence angle
nz = length( str(4).value);

for f = 1:nf


    itype = class( str(f).value );

    if f < 6

      nccreate( savename, str(f).netcdf_name, 'Dimensions', { 'lon' nj 'lat' ni 'time' nt}, 'DataType', itype,  'DeflateLevel', dflevel, 'Format', 'netcdf4_classic');

    elseif f == 6

      nccreate( savename, str(f).netcdf_name, 'Dimensions', { 'incidence_angle' nz }, 'DataType', itype,  'DeflateLevel', dflevel, 'Format', 'netcdf4_classic'    );

    else

      nccreate( savename, str(f).netcdf_name, 'Dimensions', { 'lon' nj 'lat' ni 'incidence_angle' nz 'time' nt}, 'DataType', itype,  'DeflateLevel', dflevel, 'Format', 'netcdf4_classic'    );

    end

    ncwrite( savename, str(f).netcdf_name, str(f).value );

    if ~isempty( str(f).add_offset )
      ncwriteatt( savename, str(f).netcdf_name, 'add_offset', str(f).add_offset );
    end

    if ~isempty( str(f).scale_factor )
      ncwriteatt( savename, str(f).netcdf_name, 'scale_factor', str(f).scale_factor );
    end
 
    ncwriteatt( savename, str(f).netcdf_name, 'long_name', str(f).long_name );

    ncwriteatt( savename, str(f).netcdf_name, 'standard_name', str(f).netcdf_name );

    ncwriteatt( savename, str(f).netcdf_name, 'units', str(f).units );


end



%=== writing global attributes

ncwriteatt( savename, '/', 'title', [ 'Top-of-atmosphere brigthness temperatures for geophysical scene, ', sta.date ] );
ncwriteatt( savename, '/', 'project', 'Copernicus Imaging Microwave Radiometer - SCIENTIFIC END-TO-END PERFORMANCE SIMULATION (CIMR-SCEPS) ;' );
ncwriteatt( savename, '/', 'summary', 'CIMR-SCEPS Work Package 2.1'  );
ncwriteatt( savename, '/', 'description', 'Test data set that represents a reference scenario for the evaluation of CIMR mission'  );
ncwriteatt( savename, '/', 'work_package_manager_name', 'Carlos Jimenez' );
ncwriteatt( savename, '/', 'work_package_manager_email', 'carlos.jimenez@estellus.fr' );
ncwriteatt( savename, '/', 'source', 'CIMR-SCEPS Test Data Sets' );
ncwriteatt( savename, '/', 'spatial_resolution', '1km EASE grid' );
ncwriteatt( savename, '/', 'acknowledgement' , 'Data providers are acknowledged by making their data available' ); 
ncwriteatt( savename, '/', 'creator_name', 'Carlos Jimenez, Estellus' );
ncwriteatt( savename, '/', 'creator_email', 'carlos.jimenez@estellus.fr' );
ncwriteatt( savename, '/', 'creator_url', 'www.estellus.fr' );
ncwriteatt( savename, '/', 'history', [ 'Created ', date, ' 8.30'] );


return










function [ bts, bts_harm_iso, bts_harm_cs1, bts_harm_cs2, esea, eice, elan ]  = module_parallel( co, pro, dopro, zenith_angle, azimuth_angle, Q, I, LOG );





%= Ice emis
%
%  originally loading once outside this script, but for
%  the parfor to works, the global variables need to be defined here
%  it requires the same memory x worker numbers, but still
%  not that large and 2x faster than loading every time
%  in forward_model_core

afreq	= [1.410 6.925 10.650 18.700 36.500];
nf      = length( Q.sensor_input.F );
nza     = size(zenith_angle,2);


%= hard coded switch before confirmin new SCEPS sea-ice emis codes

do_nnemis = 0;

if ~do_nnemis

for f = 1:nf


 %= angle dependence not added yet 
 if f == 1
   aangle = [40 55 ];  
 else
   aangle = [55 ];   
 end 

 for z = 1:nza

    [ ~, rfreq ]= nearest_in_vector( afreq, Q.sensor_input.F(f) );
    rfreq = floor(rfreq); 
    [ ~, rangle ]= nearest_in_vector( aangle, zenith_angle(f,z));

    edir	= [ Q.dir_input, '/EmisIce' ];
    efile	= [ edir, '/ice_emis_hd_', sprintf('%02.0f', rfreq), 'GHz_', sprintf('%0.0f',rangle), 'deg.mat'];
    eval(['global ', 'ice_emis_hd_', sprintf('%02.0f', rfreq), 'GHz_', sprintf('%0.0f',rangle) ]);
    m       = matfile(efile); 
    m = m.data;
    eval( ['ice_emis_hd_', sprintf('%02.0f', rfreq), 'GHz_', sprintf('%0.0f',rangle), '=m;'])
    clear m 

  end

end

end



p = pro(co);

nza    = size(zenith_angle,2);
naa    = size(azimuth_angle,2);
nfr    = sum(Q.sensor_input.VPOL) + sum(Q.sensor_input.HPOL)+2*length(Q.sensor_input.F);

%= initializing

bts    = nan( nza, naa, nfr );
esea   = nan( nza, naa, nfr );
eice   = nan( nza, naa, nfr );
elan   = nan( nza, naa, nfr );


bts_harm_iso = nan( nza, naa, nfr );
bts_harm_cs1 = nan( nza, naa, nfr );
bts_harm_cs2 = nan( nza, naa, nfr );


if sum(dopro==p) >0 

  clear surpre

  %=  SST is just a dummy choice as there are no inversions
  % here. But is is required to parse the forward model
  R.R_VARIABLES		  = {'SST'};

  %= atmosphere

  Q.atmos_input.P	  = I.P';
  Q.atmos_input.T	  = I.T(p,:);
  Q.atmos_input.RHO	  = I.RHO(p,:);
  Q.atmos_input.LW	  = I.LW(p,:);
  % No more checking and rescaling of colum content here
  % See previous note  
  Q.atmos_input.CWVC	  = I.CWVC(p);
  %Q.atmos_input.CLWC	  = I.CLWC(p);


  %= sea surface

  mask = I.MASK;

  if mask(p) == 1 | ( mask(p) == 2 & I.SIC(p)> 0 & I.SIC(p) < 1)
   
    % sea; for ice they are needed
    % as if 0 < SIC < 1 the scenes is ice-sea
    % mixed

    Q.surf_input.SSS	    = I.SSS(p);	
    Q.surf_input.SLP        = I.SLP(p);
    Q.surf_input.OWS        = I.OWS(p);    
    Q.surf_input.UWS        = I.UWS(p);    
    Q.surf_input.VWS        = I.VWS(p);    
    Q.surf_input.SST        = I.SST(p);    
    surpre = Q.surf_input.SLP;

  end

  %= ice surface

  if mask(p) == 2

    % ice surface, sea+ice, land+ice 

    Q.surf_input.IST      = I.IST(p);        
    Q.surf_input.T2M      = I.T2M(p);            
    Q.surf_input.SIC      = I.SIC(p);            		
    Q.surf_input.SND      = I.SND(p);      
    Q.surf_input.FDN      = I.FDN(p);          
    Q.surf_input.SIT      = I.SIT(p);          
    Q.surf_input.XYI      = I.XYI(p);            
    Q.surf_input.ILP      = I.ILP(p);            
    Q.surf_input.LAT      = I.LAT(p);        
    Q.surf_input.LON      = I.LON(p);            
    Q.surf_input.MONTH    = I.MONTH(p);        
    surpre = Q.surf_input.ILP;    

  end   

  %= land surface

  if mask(p) >= 3

    %  land including coastal pixels 

    Q.surf_input.LST      = I.LST(p);        
    Q.surf_input.LAT      = I.LAT(p);        
    Q.surf_input.LON      = I.LON(p);        
    Q.surf_input.MONTH    = I.MONTH(p);        
    Q.surf_input.LLP      = I.LLP(p);        
    surpre = Q.surf_input.LLP;

  end


  if isfield( Q.surf_input, 'LST' )

    %= Land emis
    %   Reading and interpolating in angle climatological emis
    %   here and passing them as Q.surf_input.EMI_LAN if not
    %   given externally

    EMI_LAN = emis_land_angle_interpolate( Q.surf_input.LON, Q.surf_input.LAT, Q.surf_input.MONTH, Q.sensor_input.F, zenith_angle, [], Q.dir_input);


  end


  if isfield( Q.surf_input, 'IST' )

    if 0

    t2m  = 275:-1:240;
    semi = nan(length(t2m),2,5,8);

    for tt = 1:length(t2m)

      Q.surf_input.T2M = t2m(tt);

      EMI_ICE = emis_ice_angle_interpolate( Q.surf_input.LAT, Q.surf_input.LON, Q.surf_input.MONTH, Q.surf_input.XYI, Q.surf_input.T2M, Q.surf_input.SND, Q.surf_input.SIT, Q.atmos_input.CWVC  , Q.surf_input.IST, Q.surf_input.FDN, Q.sensor_input.F, zenith_angle, [], Q.dir_input, do_nnemis);

      semi(tt,:,:,:) = EMI_ICE; 

    end

    end

    EMI_ICE = emis_ice_angle_interpolate( Q.surf_input.LAT, Q.surf_input.LON, Q.surf_input.MONTH, Q.surf_input.XYI, Q.surf_input.T2M, Q.surf_input.SND, Q.surf_input.SIT, Q.atmos_input.CWVC  , Q.surf_input.IST, Q.surf_input.FDN, Q.sensor_input.F, zenith_angle, [], Q.dir_input, do_nnemis);

  end

  if isempty(find(isnan(Q.atmos_input.T)==1))
  
    try 
  
      for za = 1:nza

        if size(zenith_angle,1) == 1
          Q.sensor_input.ZA = ones(length(Q.sensor_input.F),1) * zenith_angle(za); 
        else
          Q.sensor_input.ZA = zenith_angle(:,za); 
        end

        if exist( 'EMI_LAN' )
          Q.surf_input.EMI_LAN = EMI_LAN(:,:,za);
        end

        if exist( 'EMI_ICE' )
          Q.surf_input.EMI_ICE = EMI_ICE(:,:,za);
        end

   
        for aa = 1:naa
  
          if size(azimuth_angle,1) == 1
            Q.sensor_input.AA = ones(length(Q.sensor_input.F),1) * azimuth_angle(aa); 
          else
            Q.sensor_input.AA = azimuth_angle(:,aa); 
          end


          %= adjusting atmopsheric profiles to make them consistent with 
	  %  surface pressure. If surpre > p(1) we just extrapolate to
          %  the P(1) atmospheric quantity value. Otherwise, we linearly
	  %  interpolate with log10(P)

          Q.atmos_input = forward_model_aux_slp( Q.atmos_input, surpre );

          %= forward modelling


          [ R, bts(za,aa,:), bts_harm_iso(za,aa,:), bts_harm_cs1(za,aa,:), bts_harm_cs2(za,aa,:)] = forward_model( Q, R, [], 1 );

          %= the harmonics components are only calculaed for L-band at 50deg
          %  and reamining bands at 55deg. And they do not depend on aa
          %  the sin-cos dependence in aa is added later in simulator. So the
          %  values for all aas are always the sameand we only store them
          %  once uding the correponding switches

 
          esea(za,aa,:)   = R.emis.sea; 
          eice(za,aa,:)   = R.emis.ice; 
          elan(za,aa,:)   = R.emis.lan; 

        end

      end  


    catch me

      LOG.info( sprintf(' CELL %0.2f DID NOT PROCESS \n', p ));
      LOG.info( [ 'ERROR: ', me.message]);

    end

  else

    LOG.info( sprintf(' CELL %0.2f DID NOT PROCESS \n', p ));

  end


%else

%  aux  = nan(1,sum(Q.sensor_input.VPOL)+sum(Q.sensor_input.HPOL));
%  bts  = aux;
%  esea = aux;
%  eice = aux;
%  elan = aux;

end


return




%==========================================================================


function emis = emis_ice_nns( sit, ist, t2m, snd, tcwv, fdn, dir_input, freq)

edir	= [ dir_input, '/EmisIceNNs/NNs_seaice_500epochs' ];

switch freq

  case 1
    efile	= [ edir, '/net_seaice_propre_allyear_sic95_inputs_124567_smap_vh_trainbr_nn60_epoch500.mat'];
  case 6
    efile	= [ edir, '/net_seaice_propre_allyear_sic95_inputs_124567_amsr2_ind1_vh_trainbr_nn60_epoch500.mat'];
  case 10
    efile	= [ edir, '/net_seaice_propre_allyear_sic95_inputs_124567_amsr2_ind2_vh_trainbr_nn60_epoch500.mat'];
  case 18
    efile	= [ edir, '/net_seaice_propre_allyear_sic95_inputs_124567_amsr2_ind3_vh_trainbr_nn60_epoch500.mat'];
  case 36
    efile	= [ edir, '/net_seaice_propre_allyear_sic95_inputs_124567_amsr2_ind4_vh_trainbr_nn60_epoch500.mat'];

  end

load( efile, 'net' );

inputs  = [snd sit ist t2m tcwv fdn ];
emis    = net(inputs');



return




%==========================================================================


function emis  = emis_ice_surface_table( lat, xyi, t2m, snd, dir_input, freq)


%= table organized as T2M, SND, XYI, LAT
%it2m = 240:0.1:270;
%isnd = [0:0.01:0.6];
%ilat = [-45 45]; 
%ixyi = [1 -1];   %FYI MYI




angle = 55;


%= loading data from global variable


edir	= [ dir_input, '/EmisIce' ];
efile	= [ edir, '/ice_emis_hd_smooth_', sprintf('%02.0f', freq), 'GHz_55deg.mat'];
load(efile,'data');

it = 1+ round(10*(t2m-240));
is = 1+ round(100*snd); 

%= is limited to snd=0.6 m 
%  with the RRDP
if is>61
  is=61;
end

%= it limited to 270, same reasons
if it>301
  it =301;
end
if it<1
  it =1;
end

if lat>0
  il = 2;
else
  il = 1;
end

% XYI = 0 is FYI, XYI = 1 is MYI
% values between 0 and 1 are to smooth emis
% at edges

% data(t,s,x,l,:)


if xyi ==1

  % MYI
  emis  = squeeze(double(data( it, is, 2, il, : ))/1e3);
  if emis(1) == 0 | emis(2) == 0
    emis(1) = -999;
    emis(2) = -999;
  end

elseif xyi == 0

  % FYI
  emis  = squeeze(double(data( it, is, 1, il, : ))/1e3);
  if emis(1) == 0 | emis(2) == 0
    emis(1) = -999;
    emis(2) = -999;
  end

else

  emis0  = squeeze(double(data( it, is, :, il, 1 ))/1e3);
  emis1  = squeeze(double(data( it, is, :, il, 2 ))/1e3);

  if emis0(1) == 0 | emis0(2) == 0 | emis1(1) == 0 | emis1(2) == 0

    emis(1) = -999;
    emis(2) = -999;

  else

    emis   = [emis0(1) * (1-xyi) + emis0(2) * xyi; emis1(1) * (1-xyi) + emis1(2) * xyi];

  end

end



return




%==========================================================================

function emis  = emis_land_angle_interpolate( lon, lat, month, freqs, angles, emis, dir_input )


nf = length( freqs );
na = size( angles,2 );

%== Ancillary coefficients, for the 18 et 36 frequencies and the 10 
%   surface classes (indicated by class1 in the TELSEM2 atlases): 

a0_k0= [ 0.11509  0.091535 ;...
         0.10525  0.16627  ;...
         0.29217  0.23809  ;...
         0.17516  0.19459  ;...
         0.10521  0.12126  ;...
         0.18212  0.19625  ;...
        -0.19202  0.5411   ;...
         0.10292  0.5486   ;...
        -0.022672 0.44492  ;...
        -0.33894 -0.17621  ];
a0_k1= [ 0.61168 0.59095 ;...
         0.60271 0.69213 ;...
         0.32728 0.34334 ;...
         0.51217 0.4491  ;...
         0.48913 0.41932 ;...
         0.64474 0.30637 ;...   
         1.0405  0.17538 ;...
         0.61819 0.31298 ;... 
         0.87761 0.47583 ;...
         1.0959 0.92842 ];
a0_k2= [ 0.26726 0.32033 ;...
         0.28547 0.13592 ;...
         0.37178 0.41813 ;...
         0.30203 0.35479 ;...   
         0.40663 0.47493 ;...
         0.14811 0.52382 ;...   
         0.14286 0.27164 ;...
         0.2737 0.12001;...   
         0.13492 0.065463 ;...
         0.24905 0.25475 ];
a0_eveh=[0.9592599869E+00 0.9565299749E+00 ;...   
         0.9560700059E+00 0.9541199803E+00 ;...   
         0.9461100101E+00 0.9439799786E+00 ;...   
         0.9317600131E+00 0.9289000034E+00 ;... 
         0.9208700061E+00 0.9190599918E+00 ;...   
         0.9162799716E+00 0.8937299848E+00 ;...   
         0.9570500255E+00 0.9213600159E+00 ;...   
         0.9639400244E+00 0.9530599713E+00 ;...   
         0.9685299993E+00 0.9622600079E+00 ;...   
         0.8997200131E+00 0.9012699723E+00 ];
a1_eveh=[0.3627802414E-07 -0.7778328204E-08 ;...   
         0.2503205394E-06 0.1996262995E-06 ;...   
         0.4190530660E-06 0.3655744649E-06 ;...   
         0.5574374313E-06 0.5273076340E-06 ;...   
         0.1026844529E-05 0.9679998811E-06 ;...   
         0.3180800832E-06 0.2886778532E-06 ;...   
        -0.1118036366E-06 -0.1502856577E-06 ;...   
        -0.8410978580E-08 -0.3478669441E-07 ;...   
         0.2485776633E-06 0.1800235907E-06 ;...   
	  0.2687000915E-06 0.1740325644E-06 ];
a2_eveh=[0.3067140824E-05 0.2520012231E-05 ;...   
         0.8213598448E-05 0.7378375358E-05 ;...   
         0.1225889173E-04 0.1165553113E-04 ;...   
         0.1693615741E-04 0.1648317448E-04 ;...   
         0.2744720041E-04 0.2642072104E-04 ;...   
         0.1349592094E-04 0.1261523357E-04 ;...   
         0.2064244654E-05 0.1919016057E-06 ;...   
         0.5334760772E-05 0.4130339221E-05 ;...   
         0.6530796327E-05 0.5727014013E-05 ;...   
         0.1071246970E-04 0.9539280654E-05 ];
a3_eveh=[-0.2004991551E-07 -0.6895366056E-07 ;...   
         -0.7322448425E-07 -0.1273002681E-06 ;...   
         -0.9421125213E-07 -0.1683332300E-06 ;...   
         -0.1317753799E-06 -0.2107972250E-06 ;...   
         -0.1889465580E-06 -0.2757958271E-06 ;...   
          0.7339644004E-08 -0.4058669560E-06 ;...   
          0.6170279931E-07 -0.1998567996E-06 ;...   
         -0.1361754887E-07 -0.1765622955E-06 ;...   
         -0.3901189061E-07 -0.1305666189E-06 ;...   
         -0.2679148992E-07 -0.4441960044E-07 ];
b0_eveh=[ 0.9592599869E+00 0.9565299749E+00 ;...   
          0.9560700059E+00 0.9541199803E+00 ;...   
          0.9461100101E+00 0.9439799786E+00 ;...   
          0.9317600131E+00 0.9289000034E+00 ;...   
          0.9208700061E+00 0.9190599918E+00 ;...   
          0.9162799716E+00 0.8937299848E+00 ;...   
          0.9570500255E+00 0.9213600159E+00 ;...   
          0.9639400244E+00 0.9530599713E+00 ;...   
          0.9685299993E+00 0.9622600079E+00 ;...   
          0.8997200131E+00 0.9012699723E+00 ];
b1_eveh=[ 0.3626608347E-07 -0.7786279177E-08 ;...   
          0.2502746099E-06 0.1995944388E-06 ;...   
          0.4189516289E-06 0.3655020180E-06 ;...   
          0.5572838404E-06 0.5271903092E-06 ;...   
          0.1026605219E-05 0.9677979733E-06 ;...   
          0.3179358714E-06 0.2884899004E-06 ;...   
         -0.1118781370E-06 -0.1503948681E-06 ;...   
         -0.8455684153E-08 -0.3485171618E-07 ;...   
          0.2485595019E-06 0.1799959364E-06 ;...   
          0.2686167306E-06 0.1739760478E-06 ];
b2_eveh=[ 0.3065537157E-05 0.2518960400E-05 ;...   
          0.8209894986E-05 0.7375769655E-05 ;...   
          0.1225203869E-04 0.1165053800E-04 ;...   
          0.1692612022E-04 0.1647546378E-04 ;...   
          0.2743142431E-04 0.2640772436E-04 ;...   
          0.1348545720E-04 0.1260529825E-04 ;...   
          0.2058213340E-05 0.1860650656E-06 ;...   
          0.5330772183E-05 0.4126528893E-05 ;...   
          0.6528573977E-05 0.5725009032E-05 ;...   
          0.1070590315E-04 0.9534271157E-05 ];
b3_eveh=[-0.1370247134E-06 -0.1436897747E-06 ;...   
         -0.3118435643E-06 -0.2916583242E-06 ;...   
         -0.5048401022E-06 -0.4662823869E-06 ;...   
         -0.7210980471E-06 -0.6662896794E-06 ;...   
         -0.1110204039E-05 -0.1030801400E-05 ;...   
         -0.6330818110E-06 -0.9186441048E-06 ;...   
         -0.3242539890E-06 -0.5027602583E-06 ;...   
         -0.2747250676E-06 -0.3811997260E-06 ;...   
         -0.1994112324E-06 -0.2555484855E-06 ;...   
         -0.4413041665E-06 -0.3717419474E-06 ];


%=  emissivity at 53 deg

if isempty( emis )

  e_clim  = nan(nf,2);

  %= land emis per frequency  
  [e_clim(1,:), class1 ] = emis_land_surface_table( lon, lat, month, dir_input, floor(freqs(1))); 
  for f = 2:nf
    e_clim(f,:) = emis_land_surface_table( lon, lat, month, dir_input, floor(freqs(f))); 
  end 

end


%= better have a made up emis than a -999
%  value for the few cases where no emis
%  lan exists

e_clim( e_clim == -999 ) = 0.9;


%=== Interpolation in angle for 18 and 36

ifre = find( floor(freqs) == 18 | floor(freqs) == 36 );


emiss_interp_v=nan(na,2);
emiss_interp_h=nan(na,2);

theta0  = 0;
theta53 = 53;

for f = 1:2

  indf = ifre(f); 
  ev53 = e_clim(ifre(f),1);
  eh53 = e_clim(ifre(f),2);
            
  for a = 1:na
    
    theta  = angles(ifre(f),a);

    %=  Calculation of the emis at theta=0� with a multilinear regression
    e0 = a0_k0(class1,f)+a0_k1(class1,f)*ev53+a0_k2(class1,f)*eh53;

    %=  Reading of the polynomial coefficients for ev et eh
    a0 = a0_eveh(class1,f);
    a1 = a1_eveh(class1,f);
    a2 = a2_eveh(class1,f);
    a3 = a3_eveh(class1,f);
    b0 = b0_eveh(class1,f);
    b1 = b1_eveh(class1,f);
    b2 = b2_eveh(class1,f);
    b3 = b3_eveh(class1,f);

    %= Vertical polarization
    S1_v      = ((theta-theta53)/(theta0-theta53)) * ((e0-a0)/a0);
    em53_v    = a3*(theta53^3) + a2*(theta53^2) + a1*theta53 + a0;
    S2_v      =((theta-theta0)/(theta53-theta0))*((ev53-em53_v)/em53_v);
    S_v       = 1 + S1_v + S2_v;
    emtheta_v = a3*(theta^3) + a2*(theta^2) + a1*theta + a0;
    emiss_interp_v(a,f) = S_v * emtheta_v;     

    %= Horizontal polarization 
    S1_h      = ((theta-theta53)/(theta0-theta53)) * ((e0-b0)/b0);
    em53_h    = b3*(theta53^3) + b2*(theta53^2) + b1*theta53 + b0;
    S2_h      =((theta-theta0)/(theta53-theta0))*((eh53-em53_h)/em53_h);
    S_h       = 1 + S1_h + S2_h;
    emtheta_h = b3*(theta^3) + b2*(theta^2) + b1*theta + b0;
    emiss_interp_h(a,f) = S_h * emtheta_h; 

  end

end



%= Emissivity V has to be larger or equal to emissivity H


for f = 1:2

  for a = 1:na
    
    if (emiss_interp_v(a,f) < emiss_interp_h(a,f))
       emiss_interp_v(a,f) = (emiss_interp_v(a,f) + emiss_interp_h(a,f))/2.;
       emiss_interp_h(a,f) =  emiss_interp_v(a,f);
    end

  end

end



%= for freqs <= 18GHz using angular dependence of 18GHz
%      freqs >= 18GHz using angular dependence of 36GHz

np   = 2; % v and h
emis = nan( np,nf, na);


for a = 1:na

  for f = 1:nf

     if floor(freqs(f)) <= 18
       s = 1;
     else
       s = 2;
     end

     emis(1,f,a) = e_clim(f,1) * emiss_interp_v(a,s) /e_clim(ifre(s),1);
     emis(2,f,a) = e_clim(f,2) * emiss_interp_h(a,s) /e_clim(ifre(s),2);

  end

end



%=  Emissivity cannot be larger than 1 (this is not in TELSEM2, but for
%   large angles this can happen here with antenna integration, 
%    limits have to be set


if 0
ind = find(emiss_interp_v > 1);
emiss_interp_v(ind) = 1;

ind = find(emiss_interp_h > 1);
emiss_interp_h(ind) = 1;
end





return






%==========================================================================




function emis  = emis_ice_angle_interpolate( lat, lon, month, xyi, t2m, snd, sit, tcwv, ist, fdn, freqs, angles, emis, dir_input, do_nnemis )


nf = length( freqs );
na = size( angles,2 );

%== Ancillary coefficients, for the 18 et 36 frequencies and the 10 
%   surface classes (indicated by class1 in the TELSEM2 atlases): 

a0_k0= [ 0.11509  0.091535 ;...
         0.10525  0.16627  ;...
         0.29217  0.23809  ;...
         0.17516  0.19459  ;...
         0.10521  0.12126  ;...
         0.18212  0.19625  ;...
        -0.19202  0.5411   ;...
         0.10292  0.5486   ;...
        -0.022672 0.44492  ;...
        -0.33894 -0.17621  ];
a0_k1= [ 0.61168 0.59095 ;...
         0.60271 0.69213 ;...
         0.32728 0.34334 ;...
         0.51217 0.4491  ;...
         0.48913 0.41932 ;...
         0.64474 0.30637 ;...   
         1.0405  0.17538 ;...
         0.61819 0.31298 ;... 
         0.87761 0.47583 ;...
         1.0959 0.92842 ];
a0_k2= [ 0.26726 0.32033 ;...
         0.28547 0.13592 ;...
         0.37178 0.41813 ;...
         0.30203 0.35479 ;...   
         0.40663 0.47493 ;...
         0.14811 0.52382 ;...   
         0.14286 0.27164 ;...
         0.2737 0.12001;...   
         0.13492 0.065463 ;...
         0.24905 0.25475 ];
a0_eveh=[0.9592599869E+00 0.9565299749E+00 ;...   
         0.9560700059E+00 0.9541199803E+00 ;...   
         0.9461100101E+00 0.9439799786E+00 ;...   
         0.9317600131E+00 0.9289000034E+00 ;... 
         0.9208700061E+00 0.9190599918E+00 ;...   
         0.9162799716E+00 0.8937299848E+00 ;...   
         0.9570500255E+00 0.9213600159E+00 ;...   
         0.9639400244E+00 0.9530599713E+00 ;...   
         0.9685299993E+00 0.9622600079E+00 ;...   
         0.8997200131E+00 0.9012699723E+00 ];
a1_eveh=[0.3627802414E-07 -0.7778328204E-08 ;...   
         0.2503205394E-06 0.1996262995E-06 ;...   
         0.4190530660E-06 0.3655744649E-06 ;...   
         0.5574374313E-06 0.5273076340E-06 ;...   
         0.1026844529E-05 0.9679998811E-06 ;...   
         0.3180800832E-06 0.2886778532E-06 ;...   
        -0.1118036366E-06 -0.1502856577E-06 ;...   
        -0.8410978580E-08 -0.3478669441E-07 ;...   
         0.2485776633E-06 0.1800235907E-06 ;...   
	  0.2687000915E-06 0.1740325644E-06 ];
a2_eveh=[0.3067140824E-05 0.2520012231E-05 ;...   
         0.8213598448E-05 0.7378375358E-05 ;...   
         0.1225889173E-04 0.1165553113E-04 ;...   
         0.1693615741E-04 0.1648317448E-04 ;...   
         0.2744720041E-04 0.2642072104E-04 ;...   
         0.1349592094E-04 0.1261523357E-04 ;...   
         0.2064244654E-05 0.1919016057E-06 ;...   
         0.5334760772E-05 0.4130339221E-05 ;...   
         0.6530796327E-05 0.5727014013E-05 ;...   
         0.1071246970E-04 0.9539280654E-05 ];
a3_eveh=[-0.2004991551E-07 -0.6895366056E-07 ;...   
         -0.7322448425E-07 -0.1273002681E-06 ;...   
         -0.9421125213E-07 -0.1683332300E-06 ;...   
         -0.1317753799E-06 -0.2107972250E-06 ;...   
         -0.1889465580E-06 -0.2757958271E-06 ;...   
          0.7339644004E-08 -0.4058669560E-06 ;...   
          0.6170279931E-07 -0.1998567996E-06 ;...   
         -0.1361754887E-07 -0.1765622955E-06 ;...   
         -0.3901189061E-07 -0.1305666189E-06 ;...   
         -0.2679148992E-07 -0.4441960044E-07 ];
b0_eveh=[ 0.9592599869E+00 0.9565299749E+00 ;...   
          0.9560700059E+00 0.9541199803E+00 ;...   
          0.9461100101E+00 0.9439799786E+00 ;...   
          0.9317600131E+00 0.9289000034E+00 ;...   
          0.9208700061E+00 0.9190599918E+00 ;...   
          0.9162799716E+00 0.8937299848E+00 ;...   
          0.9570500255E+00 0.9213600159E+00 ;...   
          0.9639400244E+00 0.9530599713E+00 ;...   
          0.9685299993E+00 0.9622600079E+00 ;...   
          0.8997200131E+00 0.9012699723E+00 ];
b1_eveh=[ 0.3626608347E-07 -0.7786279177E-08 ;...   
          0.2502746099E-06 0.1995944388E-06 ;...   
          0.4189516289E-06 0.3655020180E-06 ;...   
          0.5572838404E-06 0.5271903092E-06 ;...   
          0.1026605219E-05 0.9677979733E-06 ;...   
          0.3179358714E-06 0.2884899004E-06 ;...   
         -0.1118781370E-06 -0.1503948681E-06 ;...   
         -0.8455684153E-08 -0.3485171618E-07 ;...   
          0.2485595019E-06 0.1799959364E-06 ;...   
          0.2686167306E-06 0.1739760478E-06 ];
b2_eveh=[ 0.3065537157E-05 0.2518960400E-05 ;...   
          0.8209894986E-05 0.7375769655E-05 ;...   
          0.1225203869E-04 0.1165053800E-04 ;...   
          0.1692612022E-04 0.1647546378E-04 ;...   
          0.2743142431E-04 0.2640772436E-04 ;...   
          0.1348545720E-04 0.1260529825E-04 ;...   
          0.2058213340E-05 0.1860650656E-06 ;...   
          0.5330772183E-05 0.4126528893E-05 ;...   
          0.6528573977E-05 0.5725009032E-05 ;...   
          0.1070590315E-04 0.9534271157E-05 ];
b3_eveh=[-0.1370247134E-06 -0.1436897747E-06 ;...   
         -0.3118435643E-06 -0.2916583242E-06 ;...   
         -0.5048401022E-06 -0.4662823869E-06 ;...   
         -0.7210980471E-06 -0.6662896794E-06 ;...   
         -0.1110204039E-05 -0.1030801400E-05 ;...   
         -0.6330818110E-06 -0.9186441048E-06 ;...   
         -0.3242539890E-06 -0.5027602583E-06 ;...   
         -0.2747250676E-06 -0.3811997260E-06 ;...   
         -0.1994112324E-06 -0.2555484855E-06 ;...   
         -0.4413041665E-06 -0.3717419474E-06 ];


%=  emissivity at 53 deg


if isempty( emis )

  e_clim  = nan(nf,2);

  if ~do_nnemis
    for f = 1:nf
      e_clim(f,:) = emis_ice_surface_table( lat, xyi, t2m, snd, dir_input, floor(freqs(f))); 
    end
  else
    for f = 1:nf
      e_clim(f,:)  = emis_ice_nns( sit, ist, t2m, snd, tcwv, fdn, dir_input, floor(freqs(f))); 
    end
  end

end


%=== class for interpolation

if 1

efile	= [ dir_input, '/EmisLand/Res_004/land_emis_hd_class_', sprintf('%02.0f',month), '.mat'];
m       = matfile(efile); 

flon    = 0.018:0.036:359.982;
flat    = -89.982:0.036:89.982;

if lon < 0
  lon = lon + 360;
end

[ ~, ilo ] = min( abs( flon - lon ) );
[ ~, ila ] = min( abs( flat - lat ) );


class1  = double(m.data( ilo, ila, 1));
class2  = double(m.data( ilo, ila, 2));
class1( class2 == 12 | class2 == 13 ) = 7;


%=== Interpolation in angle for 18 and 36

ifre = find( floor(freqs) == 18 | floor(freqs) == 36 );


emiss_interp_v=nan(na,2);
emiss_interp_h=nan(na,2);

theta0  = 0;
theta53 = 53;

for f = 1:2

  indf = ifre(f); 
  ev53 = e_clim(ifre(f),1);
  eh53 = e_clim(ifre(f),2);
            
  for a = 1:na
    
    theta  = angles(ifre(f),a);

    %=  Calculation of the emis at theta=0� with a multilinear regression
    e0 = a0_k0(class1,f)+a0_k1(class1,f)*ev53+a0_k2(class1,f)*eh53;

    %=  Reading of the polynomial coefficients for ev et eh
    a0 = a0_eveh(class1,f);
    a1 = a1_eveh(class1,f);
    a2 = a2_eveh(class1,f);
    a3 = a3_eveh(class1,f);
    b0 = b0_eveh(class1,f);
    b1 = b1_eveh(class1,f);
    b2 = b2_eveh(class1,f);
    b3 = b3_eveh(class1,f);

    %= Vertical polarization
    S1_v      = ((theta-theta53)/(theta0-theta53)) * ((e0-a0)/a0);
    em53_v    = a3*(theta53^3) + a2*(theta53^2) + a1*theta53 + a0;
    S2_v      =((theta-theta0)/(theta53-theta0))*((ev53-em53_v)/em53_v);
    S_v       = 1 + S1_v + S2_v;
    emtheta_v = a3*(theta^3) + a2*(theta^2) + a1*theta + a0;
    emiss_interp_v(a,f) = S_v * emtheta_v;     

    %= Horizontal polarization 
    S1_h      = ((theta-theta53)/(theta0-theta53)) * ((e0-b0)/b0);
    em53_h    = b3*(theta53^3) + b2*(theta53^2) + b1*theta53 + b0;
    S2_h      =((theta-theta0)/(theta53-theta0))*((eh53-em53_h)/em53_h);
    S_h       = 1 + S1_h + S2_h;
    emtheta_h = b3*(theta^3) + b2*(theta^2) + b1*theta + b0;
    emiss_interp_h(a,f) = S_h * emtheta_h; 

  end

end



%= Emissivity V has to be larger or equal to emissivity H


for f = 1:2

  for a = 1:na
    
    if (emiss_interp_v(a,f) < emiss_interp_h(a,f))
       emiss_interp_v(a,f) = (emiss_interp_v(a,f) + emiss_interp_h(a,f))/2.;
       emiss_interp_h(a,f) =  emiss_interp_v(a,f);
    end

  end

end



%= for freqs <= 18GHz using angular dependence of 18GHz
%      freqs >= 18GHz using angular dependence of 36GHz

np   = 2; % v and h
emis = nan( np,nf, na);


for a = 1:na

  for f = 1:nf

     if floor(freqs(f)) <= 18
       s = 1;
     else
       s = 2;
     end

     emis(1,f,a) = e_clim(f,1) * emiss_interp_v(a,s) /e_clim(ifre(s),1);
     emis(2,f,a) = e_clim(f,2) * emiss_interp_h(a,s) /e_clim(ifre(s),2);

  end

end

else

np   = 2; % v and h
emis = nan( np,nf, na);


for a = 1:na

  for f = 1:nf

     if floor(freqs(f)) <= 18
       s = 1;
     else
       s = 2;
     end

     emis(1,f,a) = e_clim(f,1);
     emis(2,f,a) = e_clim(f,2);

  end

end

end



%=  Emissivity cannot be larger than 1 (this is not in TELSEM2, but for
%   large angles this can happen here with antenna integration, 
%    limits have to be set

%ind = find(emiss_interp_v > 1);
%emiss_interp_v(ind) = 1;

%ind = find(emiss_interp_h > 1);
%emiss_interp_h(ind) = 1;






return




%==========================================================================


function [ emis, clas ] = emis_land_surface_table( lon, lat, month, dir_input, freq )
% at 004 deg resolution

edir	= [ dir_input, '/EmisLand/Res_004/' ];
efile	= [ edir, '/land_emis_hdfillint_', sprintf('%02.0f',freq), 'GHz_53deg_', sprintf('%02.0f',month), '.mat'];


flon    = 0.018:0.036:359.982;
flat    = -89.982:0.036:89.982;

m       = matfile(efile); 

%= lon to 0-360

if lon < 0
  lon = lon + 360;
end

[ ~, ilo ] = min( abs( flon - lon ) );
[ ~, ila ] = min( abs( flat - lat ) );


emisv  = double(m.data( ilo, ila, 1 ))/1e3;
emish  = double(m.data( ilo, ila, 2 ))/1e3;
emis   = [emisv emish];


if emis(1) == 0 | emis(2) == 0

  %= looking for a closer one within +-100

  gap = 200;
  io = [ilo-gap:ilo+gap];
  io = io( io>=1 & io <= 10000);
  ia = [ila-gap:ila+gap];
  ia = ia( ia>=1 & ia <= 5000);
  aux  = double(m.data( io, ia, 1 ))/1e3;   
  aux  = mean(aux(aux>0));

  if ~isnan(aux)

    emis(1) = aux;
    aux  = double(m.data( io, ia, 2 ))/1e3;   
    %aux  = double(data( io, ia, 2 ))/1e3;   
    aux  = mean(aux(aux>0));
    
    if isnan(aux)

      emis(1) = -999;
      emis(2) = -999;

    else
    
      emis(2) = aux;

    end

  else

    emis(1) = -999;
    emis(2) = -999;

  end

end


if nargout == 2

  efile	= [ dir_input, '/EmisLand/Res_004/land_emis_hd_class_', sprintf('%02.0f',month), '.mat'];
  m       = matfile(efile); 

  clas  = double(m.data( ilo, ila));

end





return




