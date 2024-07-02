%-------------------------------------------------------------------------------
%
% SUBMODULE   forward_model_core
%
%    A function to simulate vertically and horizontally polarized 
%    top-of-atmosphere brightness temperatures (BTs) at a specified 
%    frequency and observing viewing zenith angle for a given 
%    atmosphere and underlaying surface.
%    This is the "core" function to perform the radiative
%    transfer for one frequency and zenith angle, It can
%    be used alone, or it can be called by forward_model_ice_free, 
%    a wrapper around this function to provide the weigting functions
%    and a format that can be called from the inversion function oem. 
%    
%    The forward modelling is intended for sea and sea ice surfaces,
%    but it can also deal with land surfaces, although with a simplified
%    estimation of the land surface emissivity based on a climatology.
%    This is to allow the treatment of coastal areas in the framework of a 
%    sea forward modelling exercise. 
%
%    For ice-free surface a combination of land and sea BTs is possible,
%    with surf_input.LSP giving the fraction of land in the scene. For
%    ice surfaces the surf_input.SIC gives the sea ice concentration.  
%    
%    NOTICE:
%
%	(1) atmospheric profiles need to be given as the value
%	at each specfic level, i.e., not an integrated value for
%	a layer. The code internally derived an integrated value 
%	for layer i to i+1 as the average value of i and i+1.
%
%	(2) the scalar CWVC is an optional input and can be used
%	to scale the RHO profile. This is to facilitate
%	the retrieval of a CWVC scalar, assuming that the RHO
%	vertical distribution is already correct. If CRHO is the 
%	column integration of RHO, RHO is scaled by a factor CWVC/CRHO 
%	conserving the original vertical distribution of RHO. In other
%	words, RHO at layer i becomes RHO * CWVC/CRHO. If no adjustment
%	is needed, do not included a CWVC field in atmos_input.
%
%	(3) the same applies to LW and CLWC.
%
%    REFERENCES:
%
%	This RT code is based on works published by many researchers, the 
%	main references are given below. Please contact carlos.jimenez@estellus.fr
%	if you noticed that a mssing reference should be added here.
% 
%        For computing ocean surface emissivity at microwave frequencies
%        and parametrizing the corresponding look-up table
%
%        1.    Meissner, T.; F. Wentz and D. Le Vine, Aquarius Salinity Retrieval Algorithm Theoretical Basis Document (ATBD),
%              End of Mission Version; RSS Technical Report 120117; December 1, 2017;
%              Available online at ftp://podaac-ftp.jpl.nasa.gov/allData/aquarius/docs/v5/AQ-014-PS-0017_Aquarius_ATBD-EndOfMission.pdf. 
%
%        2.    Meissner, T, F. Wentz, and D, Le Vine, 2018,
%              The Salinity Retrieval Algorithms for the NASA Aquarius Version 5 and SMAP Version 3 Releases,
%              Remote Sensing 10, 1121, doi:10.3390/rs10071121.
%
%        3.    Meissner, T. and F. Wentz, The complex dielectric constant of pure and sea water from microwave satellite observations,
%              IEEE TGRS, 2004, 42(9), 1836  1849, doi:10.1109/TGRS.2004.831888.
%
%        4.     Meissner, T. and F. Wentz, The emissivity of the ocean surface between 6 and 90 GHz
%              over a large range of wind speeds and Earth incidence angles,
%              IEEE TGRS, 2012, 50(8), 3004  3026, doi: 10.1109/TGRS.2011.2179662.
%
%        5.     Meissner, T., F. Wentz, F. and L. Ricciardulli, The emission and scattering of L-band microwave radiation
%              from rough ocean surfaces and wind speed measurements from Aquarius,
%	       J. Geophys. Res. Oceans, 2014, 119, doi:10.1002/2014JC009837.
%
%	For a climatology of land emissivity
%	
%	1.	Aires, F., C. Prigent, F. Bernardo, C. Jimenez, R. Saunders, and P. Brunel, A Tool to Estimate
%		 Land-Surface Emissivities at Microwave frequencies (TELSEM) for use in numerical weather 
%		 prediction, Q. J. R. Meteorol. Soc., 137:690-699, 2011. 
%
%       For computing ice surface emissivity at microwave frequencies
%       and parametrizing the corresponding look-up table	
%    
%	1.	Pedersen, L. F., & Saldo, R. (2016). Sea ice concentration (SIC) round robin data package, 
%		sea ice climate initiative: Phase 2 (Tech. Rep. No. SICCI-RRDP-07-16 Version: 1.4). DTU: ESA.
%
%
%		 
% FORMAT   [ bts, emis ]  = forward_model_core
%				( atmos_input, surf_input, sensor_input, aux_input )
%        
% OUT   bts		row vector 		Rayleigh-Jeans brightness temperatures
%			   V-pol [K]		vertical polarization   (1st column)
%			   H-pol [K]		horizontal polarization (2nd column)
%
%	emis 	        structure array         Structure containing the emissivities
%						for the different surface types:
%
%			   lan[-]		land emissivity, row vector of similar
%						size to bts
%			   sea[-]		sea emissivity, row vector of similar
%						size to bts
%			   ice[-]		ice emissivity, row vector of similar
%						size to bts
	       		
%
% IN    atmos_input	structure array		Structure containing the
%						atmospheric profiles:
%
%			   T	[K]		Temperature
%			   P	[mbar]		Pressure
%			   RHO  [g/m3]		Humidity
%			   LW   [g/m3]		Liquid Water 
%						
%						and the OPTIONAL altitude-
%						integrated contents:
%
%			   CWVC	[kg/m2]		Column Water Vapour Content 
%			   CLWC	[kg/m2]		Cloud Liquid Water Content	
%						
%
%			   
%	surf_input	structure array		Structure containing 
%						the surface parameters for
%						a sea simulation:			
%
%			    SST	[K]		Surface Temperature
%			    SSS [g/kg]		Sea Surface Salinity
%			    OWS	[m/s]		Wind velocity
%			    UWS	[m/s]		U wind component
%			    VWS	[m/s]		V wind component
%			    SLP [mbar]		Mean Sea Level Pressure
%	
%						the surface parameters
%						for a land simulation:
%
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
%
%						the surface parameters
%						for a sea-ice simulation:
%
%			    IST	[K]		Ice Surface Temperature 
%			    T2M [K]		2-meter Air Temperature
%			    SIC [-]		Sea Ice Concentration 0-1
%			    SND [m]		Sea ice snow depth
%			    XYI [-]		Ice age
%						  <= 1 considered as first-year ice
%						  >  1 considered as multi-year ice 
%			    LAT [degrees]	Latitude, 90S to 90N
%			    ILP [mbar]          Land Level Pressure
%			    EMI_ICE_COV [0-1]   If 0 ice emis is generated without
%						adding variability, i.e., for a given set
%						of conditions only the mean ice emis
%						value is generated. If 1, the mean
%						value is randomly changed to fullfill
%						the ice emis statistics derived from
%						the ice emis database.
%
%						and optional fields
%						to pass emissivity:
%
%			    EMI_SEA [-]		Vector with sea emissivity
%						for the V-pol and H-pol
%						for the given frequency. 
%						If this field exist, the
%						internal emissivity values
%						are ignored.
%		            EMI_LAN [-]		As EMI_SEA but for land. 
%		            EMI_ICE [-]		As EMI_SEA but for ice. 
%
%					
%
%	sensor_input	structure array		Structure containing 
%						the sensor parameters:
%
%			    F	[GHz]		Frequency		
%			    ZA	[degrees]	Viewing zenith angle
%			    AA	[degrees]	Azimuth zenith angle
%			    H	[km]		Sensor height from ground
%			    TOA []		Boolean indicating if the sensor radiative
%						transfer should be done only at the surface
%						(0) or contain also the atmospheric 
%						contributions (1).	
%
%	
%	dir_input	string			Folder where the following
%						auxiliary files need to be
%						placed:		         
%			    
%			    land_emis_fXX.mat		Land surface emissivity
%							for frequency XX, placed
%							at dir_input/EmisLand/    
%
%-------------------------------------------------------------------------------
% Project:	  CIMR Algorithm Performance Evaluation
% Package:	  CIMR Scientific Work Bench
% Developer:	  Estellus 
% Contact:	  carlos.jimenez@estellus.fr 
% Initiated:	  2019-01-08
%-------------------------------------------------------------------------------


function [ bts, emis ] = forward_model_core( atmos_input, surf_input, sensor_input, dir_input )


%==== MAIN-CODE ============================================================



%==== Inputs checking


do_sea  = 0;
do_ice  = 0;
do_lan  = 0;


if isfield( surf_input, 'SST' ) 
  do_sea = 1;
end
if isfield( surf_input, 'LST' )
  do_lan = 1;
end
if isfield( surf_input, 'IST' )
  do_ice = 1;
end




if do_sea & do_lan
  if ~isfield( surf_input, 'LSP' ) 
    osfi_error('Land and sea cell, but no surf_input.LSP land percentange in cell');
  end
end

if do_sea + do_lan + do_ice > 1

  aux_var     = {'EMIV_FQ1','EMIH_FQ1','EMIV_FQ2','EMIH_FQ2','EMIV_FQ3','EMIH_FQ3','EMIV_FQ4','EMIH_FQ4','EMIV_FQ5','EMIH_FQ5'};

  for a = 1:length(aux_var)  

    if isfield( surf_input, aux_var{a} ) 
      osfi_error('The RT codes cannot be used to derive emissivity jacobians from mixed scenes');
    end
  
  end

end

%= ingesting emis or internally
 

if isfield( surf_input, 'EMI_SEA' )
  do_emi_sea = 0;
else 
  do_emi_sea = 1;
end

if isfield( surf_input, 'EMI_ICE' )
  do_emi_ice = 0;
else 
  do_emi_ice = 1;
end

if isfield( surf_input, 'EMI_LAN' )
  do_emi_lan = 0;
else 
  do_emi_lan = 1;
end




%= checking angle

if sensor_input.ZA > 89 | sensor_input.ZA < 0
  osfi_error('Zenith angle outside 0-89 range');
end



%= checking frequency separation to channel center
%    ID	    Channel frequency (GHz)
%    1	    1.41
%    2	    6.925
%    3	    10.65
%    4	    18.7
%    5	    36.5

afreq				= [1.410 6.925 10.650 18.700 36.500];
[ id_channel, rfreq, deltafreq ]= nearest_in_vector( afreq, sensor_input.F);
rfreq				= floor(rfreq);


if abs(deltafreq) > 10 & sum(do_emi_sea+do_emi_lan+do_emi_ice)>0
  osfi_error('Frequency too far from CIMR centre channels, emissivity parameterizations not valid');
end


%= checking atmos inputs
%  TBD: some basic checking for ranges


%= checking surface inputs
%  TBD: some basic checking for ranges





%===== Scenes mixing
%      assuming linear mixing in emissivity * ts

if do_sea == 1 & do_lan == 0 & do_ice == 0

  % sea
  pc_sea = 1;
  pc_lan = 0;
  pc_ice = 0;

elseif do_sea == 0 & do_lan == 1 & do_ice == 0

  % sea
  pc_sea = 0;
  pc_lan = 1;
  pc_ice = 0;

elseif do_sea == 0 & do_lan == 0 & do_ice == 1

  % sea
  pc_sea = 0;
  pc_lan = 0;
  pc_ice = 1;


elseif do_sea == 1 & do_lan == 1 & do_ice == 0

  % sea and land
  pc_lan = surf_input.LSP;
  pc_sea = 1 - pc_lan;
  pc_ice = 0;

elseif do_sea == 1 & do_lan == 0 & do_ice == 1

  % sea and ice
  pc_lan = 0;
  pc_ice = surf_input.SIC;;
  pc_sea = 1 - pc_ice;

elseif do_sea == 0 & do_lan == 1 & do_ice == 1

  % land and ice
  pc_sea = 0;
  pc_ice = surf_input.SIC;;
  pc_lan = 1 - pc_ice;


elseif do_sea == 1 & do_lan == 1 & do_ice == 1

  % sea and land
  pc_lan = surf_input.LIP;
  pc_sea = (1 - pc_lan) * ( 1 - surf_inpu.SIC);
  pc_ice = (1 - pc_lan) * surf_inpu.SIC;



else

  osfi_error('The scene mixing based on surf_input cannot be generated, check fields LST, SST and IST'); 

end



%==== Calculating atmospheric contribution
%     if TOA BTs need to be calculated

if sensor_input.TOA


  %===== Unit conversions (if needed) and from structure
  %      to original code variable names

  % Incidence Angle in degrees
  ANGLE	    = sensor_input.ZA;

  % Maximum Altitude in kms     
  if do_sea == 1
    surpre	     = surf_input.SLP; 
  elseif do_lan == 1
    surpre	     = surf_input.LLP; % 
  elseif do_ice == 1
    surpre	     = surf_input.ILP; % 
  end


  % Frquency in GHz
  FREQ	   = sensor_input.F;

  % Pressure profile in mbar
  PRES	   = atmos_input.P;

  % Number of atmospheric levels
  NLEV	   = length(PRES);

  % Sensor height
  SH	   = sensor_input.H;  % km


  % temperature profile in K
  TEMP	   = atmos_input.T;

  % humidity in g/m3
  RHO	   =  atmos_input.RHO; 

  % liquid water in g/m3
  CLD	   = atmos_input.LW; 

  % ozone number density in mol/m3
  % hard-coded to zero




  %=== altitude profile in Km

  H = p2z_barometric( PRES, TEMP, 0 );

  if H(1) > H(end)
    osfi_error('Atmospheric profiles should be given from ground up, this one looks like from top of the atmosphere to ground');
  end	



  %===== Adjustment of RHO profile to CWVC


  if isfield( atmos_input, 'CWVC' )

    %= column integrating RHO

    CRHO = RHO(1) * 0.5 * (H(2) - H(1)) ;

    for p = 2:NLEV-1

      hu = H(p+1) - (H(p+1)-H(p))/2;
      hl = H(p-1) + (H(p)-H(p-1))/2;

      CRHO = CRHO + RHO(p) * ( hu - hl );  

    end

    CRHO = CRHO + RHO(NLEV) * 0.5 * (H(NLEV) - H(NLEV-1));   %g/m2

    %= scaling RHO
    RHO = (atmos_input.CWVC/CRHO) * RHO;

  end

  

  %=== adjustment of LW profile to CLWC

  if isfield( atmos_input, 'CLWC' )

    %= column integrating LW

    CLW = CLD(1) * 0.5 * (H(2) - H(1)) ;

    for p = 2:NLEV-1

      hu = H(p+1) - (H(p+1)-H(p))/2;
      hl = H(p-1) + (H(p)-H(p-1))/2;

      CLW = CLW + CLD(p) * ( hu - hl );  

    end

    CLW = CLW + CLD(NLEV) * 0.5 * (H(NLEV) - H(NLEV-1));

    %= scaling LW

    CLD = (atmos_input.CLWC/CLW) * CLD;


  end



  %===== Cosmic background

  % including zero-point fluctuation in cosmic background compensates
  % for non-linearity of Planck function

  rcons		    = .0479923;
  efac		    = exp(rcons*FREQ/2.73);
  TBC		    = .5*rcons*FREQ*(efac+1.)/(efac-1.);




  %= checking if sensor inside the atmosphere.
  %  If yes, we reduce the atmosphere from
  %  sensor height to ground. If not, we 
  %  do RT for the whole atmosphere and 
  %  sensor height does not affect this
  %  1-D RT calculations


  if SH > H(NLEV)

    INLEV = NLEV;

  else

    INLEV   = find( H <= SH );
    [~,ind] = max( H(INLEV) );
    INLEV   = INLEV(ind);

  end



  %===== Opacity calculation
  %      moving ABSN2 and ABLIQ here instead of function to speed up
  %      processing



  TOTOP1 = 0.;
  for I=2:INLEV

    % use the 'absorption-of-averages' method to compute optical
    % depth of each slab; see M.J. Schwartz, Ph.D. thesis pp. 84-87.

    TAV = (TEMP(I) + TEMP(I-1))/2.;
    PAV = sqrt(PRES(I)*PRES(I-1));
    WVAV = (RHO(I) + RHO(I-1))/2.;
    WLAV = (CLD(I) + CLD(I-1))/2.;
    ABSCOEF =  o2abs(TAV,PAV,WVAV,FREQ) + abh2o(TAV,PAV,WVAV,FREQ) + ...
    absn2(TAV,PAV,FREQ) + abliq(WLAV,FREQ,TAV);
    OPACITY(I) = ABSCOEF*abs(H(I)-H(I-1));
    TOTOP1 = TOTOP1 + OPACITY(I);

  end



  %=====  Transmissivity, TBdown, and TBup

  tb1	= 0.;
  tb2	= 0.;
  tra	= 1.;


  CTHN	 = cos(ANGLE/57.296);
  SECANT = 1./CTHN;

  % Loops on the atmosphericlevels
  for I=2:INLEV
    TAV = (TEMP(I) + TEMP(I-1))/2.;
    % trace path 1 using integral form of RTE
    EM = tra;
    tra = tra*exp(-SECANT*OPACITY(I));
    tb1 = tb1 + TAV*(EM-tra);
    % trace path 2 using differential form of RTE
    TRAN_SLAB = exp(-SECANT*OPACITY(I));
    tb2 = TAV + TRAN_SLAB*(tb2-TAV);
  end 



  %= atmospheric TBdown
  tbdown	= tb1 + tra * TBC;

  % unpolarized atmosphere so in V H 3rd 4th base
  tbdown        = [ tbdown tbdown 0 0];


  %= atmospheric TBup
  tbup	        = tb2;

  % unpolarized atmosphere so in V H 3rd 4th base
  tbup        = [ tbup tbup 0 0];

else

  tbdown  = [0 0 0 0];
  tbup    = [0 0 0 0];
  tra     = 1;

end




%=== Surface emissivity

% some rt check while 
% fixing testing 
% emissivity 

do_rt = 1;


% surface temperature
if do_sea == 1
  SST = surf_input.SST;
  SSS = surf_input.SSS;
  OWS = surf_input.OWS;
  UWS = surf_input.UWS;
  VWS = surf_input.VWS;
else
  SST = 0;
end
if do_lan == 1
  LST = surf_input.LST;
else
  LST = 0;
end
if do_ice == 1
  IST = surf_input.IST;
else
  IST = 0;
end



if do_sea

  if do_emi_sea
    
    % owdrs  =  wind direction relative to 8 azimuth angles

    % calculating wind direction relative
    % to satellite azimuth direction:

    % Instrument simulation (SCEPS) convention on all azimuth angles is positive 
    % counterclockwise from due east, always.
    % For wind, wind direction in this convention in degrees (ie the direction 
    % TOWARD which wind is blowing), is atan2d(v,u)
    %
    % So the wind direction (downwind direction) relative to emission direction
    % (toward satellite) is npi2pi(atan2d(v,u) - hbs_eaa)
    % 
    % Now satellite LOOK direction is opposite emission direction so that
    % the wind direction (downwind direction) relative to satellite look direction 
    % (toward satellite) is npi2pi(atan2d(v,u) - (hbs_eaa+180 deg))

    OWD = mod( npi2pi(atan2d(VWS,UWS) - sensor_input.AA), 360 );    
    [emis_sea, refl_sea ] = forward_model_sea_surface( SST, SSS, OWS, OWD, rfreq, sensor_input.ZA, tra);

  else
    emis_sea = surf_input.EMI_SEA;
  end

end



if do_lan

  % no incidence or azimuthan angle variations yet
  % we use old parameterization for the moment

  if do_emi_lan

    if rfreq == 1
      aangle				= [40 53 55 57];
      [ id_angle, rangle ]= nearest_in_vector( aangle, sensor_input.ZA);
    else
      aangle				= [53 55 57];
      [ id_angle, rangle ]= nearest_in_vector( aangle, sensor_input.ZA);
    end

    [ emis_lan ] = emis_land_surface_table_angle( surf_input.LON, surf_input.LAT, surf_input.MONTH, dir_input, rfreq, rangle);   
    if emis_lan(1) == -999 | emis_lan(2) == -999
      osfi_error('A land emissivity could not be found, it looks like the location was considered as 100% sea in the database');
    end

  else

    emis_lan = surf_input.EMI_LAN;  

  end


end

if do_ice


  % no incidence or azimuthan angle variations yet
  % we use old parameterization for the moment

  % no incidence or azimuthan angle variations yet
  % we use old parameterization for the moment

  if do_emi_ice
 
    if rfreq == 1
      aangle	= [40 55];
      [ id_angle, rangle ]= nearest_in_vector( aangle, sensor_input.ZA);
    else
      rangle = 55;
    end

    if surf_input.EMI_ICE_COV
      osfi_error('surf_input.EMI_ICE_COV==1 calculations require further implememntaions');
      [ emis_ice, diagcov_ice ] = emis_ice_surface_interp( surf_input.LAT, surf_input.XYI, surf_input.T2M, surf_input.SIT, dir_input, rfreq, rangle, surf_input.EMI_ICE_COV );
    else
      [ emis_ice, diagcov_ice ] = emis_ice_surface_table( surf_input.LAT, surf_input.XYI, surf_input.T2M, surf_input.SND, dir_input, rfreq, rangle );
    end

  else

    emis_ice    =  surf_input.EMI_ICE;
    diagcov_ice = nan(size(emis_ice)); 

  end

  if emis_ice(1) == -999 | emis_ice(2) == -999 
    do_rt = 0;
  end

end




%===== TB contribtions from sea, ice, land


if do_rt

%= sea surface contribution

if do_sea

  % if sensor_input.TOA == 0 then tbdown is 0 so only
  % surface component is calcualted 

  bts_sea  = tbdown .* tra .* refl_sea' + SST * tra * emis_sea';

else

  bts_sea = [0 0 0 0];

end



%= land surface contribution


if do_lan

  bts_lan = tbdown .* tra .* ( 1 - emis_lan)  + LST * tra * emis_lan;

  % no 3rd and 4th components over land so bts to zero
  % and emis_lan at -999

  bts_lan(3:4) = 0; 

  % no azimuthal dependence for the moment so same value at all aangles


else

  bts_lan = [0 0 0 0];

end

%= ice surface contribution


if do_ice

  bts_ice = tbdown .* tra .* ( 1 - emis_ice)  + IST * tra * emis_ice;

  % no 3rd and 4th components over ice so bts to zero
  % and emis_lan at -999

  bts_ice(3:4) = 0; 

  % no azimuthal dependence for the moment so same value at all aangles

else

  bts_ice = [0 0 0 0];

end



%===== TB top of atmosphere V-pol and H-pol
% for single and mixed scenes


bts = tbup +  pc_sea * bts_sea + pc_lan * bts_lan + pc_ice * bts_ice;



% no 3rd and 4th components over land and ice so
% if no sea component we place the tbup to zero
% so we do not get in the 3rd and 4th tbup 

if pc_sea == 0
  bts(3:4) = -999;
end  

if 0
  % for debugging
  FREQ=rfreq;
  if do_ice
    fprintf('%.2f GHz Trans=%.4f Tbdown=%.2f Tbup=%.2f Ts=%.2f EmisV=%.2f EmisH=%.2f TbV=%.2f TbH=%.2f\n ',FREQ,tra,tbdown,tbup,SST,emis_ice(1),emis_ice(2),bts_ice(1),bts_ice(2));
  end
end


else

  bts = -999 * ones(1,4);;

end


emis.lan  = [];
emis.sea  = [];
emis.ice  = [];
if do_ice
  emis.ice         = emis_ice;
  emis.ice_diagcov = diagcov_ice;
end
if do_lan
  emis.lan = emis_lan;
end
if do_sea
  emis.sea = emis_sea;
end

return




%==== SUB-FUNCTION =======================================================


%=== dummy function before integration to main package
%    TO BE REMOVED later	   

function osfi_error(msgd)

error(msgd);

return




%==========================================================================

function absor=abh2o(T,P,RHO,F)

% PURPOSE- COMPUTE ABSORPTION COEF IN ATMOSPHERE DUE TO WATER VAPOR
% This version should not be used with a line list older than June 2018,
% nor the new list with an older version of this subroutine.
% 
%  CALLING SEQUENCE PARAMETERS-
%    SPECIFICATIONS
%      REAL T,P,RHO,F,ABH2O
%      NAME    UNITS    I/O  DESCRIPTON            VALID RANGE
%      T       KELVIN    I   TEMPERATURE
%      P       MILLIBAR  I   PRESSURE              .1 TO 1000
%      RHO     G/M**3    I   WATER VAPOR DENSITY
%      F       GHZ       I   FREQUENCY             
%      ABH2O   NEPERS/KM O   POWER ABSORPTION COEFFICIENT

%   Multiply ABH2O by 4.343 to obtain dB/km.
%   Line parameters will be read from file h2o_list.asc; intensities should
%   include the isotope abundance factors.
%   This version uses a line-shape cutoff.

%   REVISION HISTORY-
%     DATE- OCT.6, 1988 EQS AS PUBL.: P.W. Rosenkranz, CHAP. 2 in 
%     ATMOSPHERIC REMOTE SENSING BY MICROWAVE RADIOMETRY 
%     (M.A. Janssen, ed., 1993) (http://hdl.handle.net/1721.1/68611).
%     OCT.4, 1995  PWR- USE CLOUGH'S DEFINITION OF LOCAL LINE
%                   CONTRIBUTION,  HITRAN INTENSITIES, ADD 7 LINES.
%     OCT. 24, 95  PWR -ADD 1 LINE.
%     JULY 7, 97   PWR -SEPARATE COEFF. FOR SELF-BROADENING, 
%                       REVISED CONTINUUM.
%     Mar. 2, 2003   PWR - LINE SHIFT
%     Nov. 3, 2012 intensities at base T=296K, get line param. from file.
%     Aug. 6, 2015 read continuum param from the file also. 
%     June 19, 2018 changed file format, separate shift for self & foreign gas
%     Jan. 06, 2019 C. Prigent: change to matlab code and change in the reading
%     of the line parameters. Simpler for matlab. 
%   Initialization section
%   read line parameters as of 12/2018
NLINE = 1;
FL=[22.2351 183.3101 321.2256 325.1529 380.1974 439.1508 443.0183 448.0011 470.8890 474.6891 488.4901 556.9360 620.7008 658.0060 752.0331 916.1716];
S1=[0.1335E-13 0.2319E-11 0.7657E-13 0.2721E-11 0.2477E-10 0.2137E-11 0.4440E-12 0.2588E-10 0.8196E-12 0.3268E-11 0.6628E-12 0.1570E-08 0.1700E-10 0.9033E-12 0.1035E-08 0.4275E-10];
B2=[2.172 0.677 6.262 1.561 1.062 3.643 5.116 1.424 3.645 2.411 2.890 0.161 2.423 7.921 0.402 1.461];
W3=[0.2699E-02 0.2945E-02 0.2426E-02 0.2847E-02 0.2868E-02 0.2055E-02 0.1819E-02 0.2612E-02 0.2169E-02 0.2366E-02 0.2616E-02 0.3115E-02 0.2468E-02 0.3154E-02 0.3114E-02 0.2695E-02];
X=[0.76 0.77 0.73 0.64 0.54 0.69 0.70 0.70 0.73 0.71 0.75 0.75 0.79 0.73 0.77 0.79];
WS=[0.1329E-01 0.1478E-01 0.1065E-01 0.1395E-01 0.1440E-01 0.9060E-02 0.7960E-02 0.1301E-01 0.9700E-02 0.1124E-01 0.1358E-01 0.1424E-01 0.1194E-01 0.1384E-01 0.1358E-01 0.1355E-01];
XS=[1.20 0.78 0.54 0.74 0.89 0.52 0.50 0.67 0.65 0.64 0.72 1.00 0.75 1.00 0.84 0.48];
XH=[2.60 0.77 0.73 0.64 0.54 0.69 0.70 0.70 0.73 0.71 0.75 0.75 0.79 0.73 0.77 0.79];
XHS=[1.20 0.78 0.54 0.74 0.89 0.52 0.50 0.67 0.65 0.64 0.72 1.00 0.92 1.00 0.84 0.47];
SH=[-.3300E-04 -.7200E-04 -.1430E-03 -.1300E-04 -.7400E-04 0.5100E-04 0.1400E-03 -.1160E-03 0.6100E-04 -.2700E-04 -.6500E-04 0.1870E-03 0.0000E+00 0.1760E-03 0.1620E-03 0.0000E+00];
SHS=[0.8140E-03 0.1730E-03 0.2780E-03 0.1325E-02 0.2400E-03 0.1650E-03 -.2290E-03 -.6150E-03 -.4650E-03 -.7200E-03 -.3600E-03 -.1693E-02 0.6870E-03 -.1496E-02 -.8780E-03 0.5210E-03];
REFTLINE = 296.; %reference T for lines
%    read continuum parameters; units: Kelvin, 1/(km*mb^2*GHz^2)
REFTCON=300.;
CF=0.59E-09;
XCF=3.;
CS=0.142E-07;
XCS=7.5;
 

if (RHO>0.)        
   PVAP = RHO*T/216.68;
   PDA = P -PVAP;
   DEN = 3.344E16*RHO;  
   % CONTINUUM TERMS
   TI = REFTCON/T;
   % Xcf and Xcs include 3 for density & stimulated emission
   CON = (CF*PDA*TI^XCF + CS*PVAP*(TI^XCS))*PVAP*F*F;
   % ADD RESONANCES
   TI = REFTLINE/T;
   TI2 = TI^2.5;
   SUM = 0.;

   %= modification for CIMR work
   % leaving only the 22 GHz line

   %NLINE = 1; 
   %for I=1:NLINE

   I = 1;
       WIDTHF = W3(I)*PDA*(TI^X(I));
       WIDTHS = WS(I)*PVAP*(TI^XS(I));
       WIDTH = WIDTHF + WIDTHS;
       SHIFTF = SH(I)*PDA*(TI^XH(I));
       SHIFTS = SHS(I)*PVAP*(TI^XHS(I));
       SHIFT = SHIFTF + SHIFTS;
       WSQ = WIDTH^2;
       % line intensities include isotopic abundance
       S = S1(I)*TI2*exp(B2(I)*(1.-TI)); 
       DF(1) = F - FL(I) - SHIFT;
       DF(2) = F + FL(I) + SHIFT;
       % USE CLOUGH'S DEFINITION OF LOCAL LINE CONTRIBUTION
       BASE = WIDTH/(562500. + WSQ);
       % DO FOR POSITIVE AND NEGATIVE RESONANCES
       RES = 0.;
       for J=1:2
          if (abs(DF(J))<750.) 
             RES = RES + WIDTH/((DF(J)^2)+WSQ) - BASE;
          end  
       end
       SUM = SUM + S*RES*((F/FL(I))^2);
       absor = .3183E-4*DEN*SUM + CON;   

   %end

else

   absor = 0.;

end

return



%==========================================================================

function absor=abliq(clo,f,t)

%     ORIGINALY ABLIQ12
%     COMPUTES POWER ABSORPTION COEFFICIENT IN NEPERS/KM 
%     BY SUSPENDED CLOUD LIQUID WATER DROPLETS. MULTIPLY ABLIQ BY
%     4.343 TO CONVERT TO DB/KM.

%     ARGUMENTS (INPUT):
%     clo IN G/M**3
%     f IN GHZ     (VALID FROM 0 TO 1000 GHZ)
%     t IN KELVIN

%     REVISION HISTORY:
%     PWR 6/5/15   using dilec12 for complex dielectric constant
%     Jan. 06, 2019 Kate: change to matlab code and change in the reading
%     of the line parameters. Simpler for matlab. 
%     01/2019 transformed in matlb (C. Prigent)

if (clo<=0. | t<233.) 
    absor = 0.;
else
    eps=dilec(f,t);
    re = (eps-1.)/(eps+2.);
    absor = -.06286*imag(re)*f*clo;
    %fprintf('%.2f %.2f %.2f %f%+fj %.2f\n',f,clo,t,real(eps),imag(eps),absor);
end

return




%==========================================================================

function absor=absn2(t,p,f)
%     Copyright (c) 2002 Massachusetts Institute of Technology
%     ABSN2 = COLLISION-INDUCED POWER ABSORPTION COEFFICIENT 
%     (NEPER/KM) IN AIR ("dry continuum", mostly due to N2-N2, 
%     but also contributions from O2-N2 and O2-O2)
%     T = TEMPERATURE (K)
%     P = DRY AIR PRESSURE (MB)
%     F = FREQUENCY (GHZ)(valid 0-2000 GHz)
%       Multiply ABSN2 by 4.343 to obtain dB/km.

%     5/22/02 4/14/05 6/23/18 P.Rosenkranz
%     01/2019 transformed in Matlab (C. Prigent) 

%     References:
%     Frequency dependence based on model by A. Borysow and L. Frommhold,
%      Astrophysical Journal, v.311, pp.1043-1057 (1986).
%     See Eq. 2.6 in Thermal Microwave Radiation - Applications 
%      for Remote Sensing (C. Maetzler, ed.) London, IET, 2006.
%     Amplitude increased by 14% based on analysis by
%      M. Tretyakov and A. Zibarova, JQSRT v.216, pp. 70-75 (2018)  

th = 300./t;
fdepen=.5 + .5/(1.+((f/450.)^2));
absor = 9.95e-14*fdepen*p*p*f*f*(th^3.22);
      


return




%==========================================================================

function diel=dilec(f,t)

%   Purpose: Computes the complex dielectric constant for liquid water,
%   with a negative imaginary part representing dissipation.

%   Complex logarithm is used here. It should be defined with
%   imaginary part in the range -pi to +pi.

%   Copyright © P.W. Rosenkranz  Apr. 15, 2014
%   Creative Commons license CC BY-SA
%   Modified for matlab 01/2019 (C. Prigent). Tested and it works!

%     inputs:
%      real f  ! frequency in GHz, 
%      real t ! Kelvin temperature
%     validated for 20<f<220 GHz at 248<t<273; 1<f<1000 GHz at 273<t<330.

tc = t - 273.15;
z = complex(0.,f);
theta = 300./t;

%  static dielectric constant model from Patek et al. (J.Phys.Chem.Ref.Data. v.38(1), 21 (2009).
kappa = -43.7527*theta^0.05 +299.504*theta^1.47-399.364*theta^2.11 +221.327*theta^2.31;
%  Debye term from 
%  W. Ellison, J. Phys. Chem. Ref. Data, 36, 1-18 (2007).
delta = 80.69715*exp(-tc/226.45);
sd = 1164.023*exp(-651.4728/(tc+133.07));
kappa = kappa -delta*z/(sd+z);
%  B band from P.W. Rosenkranz, IEEE Trans. Geosci. & Remote Sens. v.53(3) pp.1387-93 (2015).
delta = 4.008724*exp(-tc/103.05);
hdelta = delta/2.;
f1 = 10.46012+0.1454962*tc+6.3267156E-02*tc^2+9.3786645E-04*tc^3;
z1 = complex(-.75,1.)*f1;
z2 = complex(-4500.,2000.);
cnorm = log(z2/z1);
chip = hdelta*log((z-z2)/(z-z1))/cnorm;
chij = hdelta*log((z-conj(z2))/(z-conj(z1)))/conj(cnorm);
dchi = chip+chij-delta;
kappa = kappa + dchi;
diel=kappa;
%fprintf('%.2f %.2f %f%+fj \n',f,t,real(diel),imag(diel));

return



%==========================================================================

function absor=o2abs(t,p,rho,freq)

%  Copyright (c) 2009 Massachusetts Institute of Technology
%     RETURNS POWER ABSORPTION COEFFICIENT DUE TO OXYGEN IN AIR,
%     IN NEPERS/KM.  MULTIPLY O2ABS BY 4.343 TO CONVERT TO DB/KM.

%      5/1/95  P. Rosenkranz 
%      11/5/97  P. Rosenkranz - 1- line modification.
%      12/16/98 pwr - updated submm freq's and intensities from HITRAN96
%      8/21/02  pwr - revised width at 425
%      3/20/03  pwr - 1- line mixing and width revised
%      9/29/04  pwr - new widths and mixing, using HITRAN intensities
%                     for all lines
%      6/12/06  pwr - chg. T dependence of 1- line to 0.8
%      10/14/08 pwr - moved isotope abundance back into intensities, 
%                     added selected O16O18 lines.
%      5/30/09  pwr - remove common block, add weak lines.
%      12/18/14 pwr - adjust line broadening due to water vapor.


%     ARGUMENTS:
%     REAL TEMP,PRES,VAPDEN,FREQ
%     NAME    UNITS    DESCRIPTION        VALID RANGE
%     TEMP    KELVIN   TEMPERATURE        UNCERTAIN, but believed to be
%                                          valid for atmosphere
%     PRES   MILLIBARS PRESSURE           3 TO 1000
%     VAPDEN  G/M**3   WATER VAPOR DENSITY  (ENTERS LINEWIDTH CALCULATION
%                      DUE TO GREATER BROADENING EFFICIENCY OF H2O)
%     FREQ    GHZ      FREQUENCY          0 TO 900
%     REFERENCES FOR EQUATIONS AND COEFFICIENTS:
%     P.W. Rosenkranz, CHAP. 2 in ATMOSPHERIC REMOTE SENSING
%       BY MICROWAVE RADIOMETRY (M.A. Janssen, ed., 1993) 
%       (http://hdl.handle.net/1721.1/68611).
%     G.Yu. Golubiatnikov & A.F. Krupnov, J. Mol. Spect. v.217, 
%       pp.282-287 (2003).
%     M.Yu. Tretyakov et al, J. Mol. Spect. v.223, pp.31-38 (2004).
%     M.Yu. Tretyakov et al, J. Mol. Spect. v.231, pp.1-14 (2005).
%     B.J. Drouin, JQSRT v.105, pp.450-458 (2007).
%     D.S. Makarov et al, J. Mol. Spect. v.252, pp.242-243 (2008).
%     M.A. Koshelev et al, JQSRT, in press (2015).
%     line intensities from HITRAN2004.
%     non-resonant intensity from JPL catalog.
%     note:
%     1. The mm line-width and mixing coefficients are from Tretyakov et al;
%        submm line-widths from Golubiatnikov & Krupnov (except 
%        234 GHz from Drouin)
%     2. The same temperature dependence (X) is used for submillimeter 
%        line widths as in the 60 GHz band: (1/T)**X 
%     Local variables:
     NL=49;
%      LINES ARE ARRANGED 1-,1+,...37-,37+ IN SPIN-ROTATION SPECTRUM;
%      BY FREQUENCY IN SUBMM SPECTRUM.
     F=[118.7503, 56.2648, 62.4863, 58.4466, 60.3061, 59.5910, ...
         59.1642, 60.4348, 58.3239, 61.1506, 57.6125, 61.8002, ...
         56.9682, 62.4112, 56.3634, 62.9980, 55.7838, 63.5685, ...
         55.2214, 64.1278, 54.6712, 64.6789, 54.1300, 65.2241, ...
         53.5958, 65.7648, 53.0669, 66.3021, 52.5424, 66.8368, ...
         52.0214, 67.3696, 51.5034, 67.9009, 50.9877, 68.4310, ...
         50.4742, 68.9603, 233.9461, 368.4982, 401.7398, 424.7630, ...
         487.2493, 566.8956, 715.3929, 731.1866, ...
         773.8395, 834.1455, 895.0710];      
     S300=[0.2906E-14,0.7957E-15,0.2444E-14,0.2194E-14, ...
         0.3301E-14,0.3243E-14,0.3664E-14,0.3834E-14, ...
         0.3588E-14,0.3947E-14,0.3179E-14,0.3661E-14, ...
         0.2590E-14,0.3111E-14,0.1954E-14,0.2443E-14, ...
         0.1373E-14,0.1784E-14,0.9013E-15,0.1217E-14, ...
         0.5545E-15,0.7766E-15,0.3201E-15,0.4651E-15, ...
         0.1738E-15,0.2619E-15,0.8880E-16,0.1387E-15, ...
         0.4272E-16,0.6923E-16,0.1939E-16,0.3255E-16, ...
         0.8301E-17,0.1445E-16,0.3356E-17,0.6049E-17, ...
         0.1280E-17,0.2394E-17, ...
         0.3287E-16,0.6463E-15,0.1334E-16,0.7049E-14, ...
         0.3011E-14,0.1797E-16,0.1826E-14,0.2193E-16, ...
         0.1153E-13,0.3974E-14,0.2512E-16];
      BE=[.010, .014, .083, 0.083, .207, 0.207, .387, .387, .621,.621, ...
         .910, .910,1.255,1.255,1.654,1.654,2.109,2.109,2.618,2.618, ...
         3.182,3.182,3.800,3.800,4.474,4.474,5.201,5.201,5.983,5.983, ...
         6.819,6.819,7.709,7.709,8.653,8.653,9.651,9.651, ...
         .019, .048, .045, .044, .049, .084, .145, .136, .141, .145, .201];
%     WIDTHS IN MHZ/MB
      WB300=.56;
      X=.8;
      W300=[1.688, 1.703, 1.513, 1.491, 1.415, 1.408, ...
         1.353, 1.339, 1.295, 1.292, 1.262, 1.263, 1.223, 1.217, ...
         1.189, 1.174, 1.134, 1.134, 1.089, 1.088, 1.037,1.038, ...
         0.996,0.996,0.955,0.955,0.906,0.906,0.858,0.858,0.811,0.811, ...
         0.764,0.764,0.717, 0.717,0.669,0.669, ...
         1.65,1.64,1.64,1.64,1.60,1.60,1.60,1.60,1.62,1.47,1.47];
      Y300=[-0.0360, 0.2547, -0.3655,  0.5495, ...
         -0.5696,  0.6181, -0.4252,  0.3517, -0.1496,  0.0430, ...
         0.0640, -0.1605,  0.2906, -0.3730,  0.4169, -0.4819, ...
         0.4963, -0.5481,  0.5512, -0.5931,  0.6212, -0.6558, ...
         0.6920, -0.7208,  0.7312, -0.7550,  0.7555, -0.7751, ...
         0.7914, -0.8073,  0.8307, -0.8431,  0.8676, -0.8761, ...
         0.9046, -0.9092,  0.9416, -0.9423, ...
         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.];
      V=[0.0079, -0.0978,  0.0844, -0.1273, ...
         0.0699, -0.0776,  0.2309, -0.2825,  0.0436, -0.0584, ...
         0.6056, -0.6619,  0.6451, -0.6759,  0.6547, -0.6675, ...
         0.6135, -0.6139,  0.2952, -0.2895,  0.2654, -0.2590, ...
         0.3750, -0.3680,  0.5085, -0.5002,  0.6206, -0.6091, ...
         0.6526, -0.6393,  0.6640, -0.6475,  0.6729, -0.6545, ...
         0.680,  -0.660,   0.685,  -0.665,  ...
         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.];
     
     TEMP=t;
     PRES=p;
     VAPDEN=rho;
     FREQ=freq;
      TH = 300./TEMP;
      TH1 = TH-1.;
      B = TH^X;
      PRESWV = VAPDEN*TEMP/216.68;
      PRESDA = PRES -PRESWV;
      DEN = .001*(PRESDA*B + 1.2*PRESWV*TH);
      DFNR = WB300*DEN;
      SUM = 1.584E-17*FREQ*FREQ*DFNR/(TH*(FREQ*FREQ + DFNR*DFNR));

      %= modification for CIMR work
      % leaving only relevant frequencies

      %ind  = find( S300 > 2e-15 & F < 100 ); 
      ind  = [3     4     5     6     7     8     9    10    11    12    13    14    16]; 

      for K=1:NL%ind
          DF = W300(K)*DEN;
          FCEN = F(K);
          Y = DEN*(Y300(K)+V(K)*TH1);
          STR = S300(K)*exp(-BE(K)*TH1);
          SF1 = (DF + (FREQ-FCEN)*Y)/((FREQ-FCEN)^2 + DF*DF);
          SF2 = (DF - (FREQ+FCEN)*Y)/((FREQ+FCEN)^2 + DF*DF);
          SUM = SUM + STR*(SF1+SF2)*(FREQ/F(K))^2;
      end
      absor = 1.6097E11*SUM*PRESDA*TH^3;
      absor = max(absor,0.);


return



%==========================================================================

function absor=o2abs_fast(TEMP,PRES,VAPDEN,FREQ)

%  Copyright (c) 2009 Massachusetts Institute of Technology
%     RETURNS POWER ABSORPTION COEFFICIENT DUE TO OXYGEN IN AIR,
%     IN NEPERS/KM.  MULTIPLY O2ABS BY 4.343 TO CONVERT TO DB/KM.

%      5/1/95  P. Rosenkranz 
%      11/5/97  P. Rosenkranz - 1- line modification.
%      12/16/98 pwr - updated submm freq's and intensities from HITRAN96
%      8/21/02  pwr - revised width at 425
%      3/20/03  pwr - 1- line mixing and width revised
%      9/29/04  pwr - new widths and mixing, using HITRAN intensities
%                     for all lines
%      6/12/06  pwr - chg. T dependence of 1- line to 0.8
%      10/14/08 pwr - moved isotope abundance back into intensities, 
%                     added selected O16O18 lines.
%      5/30/09  pwr - remove common block, add weak lines.
%      12/18/14 pwr - adjust line broadening due to water vapor.


%     ARGUMENTS:
%     REAL TEMP,PRES,VAPDEN,FREQ
%     NAME    UNITS    DESCRIPTION        VALID RANGE
%     TEMP    KELVIN   TEMPERATURE        UNCERTAIN, but believed to be
%                                          valid for atmosphere
%     PRES   MILLIBARS PRESSURE           3 TO 1000
%     VAPDEN  G/M**3   WATER VAPOR DENSITY  (ENTERS LINEWIDTH CALCULATION
%                      DUE TO GREATER BROADENING EFFICIENCY OF H2O)
%     FREQ    GHZ      FREQUENCY          0 TO 900
%     REFERENCES FOR EQUATIONS AND COEFFICIENTS:
%     P.W. Rosenkranz, CHAP. 2 in ATMOSPHERIC REMOTE SENSING
%       BY MICROWAVE RADIOMETRY (M.A. Janssen, ed., 1993) 
%       (http://hdl.handle.net/1721.1/68611).
%     G.Yu. Golubiatnikov & A.F. Krupnov, J. Mol. Spect. v.217, 
%       pp.282-287 (2003).
%     M.Yu. Tretyakov et al, J. Mol. Spect. v.223, pp.31-38 (2004).
%     M.Yu. Tretyakov et al, J. Mol. Spect. v.231, pp.1-14 (2005).
%     B.J. Drouin, JQSRT v.105, pp.450-458 (2007).
%     D.S. Makarov et al, J. Mol. Spect. v.252, pp.242-243 (2008).
%     M.A. Koshelev et al, JQSRT, in press (2015).
%     line intensities from HITRAN2004.
%     non-resonant intensity from JPL catalog.
%     note:
%     1. The mm line-width and mixing coefficients are from Tretyakov et al;
%        submm line-widths from Golubiatnikov & Krupnov (except 
%        234 GHz from Drouin)
%     2. The same temperature dependence (X) is used for submillimeter 
%        line widths as in the 60 GHz band: (1/T)**X 
%     Local variables:
     NL=49;
%      LINES ARE ARRANGED 1-,1+,...37-,37+ IN SPIN-ROTATION SPECTRUM;
%      BY FREQUENCY IN SUBMM SPECTRUM.
     F=[118.7503, 56.2648, 62.4863, 58.4466, 60.3061, 59.5910, ...
         59.1642, 60.4348, 58.3239, 61.1506, 57.6125, 61.8002, ...
         56.9682, 62.4112, 56.3634, 62.9980, 55.7838, 63.5685, ...
         55.2214, 64.1278, 54.6712, 64.6789, 54.1300, 65.2241, ...
         53.5958, 65.7648, 53.0669, 66.3021, 52.5424, 66.8368, ...
         52.0214, 67.3696, 51.5034, 67.9009, 50.9877, 68.4310, ...
         50.4742, 68.9603, 233.9461, 368.4982, 401.7398, 424.7630, ...
         487.2493, 566.8956, 715.3929, 731.1866, ...
         773.8395, 834.1455, 895.0710];      
     S300=[0.2906E-14,0.7957E-15,0.2444E-14,0.2194E-14, ...
         0.3301E-14,0.3243E-14,0.3664E-14,0.3834E-14, ...
         0.3588E-14,0.3947E-14,0.3179E-14,0.3661E-14, ...
         0.2590E-14,0.3111E-14,0.1954E-14,0.2443E-14, ...
         0.1373E-14,0.1784E-14,0.9013E-15,0.1217E-14, ...
         0.5545E-15,0.7766E-15,0.3201E-15,0.4651E-15, ...
         0.1738E-15,0.2619E-15,0.8880E-16,0.1387E-15, ...
         0.4272E-16,0.6923E-16,0.1939E-16,0.3255E-16, ...
         0.8301E-17,0.1445E-16,0.3356E-17,0.6049E-17, ...
         0.1280E-17,0.2394E-17, ...
         0.3287E-16,0.6463E-15,0.1334E-16,0.7049E-14, ...
         0.3011E-14,0.1797E-16,0.1826E-14,0.2193E-16, ...
         0.1153E-13,0.3974E-14,0.2512E-16];
      BE=[.010, .014, .083, 0.083, .207, 0.207, .387, .387, .621,.621, ...
         .910, .910,1.255,1.255,1.654,1.654,2.109,2.109,2.618,2.618, ...
         3.182,3.182,3.800,3.800,4.474,4.474,5.201,5.201,5.983,5.983, ...
         6.819,6.819,7.709,7.709,8.653,8.653,9.651,9.651, ...
         .019, .048, .045, .044, .049, .084, .145, .136, .141, .145, .201];
%     WIDTHS IN MHZ/MB
      WB300=.56;
      X=.8;
      W300=[1.688, 1.703, 1.513, 1.491, 1.415, 1.408, ...
         1.353, 1.339, 1.295, 1.292, 1.262, 1.263, 1.223, 1.217, ...
         1.189, 1.174, 1.134, 1.134, 1.089, 1.088, 1.037,1.038, ...
         0.996,0.996,0.955,0.955,0.906,0.906,0.858,0.858,0.811,0.811, ...
         0.764,0.764,0.717, 0.717,0.669,0.669, ...
         1.65,1.64,1.64,1.64,1.60,1.60,1.60,1.60,1.62,1.47,1.47];
      Y300=[-0.0360, 0.2547, -0.3655,  0.5495, ...
         -0.5696,  0.6181, -0.4252,  0.3517, -0.1496,  0.0430, ...
         0.0640, -0.1605,  0.2906, -0.3730,  0.4169, -0.4819, ...
         0.4963, -0.5481,  0.5512, -0.5931,  0.6212, -0.6558, ...
         0.6920, -0.7208,  0.7312, -0.7550,  0.7555, -0.7751, ...
         0.7914, -0.8073,  0.8307, -0.8431,  0.8676, -0.8761, ...
         0.9046, -0.9092,  0.9416, -0.9423, ...
         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.];
      V=[0.0079, -0.0978,  0.0844, -0.1273, ...
         0.0699, -0.0776,  0.2309, -0.2825,  0.0436, -0.0584, ...
         0.6056, -0.6619,  0.6451, -0.6759,  0.6547, -0.6675, ...
         0.6135, -0.6139,  0.2952, -0.2895,  0.2654, -0.2590, ...
         0.3750, -0.3680,  0.5085, -0.5002,  0.6206, -0.6091, ...
         0.6526, -0.6393,  0.6640, -0.6475,  0.6729, -0.6545, ...
         0.680,  -0.660,   0.685,  -0.665,  ...
         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.];

      TH = 300./TEMP;
      PRESWV = VAPDEN*TEMP/216.68;
      PRESDA = PRES -PRESWV;
      DEN = .001*(PRESDA*(TH^X) + 1.2*PRESWV*TH);
      SUM = 1.584E-17*FREQ*FREQ*WB300*DEN/(TH*(FREQ*FREQ + (WB300*DEN).^2));

      DF = W300*DEN;
      Y = DEN*(Y300+V.*(TH-1));

      SUM = SUM + sum( (S300.*exp(-BE*(TH-1))) .* ( ((DF + (FREQ-F).*Y)./((FREQ-F).^2 + DF.*DF)) + ((DF - (FREQ+F).*Y)./((FREQ+F).^2 + DF.*DF)) ) .* (FREQ./F).^2 );

      absor = 1.6097E11*SUM*PRESDA*TH^3;
      absor = max(absor,0.);


return





%==========================================================================

function xtbscat  = tb_sea_scattering_component( freq, tht, tc, ows, tran )

%   freq   frequency                                       GHz      [1.4, 89.0]    
%   tht    earth incidence angle                           deg        [0,65]
%   ows    ocean wind speed                                m/s        [0,40]          
%   tran   atmospheric trnasmittance                                   [0,1]    
%   tbdw   downwelling atmospheric brightness temperature  Kelvin       >=0                                                  optional
%   tc     cold space temperature                          Kelvin       >=0    
%   etot   sea surface emissivity  (1)=v, (2)=h                        [0,1]    


    costht    = cosd(tht);
    path      = 1.00035/sqrt(costht*costht+7.001225e-4); 
    opacty    = -log(tran)/path;
    xscat     = fd_scatterm_all(freq,tht,ows,opacty);
    xtbscat(1)  = xscat(1,1,1,1,1);
    xtbscat(2)  = xscat(1,1,1,1,2);

return



%==========================================================================

function xtbscat  = tb_sea_scattering_playing( freq, tht, tc, tbdw, ows, tran, etot )

%   freq   frequency                                       GHz      [1.4, 89.0]    
%   tht    earth incidence angle                           deg        [0,65]
%   ows    ocean wind speed                                m/s        [0,40]          
%   tran   atmospheric trnasmittance                                   [0,1]    
%   tbdw   downwelling atmospheric brightness temperature  Kelvin       >=0                                                  optional
%   tc     cold space temperature                          Kelvin       >=0    
%   etot   sea surface emissivity  (1)=v, (2)=h                        [0,1]    
     
    costht    = cosd(tht);
    path      = 1.00035/sqrt(costht*costht+7.001225e-4); 
    opacty    = -log(tran)/path;
    xscat     = fd_scatterm_all(freq,tht,ows,opacty);

    xtbscat(1)=xscat(1,1,1,1,1);
    xtbscat(2)=xscat(1,1,1,1,2);

    if 0 
    xomega(1:2)  = xscat(1:2) /(tbdw + tran*tc - tc); %! [MW 2012] eq. (21)
    xrtot(1:2)   = 1.0-etot(1:2);
    %xtbscat(1:2) = xomega(1:2).*(tbdw+tran*tc) - xomega(1:2).*tc;
    %xtbscat(1:2) = xomega(1:2).*tbdw + xomega(1:2) * tran*tc - xomega(1:2).*tc;
    xtbscat(1:2) = xomega(1:2).*tbdw + xomega(1:2) *tc * ( tran -1 );
    end


return




%==========================================================================

function xtbscat  = tb_sea_scattering_component_orig( freq, tht, tc, tbdw, ows, tran, etot )

%   freq   frequency                                       GHz      [1.4, 89.0]    
%   tht    earth incidence angle                           deg        [0,65]
%   ows    ocean wind speed                                m/s        [0,40]          
%   tran   atmospheric trnasmittance                                   [0,1]    
%   tbdw   downwelling atmospheric brightness temperature  Kelvin       >=0                                                  optional
%   tc     cold space temperature                          Kelvin       >=0    
%   etot   sea surface emissivity  (1)=v, (2)=h                        [0,1]    


    costht=cosd(tht);
    path=1.00035/sqrt(costht*costht+7.001225e-4); 
    opacty=-log(tran)/path;
    xscat=fd_scatterm_all(freq,tht,ows,opacty);
    xomega(1:2) = xscat(1:2) /(tbdw + tran*tc - tc); %! [MW 2012] eq. (21)
    xrtot(1:2) = 1.0-etot(1:2);
    xtbscat(1:2) = ((1.0+xomega(1:2)).*(tbdw+tran*tc) - xomega(1:2).*tc).*xrtot(1:2);
    %xtbscat(1:2) = ( tbdw + tran*tc + [xscat(1,1,1,1,1) xscat(1,1,1,1,2)] ).*xrtot(1:2);


return




%==========================================================================

function [ xscat ] = fd_scatterm_all(freq,tht,ows,opacty)
%function [ xscat ] = fd_scatterm_all(freq,tht,ows,opacty,filetable)
%   freq   frequency                 GHz      [1.4, 89.0]    
%   tht    earth incidence angle     deg        [0,65]
%   ows    ocean wind speed          m/s        [0,40]  
%   opacty opacity                               >=0
% [MW 2012], section V.


    global scatterm

    %= original
    %fileID=fopen('mk_scatterm_table_all.dat');
    %data=fread(fileID,Inf,'real*4');
    %scatterm=reshape(data,91,50,26,13,2);
    %= .mat saved 

    xlog_freq=log10(freq);
    
%   multi-linear interpolation from table values
    
    brief=tht;
    if(brief >= 89.99)
        brief=89.99;
    end
    i1=floor(1+brief);
    i2=i1+1;
    a1=i1-brief;
    a2=1.-a1;
    
    brief=ows;
    if(brief >= 24.99) 
        brief=24.99;
    end
    j1=floor(1+brief);
    j2=j1+1;
    b1=j1-brief;
    b2=1-b1;
    
    brief=xlog_freq/0.2;
    if(brief >= 11.99) 
        brief=11.99;
    end
    k1=floor(1+brief);
    k2=k1+1;
    c1=k1-brief;
    c2=1-c1;
    
    brief=opacty/0.025;
    if(brief >= 48.99) 
        brief=48.99;
    end
    l1=floor(1+brief);
    l2=l1+1;
    d1=l1-brief;
    d2=1-d1;

    xscat1= a1.*b1.*(c1.*scatterm(i1,l1,j1,k1,:)+c2.*scatterm(i1,l1,j1,k2,:))+a1.*b2.*(c1.*scatterm(i1,l1,j2,k1,:)+c2.*scatterm(i1,l1,j2,k2,:))+a2.*b1.*(c1.*scatterm(i2,l1,j1,k1,:)+c2.*scatterm(i2,l1,j1,k2,:))+a2.*b2.*(c1.*scatterm(i2,l1,j2,k1,:)+c2.*scatterm(i2,l1,j2,k2,:));

    xscat2= a1.*b1.*(c1.*scatterm(i1,l2,j1,k1,:)+c2.*scatterm(i1,l2,j1,k2,:))+a1.*b2.*(c1.*scatterm(i1,l2,j2,k1,:)+c2.*scatterm(i1,l2,j2,k2,:))+ a2.*b1.*(c1.*scatterm(i2,l2,j1,k1,:)+c2.*scatterm(i2,l2,j1,k2,:))+a2.*b2.*(c1.*scatterm(i2,l2,j2,k1,:)+c2.*scatterm(i2,l2,j2,k2,:));
    
    xscat=d1.*xscat1 + d2.*xscat2;

return



%==========================================================================

% Z2P_SIMPLE   Simple conversion from altitudes to pressures
%
%    Altitudes are converted to pressures by assuming that the pressure
%    at 0 m is 1000 hPa, and that the pressure drops with a factor of 10
%    for each 16 km increase in altitude.
%
% FORMAT   p = z2p_simple( z )
%        
% OUT   p   Pressure [Pa]
% IN    z   Altitudes [m]

% 2004-09-26   Created by Patrick Eriksson.
%              from atmlab package http://www.radiativetransfer.org/tools/

function p = z2p_simple( z )

p = 10.^( 5 - z/16e3 );

return




%==========================================================================

% P2Z_SIMPLE   Simple conversion from pressures to altitudes
%
%    Pressures are converted to altitudes by assuming that the pressure
%    at 0 m is 1000 hPa, and that the pressure drops with a factor of 10
%    for each 16 km increase in altitude.
%
% FORMAT   z = p2z_simple( p )
%        
% OUT   z   Altitudes [Km]
% IN   p   Pressure [Pa]

% 2005-05-11   Created by Patrick Eriksson.
%              from atmlab package http://www.radiativetransfer.org/tools/

function z = p2z_simple( p )

z = 16 * ( 5 - log10(p) );

return


%==========================================================================

% P2Z_BAROMETRIC   Conversion from pressures to relative altitudes
%
%    Pressures are converted to relative altitudes 
%    using the barometric formula. 
%
% FORMAT   z = p2z_vbarometric( t, p, zo )
%        
% OUT   z   Relative altitudes [Km]
% IN    p   Pressure [Pa]
%       t   Temperature [K]
%       zo  Reference altitude of p(1) [Km]

function z = p2z_barometric( p, t, zo )

Mdry = 28.9644;     % molar mass of dry air: g.mol-1
R    = 8.3144598;   % universal gas constant: J/(molÂ·K)
g    = 9.80665;     % gravitational acceleration: m/s2

% assuming constant gravitational acceleration and molar mass
z    = R.*(t-t(1))./Mdry./g.*log(p./p(1))./log(t(1)./t); 


% reference altitude
z(1) = zo;

return



%==========================================================================

function emis = emis_land_surface_table_at_005km( lon, lat, month, dir_input, freq, angle )
% at 005 deg resolution

edir	= [ dir_input, '/EmisLand/Res_005/' ];
efile	= [ edir, '/land_emis_hd_', sprintf('%02.0f',freq), 'GHz_', sprintf('%0.0f',angle), 'deg_', sprintf('%02.0f',month), '.mat'];
emis  = load_data_single( efile );

efile	= [ edir, '/land_emis_hd_lon_lat.mat' ];
lonlat  = load_data_single( efile );



%= lon to 0-360

if lon < 0
  lon = lon + 360;
end


ilon  = lonlat(1:360*20,1);
ilo   = nearest_in_vector( ilon, lon );
ilat  = lonlat(ilo:360*20:end,2);
ila   = nearest_in_vector( ilat, lat );


%= the climatological value

%= lon-lat matching    

ind   = ilo+ila*360*20; 

%= month and location selection

emisv  = emis( ind, 1 );
emish  = emis( ind, 2 );



%= if emis is nan it means we are at sea, so there is a
%  coast resolution inconsistency problem; we take the closest
%  land at our ~25 km resolution


if isnan(emisv) | isnan(emish)

  jlon = [-1 0 1];
  jlat = [-1 0 1];

  gemisv = zeros(9,1);
  gemish = zeros(9,1);
  glon	 = zeros(9,1);
  glat	 = zeros(9,1);
  dist	 = zeros(9,1);
 
  co     = 0;

  for alo = 1:3
    for ala = 1:3

      co = co+1;

      glon(co) = lon + jlon(alo);
      ilo      = nearest_in_vector( ilon, glon(co) );

      ilat     = emis(ilo:360*4:end,2);
      glat(co) = lat + jlat(ala);
      ila      = nearest_in_vector( ilat, glat(co) );

      idist(co) = sqrt( (glon(co)-lon).^2 + (glat(co)-lat).^2 );

      ind      = ilo+ila*360*4;   

      gemisv(co)  = emis( ind, 1 );
      gemish(co)  = emis( ind, 2 );


    end
  end 

  % sorting by distance

  [ ~, ind ] = sort( idist );

  emisv = -999;
  emish = -999;
 
  for d = ind

    if gemisv(d) >0 & gemish(d) > 0
      emisv = gemisv(d);
      emish = gemish(d);
      break
    end

  end

end

emis = [emisv emish];


return




%==========================================================================


function emis = emis_land_surface_table_angle( lon, lat, month, dir_input, freq, angle )
% at 004 deg resolution

edir	= [ dir_input, '/EmisLand/Res_004/' ];
efile	= [ edir, '/land_emis_hdfillint_', sprintf('%02.0f',freq), 'GHz_', sprintf('%0.0f',angle), 'deg_', sprintf('%02.0f',month), '.mat'];


flon    = 0.018:0.036:359.982;
flat    = -89.982:0.036:89.982;

m       = matfile(efile); 
%load(efile)

%= lon to 0-360

if lon < 0
  lon = lon + 360;
end

[ ~, ilo ] = min( abs( flon - lon ) );
[ ~, ila ] = min( abs( flat - lat ) );


emisv  = double(m.data( ilo, ila, 1 ))/1e3;
emish  = double(m.data( ilo, ila, 2 ))/1e3;
%emisv  = double(data( ilo, ila, 1 ))/1e3;
%emish  = double(data( ilo, ila, 2 ))/1e3;
emis   = [emisv emish];


if emis(1) == 0 | emis(2) == 0

  %= looking for a closer one within +-10

  io = [ilo-10:ilo+10];
  io = io( io>=1 & io <= 10000);
  ia = [ila-10:ila+10];
  ia = ia( ia>=1 & ia <= 5000);
  aux  = double(m.data( io, ia, 1 ))/1e3;   
  %aux  = double(data( io, ia, 1 ))/1e3;   
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

emis = [ emis -999 -999];



return






%==========================================================================


function emis = emis_land_surface_table( lon, lat, month, dir_input, freq, angle )
% at 004 deg resolution

edir	= [ dir_input, '/EmisLand/Res_004/' ];
efile	= [ edir, '/land_emis_hdfillint_', sprintf('%02.0f',freq), 'GHz_', sprintf('%0.0f',angle), 'deg_', sprintf('%02.0f',month), '.mat'];


flon    = 0.018:0.036:359.982;
flat    = -89.982:0.036:89.982;

m       = matfile(efile); 
%load(efile)

%= lon to 0-360

if lon < 0
  lon = lon + 360;
end

[ ~, ilo ] = min( abs( flon - lon ) );
[ ~, ila ] = min( abs( flat - lat ) );


emisv  = double(m.data( ilo, ila, 1 ))/1e3;
emish  = double(m.data( ilo, ila, 2 ))/1e3;
%emisv  = double(data( ilo, ila, 1 ))/1e3;
%emish  = double(data( ilo, ila, 2 ))/1e3;
emis   = [emisv emish];


if emis(1) == 0 | emis(2) == 0

  %= looking for a closer one within +-10

  io = [ilo-10:ilo+10];
  io = io( io>=1 & io <= 10000);
  ia = [ila-10:ila+10];
  ia = ia( ia>=1 & ia <= 5000);
  aux  = double(m.data( io, ia, 1 ))/1e3;   
  %aux  = double(data( io, ia, 1 ))/1e3;   
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

emis = [ emis -999 -999];



return




%==========================================================================

function emis = emis_sea_surface_interp( sst, sss, ows,  dir_input, freq, angle );



%= table organized as OWS, SST, and SSS
%       OWS= 0:1:40
%       SST= 270:1:305
%       SSS= 0:1:40

edir	= [ dir_input, '/EmisSea' ];
efile	= [ edir, '/sea_emis_', sprintf('%02.0f',freq), 'GHz_', sprintf('%0.0f',angle), 'deg.mat'];
table   = load( efile, 'data','ows','sst','sss' );


%= linear interpolation


% V-pol
emis(1) = interpn( table.ows, table.sst, table.sss, squeeze(table.data(:,:,:,1)), ows, sst, sss, 'linear' ); 


% H-pol
emis(2) = interpn( table.ows, table.sst, table.sss, squeeze(table.data(:,:,:,2)), ows, sst, sss, 'linear' ); 

return



%==========================================================================

function emis = emis_sea_surface_table( sst, sss, ows, owd, freq, angle);


%= table organized as OWS, SST, SSS, and OWS
iows= 0:0.1:40;
isst= 270:0.1:305;
isss= 0:0.1:40;
iowd= 0:18:360;


eval(['global ', 'sea_emis_surfem_ocean_', sprintf('%02.0f', freq), 'GHz_', sprintf('%0.0f',angle) ]);
eval( [' data = sea_emis_surfem_ocean_', sprintf('%02.0f', freq), 'GHz_', sprintf('%0.0f',angle), ';'])


io = 1+round(10*ows);
it = 1+ round(10*(sst-270));
is = 1+round(10*sss);
id = 0*owd;
for a = 1:length(owd)
  [~, id(a)] = min(abs(owd(a)-iowd));
end 


if io>401
  io=401;
end
if io<1
  io=1;
end
if it>351
  it=351;
end
if it<1
  it = 1;
end
if is>401
  is=401;
end
if is<1
  is=1;
end

if 0
emisv  = double(data( io, it, is, 1 ))/1e3;
emish  = double(data( io, it, is, 2 ))/1e3;
emis   = [emisv emish];



if emis(1) == 0 | emis(2) == 0
  emis(1) = -999;
  emis(2) = -999;
end
end

emis = squeeze(data( io, it, is, id, : ));


return



%==========================================================================


function [ emis, diagcov ]  = emis_ice_surface_table( lat, xyi, t2m, snd, dir_input, freq, angle)


%= table organized as T2M, SND, XYI, LAT
%it2m = 240:0.1:270;
%isnd = [0:0.01:0.6];
%ilat = [-45 45]; 
%ixyi = [1 -1];   %FYI MYI



%= hack using 55deg emis at all freqs

if freq ~= 1
  angle = 55;
else
  if angle < 35
    angle = 35;
  elseif angle> 55
    angle = 55;
  end
end



%= loading data from global variable

eval(['global ', 'ice_emis_hd_', sprintf('%02.0f', freq), 'GHz_', sprintf('%0.0f',angle) ]);
eval( [' data = ice_emis_hd_', sprintf('%02.0f', freq), 'GHz_', sprintf('%0.0f',angle), ';']);



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

emis    = [ emis' -999 -999];
diagcov = [nan nan 0 0];


return





%==========================================================================


function [ emis, diagcov ]  = emis_ice_surface_table_xyi( lat, xyi, t2m, sit, freq, angle)



%= table organized as T2M, SIT, XYI, LAT
it2m = 240:0.1:270;
isit = [0:0.01:0.6];
ilat = [-45 45]; 
ixyi = [1 -1];   %FYI MYI


eval(['global ', 'ice_emis_hd_', sprintf('%02.0f', freq), 'GHz_', sprintf('%0.0f',angle) ]);
eval( [' data = ice_emis_hd_', sprintf('%02.0f', freq), 'GHz_', sprintf('%0.0f',angle), ';'])


it = 1+ round(10*(t2m-240));
is = 1+ round(100*sit); 

%= is limited to sit=0.6 as there was
%  no larger SIT when constructing the table
%  with the RRDP
if is>61
  is=61;
end

%= it limited to sit=270, same reasons
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

emis  = squeeze(double(data( it, is, :, il, 1 ))/1e3);
emisv = emis(1) * (1-xyi) + emis(2) * xyi;

emis  = squeeze(double(data( it, is, :, il, 2 ))/1e3);
emish = emis(1) * (1-xyi) + emis(2) * xyi;

emis   = [emisv; emish];

if isnan(emis(1)+emis(2))
  emis(1) = -999;
  emis(2) = -999;
end


diagcov = [nan;nan];


return




%==========================================================================


function [ emis, diagcov ]  = emis_ice_surface_table_no_xyi_interp( lat, xyi, t2m, sit, dir_input, freq, angle)



%= table organized as T2M, SIT, XYI, LAT
it2m = 240:0.1:270;
isit = [0:0.01:0.6];
ilat = [-45 45]; 
ixyi = [1 -1];   %FYI MYI


eval(['global ', 'ice_emis_hd_', sprintf('%02.0f', freq), 'GHz_', sprintf('%0.0f',angle) ]);
eval( [' data = ice_emis_hd_', sprintf('%02.0f', freq), 'GHz_', sprintf('%0.0f',angle), ';'])


if 0
edir	= [ dir_input, '/EmisIce' ];
efile	= [ edir, '/ice_emis_hd_', sprintf('%02.0f',freq), 'GHz_', sprintf('%0.0f',angle), 'deg.mat'];
m       = matfile(efile); 
end


it = 1+ round(10*(t2m-240));
is = 1+ round(100*sit); 

%= is limited to sit=0.6 as there was
%  no larger SIT when constructing the table
%  with the RRDP
if is>61
  is=61;
end

%= it limited to sit=270, same reasons
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
if xyi > -0.025
  ix = 1; %FYI
else
  ix = 2; %MYI
end

%emisv  = double(m.data( it, is, ix, il, 1 ))/1e3;
%emish  = double(m.data( it, is, ix, il, 2 ))/1e3;
emisv  = double(data( it, is, ix, il, 1 ))/1e3;
emish  = double(data( it, is, ix, il, 2 ))/1e3;
emis   = [emisv; emish];

if isnan(emis(1)+emis(2))
  emis(1) = -999;
  emis(2) = -999;
end



diagcov = [nan;nan];


return


%==========================================================================


function [ emis, diagcov ]  = emis_ice_surface_nointerp( lat, xyi, t2m, sit, dir_input, freq, do_cov )


global DATADIR_MAIN

edir	= [ dir_input, '/EmisIce' ];
efile	= [ edir, '/ice_emis_SICCI-SMOS.mat' ];      
load( efile, 'str' );


%=== frequencies stored for look up table
%    selecting the index in emis

ifreq   = [6.925 10.65 18.7 36.5 89.0 1.4];
id_emis = nearest_in_vector( ifreq, freq);
id_emis = [ 2 * (id_emis -1) + 1 2 * id_emis ];


%=== finding entry in look up table

it2m = find_in_range( str.index_t2m, t2m);
isit = find_in_range( str.index_sit, sit);
ilat = find_in_range( str.index_lat, lat);
ixyi = find_in_range( str.index_xyi, xyi);




%=== getting mean value and covariance

mea_emis      = squeeze( str.mea_emis( ilat, ixyi, it2m, isit,  :) );
min_emis      = squeeze( str.min_emis( ilat, ixyi, it2m, isit,  :) );
max_emis      = squeeze( str.max_emis( ilat, ixyi, it2m, isit,  :) );
cov_emis_A    = squeeze( str.cov_emis_A( ilat, ixyi, it2m, isit, :, :) );
cov_emis_S    = squeeze( str.cov_emis_S( ilat, ixyi, it2m, isit, :, :) );



%=== if a case is found in the table

if ~isnan( sum( mea_emis) )

  %=== generating channels with mean and covariance matrix

  aux = size( cov_emis_A, 1 );

  if max(id_emis) <= 10  

    if do_cov

      %= AMSRE
      na  = 1:aux;
      emis = rnd_emis( mea_emis(na), cov_emis_A, min_emis(na), max_emis(na) );

      %= selecting for the given frequency
      emis    = emis( id_emis );

    else

      emis = mea_emis( id_emis );

    end

    diagcov = diag(cov_emis_A(id_emis,id_emis));

    else

      if do_cov

        %= SMOS
        ns =  aux+1:length(mea_emis);
        emis = rnd_emis( mea_emis(ns), cov_emis_S, min_emis(ns), max_emis(ns) );

      else
 
        emis = mea_emis( id_emis );

      end
        
      diagcov = diag(cov_emis_S);

    end


  %= reversing to produce V-pol and H-pol
  emis = emis([2 1]);

else

  %=== no case found, no extrapolation so far
  emis = [];

end


%=== limiting emis to one removed as emis > 1 
%    likely related to emitting temperature 
%    warmer than ice skin temperature. As 
%    emis will be again mutiplied by skin
%    temperature to derive BTs we need this
%    "effective emissivity" even if it looks
%    unphysical.



return


%==========================================================================


function [ emis, diagcov ]  = emis_ice_surface_interp( lat, xyi, t2m, sit, dir_input, freq, ang, do_cov )


global DATADIR_MAIN

edir	= [ dir_input, '/EmisIce' ];
efile	= [ edir, '/ice_emis_SICCI-SMOS.mat' ];      
load( efile, 'str' );


%=== frequencies stored for look up table
%    selecting the index in emis

ifreq   = [6.925 10.65 18.7 36.5 89.0 1.4];
nf      = 2+2*length(ifreq);
id_emis = nearest_in_vector( ifreq, freq);


if id_emis ~= 6

  id_emis = [ 2 * (id_emis -1) + 1   2 * id_emis ];

else

  %=== only for 1.4 there is angle selection
  %    for the other freqs only estimates at
  %    55 deg are used independent of angle
 
  if ang > 45
    id_emis = [ 2 * (id_emis -1) + 1   2 * id_emis ];    
  else
    id_emis = [ 2 * (id_emis) + 1   2 * (id_emis+1) ];    
  end

end 



%=== finding entry in look up table

% temp limited in table to 240 and 27e
t2m(t2m<240) = 240.1;
t2m(t2m>273) = 272.9;
t2m = double(t2m);

% sit limited to 0.1 and 0.6
sit(sit<0.1) = 0.11;
sit(sit>0.6) = 0.59;
sit = double(sit);

it2m = find_in_range( str.index_t2m, t2m);
isit = find_in_range( str.index_sit, sit);
ilat = find_in_range( str.index_lat, lat);
ixyi = find_in_range( str.index_xyi, xyi);


%=== getting mean value and covariance

if 0

  % no interp 
  mea_emis      = squeeze( str.mea_emis( ilat, ixyi, it2m, isit,  :) );

else

  % interp in T2M and SIT
  tab_emis      = squeeze( str.mea_emis( ilat, ixyi, :, :,  :) );
  mea_emis      = nan(nf,1);
  for f = 1:nf
    mea_emis(f) = interpn( str.mid_t2m, str.mid_sit, squeeze(tab_emis(:,:,f)), t2m, sit, 'linear' ); 
  end

end

min_emis      = squeeze( str.min_emis( ilat, ixyi, it2m, isit,  :) );
max_emis      = squeeze( str.max_emis( ilat, ixyi, it2m, isit,  :) );
cov_emis_A    = squeeze( str.cov_emis_A( ilat, ixyi, it2m, isit, :, :) );
cov_emis_S    = squeeze( str.cov_emis_S( ilat, ixyi, it2m, isit, :, :) );



%=== if a case is found in the table

if ~isnan( sum( mea_emis) )

  %=== generating channels with mean and covariance matrix

  aux = size( cov_emis_A, 1 );

  if max(id_emis) <= 10  

    if do_cov

      %= AMSRE
      na  = 1:aux;
      emis = rnd_emis( mea_emis(na), cov_emis_A, min_emis(na), max_emis(na) );

      %= selecting for the given frequency
      emis    = emis( id_emis );

    else

      emis = mea_emis( id_emis );

    end

    diagcov = diag(cov_emis_A(id_emis,id_emis));

    else

      if do_cov

        %= SMOS
        ns =  aux+1:length(mea_emis);
        emis = rnd_emis( mea_emis(ns), cov_emis_S, min_emis(ns), max_emis(ns) );

      else
 
        emis = mea_emis( id_emis );

      end
        
      diagcov = diag(cov_emis_S);

    end


  %= reversing to produce V-pol and H-pol
  emis = emis([2 1]);

else

  %=== no case found, no extrapolation so far
  emis = [];

end




%=== limiting emis to one removed as emis > 1 
%    likely related to emitting temperature 
%    warmer than ice skin temperature. As 
%    emis will be again mutiplied by skin
%    temperature to derive BTs we need this
%    "effective emissivity" even if it looks
%    unphysical.


return




%==========================================================================


function ind = find_in_range(z,zi);

[ zn, ind ] = min( abs( z -zi ) );
zd           = zi - z(ind); 
if zd < 0
  ind = ind-1;
end
ind = ind+1;

return




%==========================================================================

function  emis = rnd_emis( mea_emis, cov_emis, min_emis, max_emis )


%= choleski decomposition
[R,p] = chol(cov_emis);

%= making cov matrix positive definite

if p ~= 0

  [V,D] = eig( cov_emis );       % Calculate the eigendecomposition of your matrix (A = V*D*V') 
  	                         % where "D" is a diagonal matrix holding the eigenvalues of your matrix "A"
   d= diag(D);       		 % Get the eigenvalues in a vector "d" 
   d(d <= 1e-7) = 1e-7;  	 % Set any eigenvalues that are lower than threshold "TH" ("TH" here being 
                        	 % equal to 1e-7) to a fixed non-zero "small" value (here assumed equal to 1e-7)
   D_c = diag(d);        	 % Built the "corrected" diagonal matrix "D_c"
   cov_emis = V*D_c*V';      		 % Recalculate your matrix "A" in its PD variant "A_PD"

end


%=== rnd generation of 1000 cases

n = size(cov_emis,1);
r = 10000;
X = randn( r, n );
X = X - mean(X);
X = X * inv(chol(cov(X)));
X = X * chol(cov_emis); 


mea  = repmat(mea_emis', r, 1);
mea  = mea + X;


%=== outputting first case where min and max values
%   are not overpassed

emis = '';

for f = 1:r

  co(f)  = sum( mea(f,:) > min_emis' ) +  sum( mea(f,:) < max_emis' );
  if sum( mea(f,:) > min_emis' & mea(f,:) < max_emis' ) == n

    emis = mea(f,:); 
    break

  end  
  
end



if isempty(emis)

  [ dum, ind] = max(co);
  disp('generated emis not within min max, mean as output')
  emis = mea(ind,:);

end



return


% SUB MODULE nearest_in_vector
%
%          Find the nearest value in a vector z to a given value zi
%
% FORMAT   [ind,zn,zd] = nearest_in_vec(z,zi)
%
% OUT      ind  the index in vector z of closest value to zi
%          zn   the closest value in z to zi
%          zd   the distance between zn and zi
% IN       z    the vector
%          zi   the value

% 2019-01-08   Created by carlos.jimenez@obspm.fr

function [ind,zn,zd] = nearest_in_vector(z,zi);

[ aux, ind ] = min( abs( z -zi ) );
zn           = z( ind );
zd           = zi - zn; 

return






