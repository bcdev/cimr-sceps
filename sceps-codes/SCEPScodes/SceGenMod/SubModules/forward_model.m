%-------------------------------------------------------------------------------
%
% SUBMODULE   forward_model
%
%    A function to simulate vertically and horizontally polarized 
%    top-of-atmosphere brightness temperatures (BTs) at a specified set 
%    of frequencies and polarizations for a given observing zenith 
%    angle, atmosphere and underlaying ice-free, sea, or coastal 
%    surface. 
%
% FORMAT   [R,y,J] = forward_model( Q, R, x, iter )
%
% OUT   R		structure 		Retrieval data structure. Identical
%						to the input R, just here to comply
%						with expected format when used by
%						the inversion function oem.
%
%	y               vector			BTs column vector, length the number
%						of frequencies x number of polarizations.
%						E.g., [ freq1-Vpol freq1-Hpol freq2-Hpol]
%						for Q.sensor_input.F = [freq1 freq2], 
%						Q.sensor_input.Vpol = [1 1], and
%						Q.sensor_input.Hpol = [0 1].
%
%       J	        matrix			Weighting functions (Jacobians) matrix, 
%						size (length of y) x (length of x).
%
%												
% IN    Q		Structure		Nested structure containing the atmosphere, 
%						surface, and sensor  for the forward
%						modelling:
%
% 				atmos_input	structure array containing the
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
%			   CWVC	[Kg/m2]		Column Water Vapour Content 
%			   CLWC	[Kg/m2]		Cloud Liquid Water Content	
%						
%
%			   
%				surf_input	structure containing 
%						the surface parameters for
%						a sea simulation:			
%
%			    SST	[K]		Surface Temperature
%			    SSS [g/Kg]		Sea Surface Salinity
%			    OWS	[m/s]		Wind velocity
%			    UWS	[m/s]		U wind component
%			    VWS	[m/s]		V wind component
%			    SLP[mbar]		Mean Sea Level Pressure
%	
%						and/or the surface parameters
%						for a land simulation:
%
%			    LST	[K]		Land Surface Temperature
%			    LAT [degrees]	Latitude, 90S to 90N
%			    LON [degrees]	Longitude, 0 to 360
%			    MONTH []		Month of the year
%			    LLP  [mbar]         Land Level Pressure 
%			    LSP [0-1]		Land Surface Percentage,
%						optional, required for a
%						mixed land and sea scene
%
%						and/or the surface parameters
%						for an ice simulation:
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
%			    LIP [0-1]		Land Surface Percentage,
%						optional, required for a
%						mixed land and ice scene
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
%				sensor_input	structure containing 
%						the sensor parameters:
%
%			    F	 [GHz]		Frequency		
%			    VPOL []		Vector the size of F indicating
%						if V-polarization for the corresponding
%						freq is calculated (1) or not (0).	
%			    VPOL []		Vector the size of F indicating
%						if H-polarization for the corresponding
%						freq is calculated (1) or not (0).	
%			    ZA	 [degrees]	Viewing zenith angle
%			    AA	 [degrees]	Viewing azimuth angle
%			    H	 [m]		Sensor height from ground
%			    TOA  []		Boolean indicating if the sensor radiative
%						transfer should be done only at the surface
%						(0) or contain also the atmospheric 
%						contributions (1).

%			    	
%	
%				dir_input	string indicating folder where the
%						following auxiliary files may to be
%						placed:		         
%			    
%			    land_emis_fXX.mat	Land surface emissivity
%						for frequency XX, placed
%						at dir_input/EmisLand/    
%						
%					J_DO	and a boolean flag indicating
%						if the Jacobians need to be 
%						calculated:
%
%						true for jacobian calculation
%						false for no jacobian calculation
%						 						       	
% 
%
%       R	structure array			Structure containing retrieval
%						settings to be used with oem.m:
%			    			
%			    R_VARIABLES		cell array with the variables to 
%						be retrieved, e.g.,
%						{'SST'  'SSS'  'OWS'  'CWVC'  'CLWC'}
%			    			or a selection of those ones. It is
%						also possible to retrieve emissivity. 
%						The variables should be named EMIX_FQY
%						where X is V or H depending on the 
%						polarization and Y is a number between
%						1 and the number of frequencies in 
%						sensor_input. For instance
%						{'LST','EMIV_FQ3'} will make oem.m
%						retrieve LST and the V-pol emis of the
%						third given frequency. 
%						
%	
%		            J_PERTURBATION	column vector of length of x 
%						indicating the relative perturbation
%						in percentageto numerically calculate
%						the Jacobians, e.g. [10 15 5 10 5]
%
%	x	vector				Column vector of length the number of 
%						retrieval variables containing their value.
%						NOTICE that the values here overwrite the
%						values passed in the Q structure, but both 
%						values are input to facilitate adapting 
%						the code to the oem inversions. If no
%						jacobians are needed becuase no retrieval is
%						done, i.e., J_DO=0, it can be set with dummy
%						values, or with [], x will not be used.
%
%						
%	iter	scalar				Variable to indicate the linearization
%						point for the Jacobians. So far it has to
%						be set to "1" to signal that the linearization
%						point is the vector x, with relative
%						perturbations as indicated by R.J_PERTURBATION, 
%						with no other possibilities implemented.
%
%-------------------------------------------------------------------------------
% Project:	  CIMR Algorithm Performance Evaluation
% Package:	  CIMR Scientific Work Bench
% Developer:	  Estellus 
% Contact:	  carlos.jimenez@estellus.fr 
% Initiated:	  2019-01-08
%-------------------------------------------------------------------------------



%==== MAIN-CODE =====================================================


function [R,y,y_sea_harm_iso, y_sea_harm_cs1, y_sea_harm_cs2,J] = forward_model( Q, R, x, iter )


%=== inputs checking


if isfield( Q.surf_input, 'SIC' ) & (~isfield( Q.surf_input, 'IST' ) | ~isfield( Q.surf_input, 'SST' ) )
  if Q.surf_input.SIC < 1
    osfi_error('If SIC is an input, we need both IST and SST for the mixed scene')
  elseif ~isfield( Q.surf_input, 'IST' )
    osfi_error('If SIC is an input, we need IST')
  end

end

if isfield( Q.surf_input, 'LSP' ) & (~isfield( Q.surf_input, 'LST' ) | ~isfield( Q.surf_input, 'SST' ) )
  osfi_error('If LSP is an input, we need both LST and SST for the mixed scene')
end



%=== retrieval variables checking
  
aux_var     = {'SST','SSS','OWS','IST','SIC','LST','EMIV_FQ1','EMIH_FQ1','EMIV_FQ2','EMIH_FQ2','EMIV_FQ3','EMIH_FQ3','EMIV_FQ4','EMIH_FQ4','EMIV_FQ5','EMIH_FQ5','CWVC','CLWC'};

%= to modify Q in perturbations, 1 for surface, 2 for atmosphere

if Q.J_DO
  aux_type    = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2];
end

%= variable limits so perturbations
%  do not exceed these values. Ordering
%  as aux_var. Also important to avoid 
%  extrapolation in the emis lookup table

min_per = [270 0   0  150 0 100 0 0 0 0 0 0 0 0 0 0 0   0]; 
max_per = [305 40 50 273 1 370 1 1 1 1 1 1 1 1 1 1 100 100];


nx	    = length( R.R_VARIABLES );
na	    = length( aux_var );
ind_ret     = zeros(nx,1);

for f = 1:nx

  for a = 1:na

    if strcmp( aux_var{a}, R.R_VARIABLES{f} );
      ind_ret(f) = a;
    end 

  end

  if ind_ret(f) == 0
    osfi_error([ R.R_VARIABLES{f}, ' cannot be retrieved for the moment']);
  end

end



if Q.J_DO
  aux_type    = aux_type(ind_ret);
  min_per     = min_per(ind_ret);
  max_per     = max_per(ind_ret);
end




%=== Input values checking
%    NOTE: Some variables have identified ranges that cannot
%    be exceeded. They are contraing here with a message to
%    inform of that. This can be of importance when using the
%    FM for retrievals. For instance, SIC is between 0 and 1.
%    If the inversion module inputs a value outside 0-1, the
%    FM may be operating anomously, e.g., with a negative SIC.
%    Only checked for the retrieval variables as they can be
%    modified by the inversion. For the others no check, so
%    user responsability.
    

if Q.J_DO

  for f = 1:nx
  
    if x(f) < min_per(f)
      x(f) = min_per(f);
      %disp([ R.R_VARIABLES{f}, ' changed to minimum value' ]);
    elseif x(f) > max_per(f)
      x(f) = max_per(f);
      %disp([ R.R_VARIABLES{f}, ' changed to maximum value' ]);
    end

  end

end



%=== overwriting common Q inputs with 
%    surface retrieved state



if Q.J_DO

  if length(x) ~= nx
    osfi_error('R and x are not consistent!');
  end

  for f = 1:nx
    if aux_type(f) == 1
      Q.surf_input = setfield( Q.surf_input, R.R_VARIABLES{f}, x(f) );
    elseif aux_type(f) == 2
      Q.atmos_input = setfield( Q.atmos_input, R.R_VARIABLES{f}, x(f) );
    end
  end

end




%=== forward modelling


[ y, R.emis, y_sea_harm_iso, y_sea_harm_cs1, y_sea_harm_cs2 ] = forward_model_freq( Q.atmos_input, Q.surf_input, Q.dir_input, Q.sensor_input.F, Q.sensor_input.ZA, Q.sensor_input.AA, Q.sensor_input.H, Q.sensor_input.VPOL, Q.sensor_input.HPOL, Q.sensor_input.TOA);


%=== Jacobians 

J  = [];
ny = length(y);






if Q.J_DO

  J = nan( ny, nx );
  
  for  f = 1:nx

    Qj.surf_input  = Q.surf_input;
    Qj.atmos_input = Q.atmos_input;

    %= passing of emissivities from previous RT to be sure that the 
    %  emis do not change. This is only relevant for ice RT as the 
    %  RT code internally generate random emissivities based on an
    %  statistical analysis of the SICCI database, so they change
    %  every time the model is run. For land and sea the emis do
    %  not change internally for a given surface input, but we
    %  pass them to facilitate the coding. Notice that we pass
    %  all emis at all frequencies

    Qj.surf_input.FEMI_ICE = R.emis.ice; 
    if 0
      Qj.surf_input.FEMI_SEA = R.emis.sea; 
      Qj.surf_input.FEMI_LAN = R.emis.lan;   
    end

    %= positive perturbation

    pertu_pos      = x(f) * ( 1 + R.J_PERTURBATION(f)/100 );

    if pertu_pos > max_per(f)
      pertu_pos = max_per(f);
    end

    if aux_type(f) == 1
      Qj.surf_input  = setfield( Qj.surf_input, R.R_VARIABLES{f}, pertu_pos  );
    elseif aux_type(f) == 2
      Qj.atmos_input = setfield( Qj.atmos_input, R.R_VARIABLES{f}, pertu_pos  );
    end

    y_pos = forward_model_freq( Qj.atmos_input, Qj.surf_input, Q.dir_input, Q.sensor_input.F, Q.sensor_input.ZA, Q.sensor_input.AA, Q.sensor_input.H, Q.sensor_input.VPOL, Q.sensor_input.HPOL, Q.sensor_input.TOA);


    %= negative perturbation

    pertu_neg      = x(f) * ( 1 - R.J_PERTURBATION(f)/100 );

    if pertu_neg < min_per(f)
      pertu_neg = min_per(f);
    end

    if aux_type(f) == 1
      Qj.surf_input  = setfield( Qj.surf_input, R.R_VARIABLES{f}, pertu_neg  );
    elseif aux_type(f) == 2
      Qj.atmos_input = setfield( Qj.atmos_input, R.R_VARIABLES{f}, pertu_neg  );
    end

    y_neg = forward_model_freq( Qj.atmos_input, Qj.surf_input, Q.dir_input, Q.sensor_input.F, Q.sensor_input.ZA, Q.sensor_input.AA, Q.sensor_input.H, Q.sensor_input.VPOL, Q.sensor_input.HPOL, Q.sensor_input.TOA);

    %= dy/dx for J matrix

    if sum(y_pos-y_neg) ~= 0
      J(:,f) = (y_pos - y_neg) / (pertu_pos - pertu_neg); 
    else
      J(:,f) = 0;
    end

  end

end



return



%==== SUB-FUNCTION =================================================

function [ y, emis, y_harm_iso, y_harm_cs1, y_harm_cs2 ] = forward_model_freq( atmos_input, surf_input, dir_input, freqs, zangle, aangle, height, do_vpol, do_hpol, do_toa )



nf        = length( freqs );
ny        = sum(do_vpol) + sum(do_hpol)+2*nf;
na        = size(aangle,2);


%y         = nan( ny, na);
%emis.lan  = nan( ny, na);
%emis.sea  = nan( ny, na);
%emis.ice  = nan( ny, na);
%emis.ice_diagcov  = nan( ny, 1);

%= bts for the given eia and aa
y         = -999*ones( ny, na);

%= bts for harmonics calculation

y_harm_iso = y; 
y_harm_cs1 = y; 
y_harm_cs2 = y; 

%= emis
emis.lan  = y;
emis.sea  = y;
emis.ice  = y;
emis.ice_diagcov  = -999*ones( ny, 1);





%=== loop in frequency

co = 0;

sensor_input.H    = height;


for f = 1:nf
  
  sensor_input.F   = freqs(f);
  sensor_input.ZA  = zangle(f);
  sensor_input.AA  = aangle(f);
  sensor_input.TOA = do_toa(f);


  %=== dealing with emis
  %    case 1 ) we are doing a retrieval of emis
  %             we read from the latest EMIV_FQ
  %		at the corresponding frequency
  
  auv = [ 'EMIV_FQ', num2str(f) ];
  auh = [ 'EMIH_FQ', num2str(f) ];


  if isfield( surf_input, auv ) & isfield( surf_input, auh )

    %= the R.EMI values are passed to the 3 possible emis types
    %  for consistency, but forward_model_core will only allow
    %  to be used for single scenes, i.e., only land, sea, or
    %  ice, so the 2 non-concerned emis types are just
    %  dummy values not in use

    surf_input.EMI_SEA = [ getfield( surf_input, auv) getfield( surf_input, auh)];
    surf_input.EMI_LAN = [ getfield( surf_input, auv) getfield( surf_input, auh)];
    surf_input.EMI_ICE = [ getfield( surf_input, auv) getfield( surf_input, auh)];

  elseif isfield( surf_input, 'FEMI_ICE') 

  %=== dealing with emis
  %    case 2 ) we are not retrieving emis but doing the 
  %		the Jacobians so we want to get the 
  %		emissivities used for the previous RT
  %		Only relevant for ice as there are the 
  %		only ones with a noise component so 
  %		far and therefore changing	

    ind = (2*(f-1)+1):2*f;
    if 0
      surf_input.EMI_SEA  = surf_input.FEMI_SEA(ind);
      surf_input.EMI_LAN  = surf_input.FEMI_LAN(ind);
    end
    surf_input.EMI_ICE  = surf_input.FEMI_ICE(ind);

  else

    if isfield( surf_input, 'EMI_LAN') 

       %=== so far used to pass an area averaged emis lan
       %    for the simulations. It may need further thoughts
       %    if doing retrievals when Q.J_DO == 1
       %    All freq emis passed, we leave the v & h
       %    of the corresponding frequency

       if f == 1
         % we save the emis as we overwrite
         lemi = surf_input.EMI_LAN;
       end 
 
       surf_input.EMI_LAN  = [ squeeze(lemi(:,f))' -999 -999];

    end

    if isfield( surf_input, 'EMI_ICE') 

       %=== so far used to pass an area averaged emis lan
       %    for the simulations. It may need further thoughts
       %    if doing retrievals when Q.J_DO == 1
       %    All freq emis passed, we leave the v & h
       %    of the corresponding frequency

       if f == 1
         % we save the emis as we overwrite
         iemi = surf_input.EMI_ICE;
       end 
 
       surf_input.EMI_ICE  = [ squeeze(iemi(:,f))' -999 -999];

    end  
       
  end

  [ bts, emi, bts_harm_iso, bts_harm_cs1, bts_harm_cs2  ]  = forward_model_core( atmos_input, surf_input, sensor_input, dir_input ); 

  %= V component

  if do_vpol(f)
    co        = co+1;

    y(co)     = bts(1);
    y_harm_iso(co) = bts_harm_iso(1); 
    y_harm_cs1(co) = bts_harm_cs1(1); 
    y_harm_cs2(co) = bts_harm_cs2(1); 

    if ~isempty(emi.lan)
      emis.lan(co)  = emi.lan(1);
    end
    if ~isempty(emi.sea)
      emis.sea(co)  = emi.sea(1);
    end
    if ~isempty(emi.ice)
      emis.ice(co)  = emi.ice(1);
      emis.ice_diagcov(co)  = emi.ice_diagcov(1);
    end


  end

  %= H component

  if do_hpol(f)
    co        = co+1;
    y(co)     = bts(2);
    if ~isempty(emi.lan)
      emis.lan(co)  = emi.lan(2);
    end
    if ~isempty(emi.sea)
      emis.sea(co)  = emi.sea(2);
    end
    if ~isempty(emi.ice)
      emis.ice(co)  = emi.ice(2);
      emis.ice_diagcov(co)  = emi.ice_diagcov(2);
    end
    y_harm_iso(co) = bts_harm_iso(2); 
    y_harm_cs1(co) = bts_harm_cs1(2); 
    y_harm_cs2(co) = bts_harm_cs2(2); 

  end

  %= 3rd and 4th stokes

    co = co+1;
    y(co)     = bts(3);
    if ~isempty(emi.lan)
      emis.lan(co)  = emi.lan(3);
    end
    if ~isempty(emi.sea)
      emis.sea(co)  = emi.sea(3);
    end
    if ~isempty(emi.ice)
      emis.ice(co)  = emi.ice(3);
      emis.ice_diagcov(co)  = emi.ice_diagcov(3);
    end
    y_harm_iso(co) = bts_harm_iso(3); 
    y_harm_cs1(co) = bts_harm_cs1(3); 
    y_harm_cs2(co) = bts_harm_cs2(3); 

   co = co+1;
    y(co)     = bts(4);
    if ~isempty(emi.lan)
      emis.lan(co)  = emi.lan(4);
    end
    if ~isempty(emi.sea)
      emis.sea(co)  = emi.sea(4);
    end
    if ~isempty(emi.ice)
      emis.ice(co)  = emi.ice(4);
      emis.ice_diagcov(co)  = emi.ice_diagcov(3);
    end
    y_harm_iso(co) = bts_harm_iso(4); 
    y_harm_cs1(co) = bts_harm_cs1(4); 
    y_harm_cs2(co) = bts_harm_cs2(4); 

  if isfield( surf_input, auv ) & isfield( surf_input, auh )
    surf_input = rmfield( surf_input, 'EMI' );
  end

end


return

