%============================================================================
%
% FORMAT   [R,y,J] = fom_atmos_ocean_seaice_land_tbs_retrieval( Q, R, x, iter )
%
% DESCRIPTION
%
%    A function to simulate vertically and horizontally polarized 
%    top-of-atmosphere brightness temperatures (Tbs) at a specified set 
%    of frequencies and polarizations for a given observing zenith 
%    angle, atmosphere and underlaying surface.
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
% OUT   R		structure 		Retrieval data structure. Identical
%						to the input R, just here to comply
%						with expected format when used by
%						the inversion function oem.
%
%	y               vector			Tbs column vector, length the number
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
% REFERENCES
%
%               This function is implemented as a general forward 
%		model for different retrieval aplications, and it has
%		been developed as part of the Study on CIMR synergistic 
%		global ocean and atmospheric products, contract number
%		ITT 23/224353.
%		
%
% HISTORY
%		2025-03-10    Started by carlos.jimenez@estellus.fr
%		
%		2025-11-10    Minor fixes and adding documentation. 	
%
%==========================================================================


function [R,y,J] = fom_atmos_ocean_seaice_land_tbs_retrieval( Q, R, x, iter )


%=== if not Q.J_DO default is making them

if ~isfield( Q, 'J_DO' )
  % doing jacobians
  Q.J_DO = 1;
end

if ~isfield( Q, 'ISO_DO' )
  % only the isotropic omponent on Tb
  % i.e., the term that it does not depend
  % on relative wind direction
  Q.ISO_DO = 1;
end

%=== retrieval variables checking
  
aux_var     = {'SST','SSS','OWS','TCWV','TCLRW'};

%= to modify Q in perturbations, 1 for surface, 2 for atmosphere

if Q.J_DO

  aux_type    = [1 1 1 2 2];

  %= variable limits so perturbations
  %  do not exceed these values. Ordering
  %  as aux_var. Zeros to be avoided as
  %  x(f) * ( 1 + R.J_PERTURBATION(f)/100 ) and
  %  x(f) * ( 1 + R.J_PERTURBATION(f)/100 ) 
  %  will result in the same x value and same
  %  y value, so the jacobian will be zero. 

  min_value = [270 0.01   0.01  0.01  0.0001]; 
  max_value = [305 40 50 100 100];

end



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
    error([ R.R_VARIABLES{f}, ' cannot be retrieved for the moment']);
  end

end



if Q.J_DO
  aux_type    = aux_type(ind_ret);
  min_value     = min_value(ind_ret);
  max_value     = max_value(ind_ret);
end




%=== Input values checking
%    NOTE: Some variables have identified ranges that cannot
%    be exceeded in the RT calculation (e.g. TCWV < 0 ). 
%    They are contraining here in case the retrieval solution
%    exceeds those ranges.

if Q.J_DO

  for f = 1:nx
  
    if x(f) < min_value(f)
      x(f) = min_value(f);
      %disp([ R.R_VARIABLES{f}, ' changed to minimum value' ]);
    elseif x(f) > max_value(f)
      x(f) = max_value(f);
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

[ y ] = fom_atmos_ocean_seaice_land_tbs_core( Q.atmos_input, Q.surf_input, Q.sensor_input, Q.dir_input, Q.ISO_DO  );


%=== Jacobians 

J  = [];
ny = length(y);



if Q.J_DO

  J = nan( ny, nx );
  
  for  f = 1:nx

    Qj = Q;

    %= positive perturbation

    pertu_pos      = x(f) * ( 1 + R.J_PERTURBATION(f)/100 );

    if pertu_pos > max_value(f)
      pertu_pos = max_value(f);
    end

    if aux_type(f) == 1
      Qj.surf_input  = setfield( Qj.surf_input, R.R_VARIABLES{f}, pertu_pos  );
    elseif aux_type(f) == 2
      Qj.atmos_input = setfield( Qj.atmos_input, R.R_VARIABLES{f}, pertu_pos  );
    end

    [y_pos] = fom_atmos_ocean_seaice_land_tbs_core( Qj.atmos_input, Qj.surf_input, Qj.sensor_input, Qj.dir_input, Q.ISO_DO  );

    %= negative perturbation

    pertu_neg      = x(f) * ( 1 - R.J_PERTURBATION(f)/100 );

    if pertu_neg < min_value(f)
      pertu_neg = min_value(f);
    end

    if aux_type(f) == 1
      Qj.surf_input  = setfield( Qj.surf_input, R.R_VARIABLES{f}, pertu_neg  );
    elseif aux_type(f) == 2
      Qj.atmos_input = setfield( Qj.atmos_input, R.R_VARIABLES{f}, pertu_neg  );
    end

    [y_neg] = fom_atmos_ocean_seaice_land_tbs_core( Qj.atmos_input, Qj.surf_input, Qj.sensor_input, Qj.dir_input, Q.ISO_DO  );

    %= dy/dx for J matrix

    if sum(y_pos-y_neg) ~= 0
      J(:,f) = (y_pos - y_neg) / (pertu_pos - pertu_neg); 
    else
      J(:,f) = 0;
    end

  end

end


return
