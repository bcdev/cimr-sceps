%============================================================================================
%
% FORMAT 
%	 [ e, cagb ] = fom_snowland_emis_surfem(freq, eia, lat, lon, mm, dd, lt, ...
%		   dircli, cl, sdor, asn, rsn, dsn, skt, tsn, stl1, smlt, cemis, cagb)
%
%
% DESCRIPTION
%
%	SURface Fast Emissivity Model for xontinental Snow Land (SURFEM-SnowLand) developed
%	developed by Iris de Gelis, Catherine Prigent, and Carlos Jimenez
%	as part of the ESA CIMR SCEPS project.
%
%	It provides the V-pol and H-pol emissivity (emis) for a given lat-lon location
%       and a number of frequencies and Earth incidence angles. There is no frequency
%	interpolation and the emis corresponds to the microwave imager channels 
%	1.41 6.925 10.65 18.7 or 36.5 GHz. There is an interpolation in frequency
%       and, althoug the emis parameterization is derived mainly using observations
%       at the SMAP (40deg) and AMSR2 (55deg) incidence angles, any incidence angle
%	between 0 and 90 can be given.
%
%  It  works in 6 modes: 
%
%	    (1) To output a climatological emis value from a file, derived from the
%	    SMAP and AMSR2 instruments at their respective local times and incidence
%	    angles of ~6AM/PM and 40deg, and 1.30AM/PM and 55 deg, respectivley. The 
%           local time is used to select the closest SMAP and AMSR2 overpass, i.e., 
%           to select the descending or ascending orbit of the intrument. For instance,
%           if lt = 11, the climatological emis is selected from the 6AM SMAP and
%           1.30PM AMSR2 overpass, depending on the choise of frequency. The variables 
%           lat, lon, mm, dd, and lt can be vectors with one element per location. This
%	    mode is activated by callng the function as
%
%              e = fom_snowland_emis_surfem( freq, [], lat, lon, mm, dd, lt, dircli  );
%	
%	    (2) To output a "climatological" emis value as from mode (1), but the 
%	    climatological emissivity is an average of the descending/ascending SMAP
%	    and/or AMSR2 overpasses. Local time is then not needed. This mode is
%           activated by calling the function as
%
%              e = fom_snowland_emis_surfem( freq, [], lat, lon, mm, dd, [], dircli  );
%
%           (3) To output a climatological emis value, only for one lat-lon location, but for 
%	    any given incidence angle, i.e., there is an internal interpolation from the
%           SMAP and AMSR2 observation angle to the selected incidence angle. This mode is
%	    activated by callng the function as
%	    
%              e = fom_snowland_emis_surfem( freq, eia, lat, lon, mm, dd, lt, dicli  );
%
%            (4) To output an actual emis value driven by the geo-physical inputs, 
%	     for one lat-lon location. Climatological emis, AGB values, are also 
%	     required for the actual emis derivation and are calculated internally.
%	     This mode is activated by calling the function as  	     
%
%	       e  = fom_snowland_emis_surfem(freq, eia, lat, lon, mm, dd, lt, dircli, cl,
%	                  sdor, asn, rsn, dsn, skt, tsn, stl1, smlt);	
%
%            (5) To output an actual emis value driven by the geo-physical inputs, 
%	     for one lat-lon location as for mode (4), but with climatological emis and AGB
%	     values passed as an input instead of being derived internally ans in mode (4).
%	     This mode is activated by calling the function as 
%
%	       e  = fom_snowland_emis_surfem(freq, eia, lat, lon, mm, dd, lt, dircli, cl,
%	                  sdor, asn, rsn, dsn, skt, tsn, stl1, smlt, cemis, agb);
%
%            (6) To output an emis value derived from the climatological value passed in
%	     cemis updated with the angular dependence. This mode is activated by calling 
%	     the function with an empty value in any of the cl, sdor, asn, rsn, dsn, skt, 
%            dsn, stl1, and smlt parameters and can be used to revert to an angle-interpolated
%	     external climatoloical emis if some of the inputs are detected as unrealistic 
%	     for some cases of the simulation, e.g.   
%
%	       e  = fom_snowland_emis_surfem(freq, eia, lat, lon, mm, dd, lt, dircli, cl,
%	                  sdor, asn, rsn, dsn, skt, [], stl1, smlt, cemis, agb);
%
%	     Notice that this mode is different from (3) as in (3) cemis and agb are not passed
%	     but internally calculated. 	
%
%
% OUT  
%	e                 emissivty values as a matrix of dimensions: 
%
%			    mode (1) and (2): 
%				 2 (v-pol and h-pol) x number of lat-lons x number of eia
%			    mode (3),  (4), (5), and (6): 
%				 2 (v-pol and h-pol) x number of freqs x number of eia
%
%
%
% IN    
%	freq	[GHz]     Frequency, it should be one of the CIMR frequencies
%			  (1.41 6.925 10.65 18.7 or 36.5), it can be passed
%			  as just one frequency or a vector with a number of them. 
%			  Notice that no interpolation is done in frequency, so the 
%			  frequencies have to match the given list.
%
%	eia	[deg]     Earth incidence angle between 0 and 90, it can be passed as 
%			  a evtor of any length.
%
%	lat     [deg]	  Latitude, passed as a vector of any lenght for modes (1) 
%			  and (2), or a single value for remaining modes.
% 
%	lon     [deg]	  Longitude between -180 and +180, passed as a vector matching  
%			  the lon dimensions for modes (1) and (2), or a single value 
%			  for remaining modes.
%
%	mm                Numeric month of the year, single value.
%
%	dd                Numeric day of the month, single value.
%
%	lt      [hour]    Local time in decimal hours, passed as a vector matching
%			  the lat dimensions for modes (1), as [] for mode (2), or 
%			  a single value for remaining modes.
%
%	dircli		  Folder where the ancillary files where the climaological
%			  emis are found. 
%			  
%			  Only relevant for modes (4), and (5), single values
%
%	cl	[-]	  Lake Cover
%			      (ERA5 "cl"  used in paremeterization derivation)
%	sdor	[m]       Orography
%			      (ERA5 "sdor"  used in paremeterization derivation)
%	asn	[-]	  Snow Albedo
%			      (ERA5 "asn"  used in paremeterization derivation)
%	rsn	[Kg/m3]	  Snow Density
%			      (ERA5 "rsn"  used in paremeterization derivation)
%	dsn	[m]       Snow depth 
%			      (Product of rho_water * sd ./ rsn where rho_water
%			      is the water desnity, sd is ERA5 snow depth in meters of water 
%			      equivalent and rsn is ERA5 snow density, this dns used 
%			      in paremeterization derivation)
%	skt	[K]       Skin Temperature
%			      (ERA5 "skt"  used in paremeterization derivation)
%	tsn	[K]	  Snow Temperature
%			      (ERA5 "tsn"  used in paremeterization derivation)
%	stl1	[K]	  Soil (Layer 1) Temperature
%			      (ERA5 "stl1"  used in paremeterization derivation)
%	smlt	[m]	  Snow Melt  
%			      (ERA5 "smlt"  used in paremeterization derivation)
%
%			  Only relevant for modes (5) and (6)	
%
%	cemis	[-]	  Climatological emissivity
%			      (it should have been extracted from the files
%			      in dircli and be passed as a matrix 
%			      2 (v-pol and h-pol) x number of freqs)
%	cagb	[Mg/ha]	  Climatological AGB
%			      (it should have been extracted from the files
%			      in dircli, single value)
%
% REFERENCES
%
%        Kilic, L.,, Prigent, C., Jimenez, C., and Sandells, 
%	 Forward modelling of passive microwave emissivities over snow-covered
%	 areas at continental scale, https://doi.org/10.1016/j.rse.2025.114821DOI
%
%	 Funded by ESA contract 4000138129/22/NLIA. 
% 
% HISTORY
%
%	2024-12-10	iris.degelis(at)estellus.fr: coded as part of the ESA CIMR SCEPS project, 
%
%	2025-09-01	carlos.jimenez(at)estellus.fr: adding some code documentation and minor fixes
%
%	2025-10-15	carlos.jimenez(at)estellus.fr: removing sunza as input
%
%============================================================================================


function [ e, cagb ] = fom_snowland_emis_surfem(freq, eia, lat, lon, mm, dd, lt, dircli, cl, sdor, asn, rsn, dsn, skt, tsn, stl1, smlt, cemis, cagb)

nin = length(eia);
len = length(freq);
nca = length(lat);

%=== checkin local time
%    local time can be given wrt to a given day, i.e., it can be negative or > 24h
%    e.g. -2 it means 22h of day before, here it is needed in 0-24


ilo = find( lt < 0);
iho = find( lt >= 24 );
lt(ilo) = lt( ilo ) +24;
lt(iho) = lt( iho ) -24;


%=== Get climatological emissivities and agb for multiple location and
%    output them for further processing. They can be then input to calculate
%    the dynamical emis and avoid multiple loading of the files per location.
%    This will be faster for many emis calculations as the climato emis
%    files are only read one. These emis are at the SMAP and AMSR2 incidence
%    angles, no frequency dependence added, as this is the input to the basic
%    emis calculation. To output the climato emis with the angular dependence
%    use the standard mode with nargin == 8. 

if isempty( eia )

  no = length( lon );
  e   = nan( 2, no, len);
  cagb = nan( 1, no );

  for k = 1:len

    ifreq = freq(k);
    ico   = nan(1,no);
    iro   = nan(1,no);

    %=== loading climatological emis and agb 
    %    using 0LT for desc and 17LT for asc
    %    with the first lat and lon to decide 
    %    EASE hemisphere 

    [doy_clim, grid, sat, group, orbit, ~, ~, idx_v, idx_h ] = which_climato(mm, dd, 0, lat(1), lon(1), ifreq);
    path_clim_file = fullfile(dircli, 'Climatologies', grid, sat, sprintf("Emissivity_for_%s_%s_%s_%s_%03d.mat", sat, group, orbit, grid, doy_clim));
    climdes = load(path_clim_file).data_desc;
    climdes = squeeze( climdes(1,:,:, [ idx_v idx_h ] ) );

    [doy_clim, grid, sat, group, orbit, ~, ~, idx_v, idx_h ] = which_climato(mm, dd, 17, lat(1), lon(1), ifreq);
    path_clim_file = fullfile(dircli, 'Climatologies', grid, sat, sprintf("Emissivity_for_%s_%s_%s_%s_%03d.mat", sat, group, orbit, grid, doy_clim));
    climasc = load(path_clim_file).data_asc;
    climasc = squeeze( climasc(1,:,:, [ idx_v idx_h ] ) );

    if k == 1 
      path_agb_file  = fullfile(dircli, sprintf("CCI_AGB_%s.mat", grid));
      iagb = load(path_agb_file).gridded_tab;
    end  
  
    %=== finding index locations and asc or des as function of LT

    [iro, ico] = latlon_to_ease2(strrep(grid, 'p', '.'), lat, lon);


    %=== climato using the local time to decide if taking from asc or desc files

    if ~isempty( lt )

      if strcmp( group, 'SMAP' )

        % if lt<12
        % 6.00 AM
        % orbit = 'DESC';
        % else
        % 6.00 PM
        % orbit = 'ASC';
        % end
        ult = 12; 
        for o = 1:no
          if lt(o) < ult
            e(:,o,k) = climdes( iro(o),ico(o),:);
          else 
            e(:,o,k) = climasc( iro(o),ico(o),:);
          end
          cagb(o)   = iagb(iro(o),ico(o)); 

        end
 
      elseif strcmp( group, 'GCOM-W1' )

        % if lt<7.5 || lt>=19.5
        % 1.30 AM
        % orbit = 'DESC';
        % else
        % 1.30 PM
        % orbit = 'ASC';
        llt = 7.5;
        ult = 19.5; 
        for o = 1:no
            if lt(o) < llt | lt(o) >= ult
              e(:,o,k) = climdes( iro(o),ico(o),:);
            else 
              e(:,o,k) = climasc( iro(o),ico(o),:);
            end
            cagb(o)   = iagb(iro(o),ico(o)); 
        end

      end

    elseif isempty( lt )

      %=== climato as an average or asc and desc 

      clim = (climasc + climdes)/2; 

      for o = 1:no
        e(:,o,k) = clim( iro(o),ico(o),:);
        cagb(o)   = iagb(iro(o),ico(o)); 
      end

    end 

    %=== filling cimato gaps before doing this
    %    at original files 

    for r = 1:3

      for p = 1:2

        aux = squeeze( e(p,:,k) );
        ind = find( isnan(aux) == 1 );

        if ~isempty(ind) 
          for i = ind
            if ind > 1 & ind < nca
              aux(i) = nanmean( [ aux(i-1) aux(i+1) ]);
            end
          end
        end 

        e(p,:,k) = aux;

      end  

    end 

  end  %k

  return

end



%=== Get climatological emissivities for a given location to be used as
%    input of the internal derivation of dyanmical emis or used a 
%    pre-calculated value from the inputs
 
if ~exist( 'cemis', 'var' )

  clim_v = nan( 1,len);
  clim_h = nan( 1,len);

  for k = 1:len

    %=== identifying climatological emis to load

    [doy_clim, grid, sat, group, orbit, idx_col, idx_row, idx_v, idx_h] = which_climato(mm, dd, lt, lat, lon, freq(k));
 
    %=== emis file

    path_clim_file = fullfile(dircli, 'Climatologies', grid, sat, sprintf("Emissivity_for_%s_%s_%s_%s_%03d.mat", sat, group, orbit, grid, doy_clim));
    
    %=== loading emis from file

    if ~isempty( strfind( path_clim_file, 'ASC') )
      clim = load(path_clim_file).data_asc;
    else
      clim = load(path_clim_file).data_desc;
    end
    clim_v(k) = clim(1, idx_row,idx_col, idx_v);
    clim_h(k) = clim(1, idx_row,idx_col, idx_h);

  end


else

  clim_v = cemis(1,:);
  clim_h = cemis(2,:);

end




%=== dynamical or cimatological emis

if nargin > 8 & ~isempty( cl + sdor + asn + rsn + dsn + skt + tsn + stl1 + smlt)

  %=== calculate dynamical emissivities if inputs available

  %= get agb if not input

  if ~exist( 'cagb', 'var' )

    [doy_clim, grid, sat, group, orbit, idx_col, idx_row, idx_v, idx_h] = which_climato(mm, dd, lt, lat, lon, freq(end));
    path_agb_file  = fullfile(dircli, sprintf("CCI_AGB_%s.mat", grid));
    cagb = load(path_agb_file).gridded_tab;
    cagb = cagb(idx_row,idx_col);
   
  end

  %= emis per frequency and angle

  e_v = nan( nin,len);
  e_h = nan( nin,len);

  for k = 1:len

    inputs = [ cl, sdor, cagb, asn, rsn, dsn, skt-stl1, tsn, stl1, smlt, clim_v(k), clim_h(k)]';
    [e_v(:,k),e_h(:,k)] = get_emis_freq_incident_angle(inputs, freq(k), eia);

  end

else

  %=== adding angular dependence to climato emis

  e_v = nan( nin,len);
  e_h = nan( nin,len);

  for k = 1:len

    %= adding angular dependence 

    [e_v(:,k),e_h(:,k)] = get_emis_freq_incident_angle([], freq(k), eia, [clim_v(k);clim_h(k)]);

  end

end


%= emis returned as v-h x n freqs x n ang

e = nan( 2, len, nin );
e(1,:,:) = e_v';
e(2,:,:) = e_h';

return


%------------------------------------------------------------------------


function [e_v,e_h] = get_emis_freq_incident_angle(inputs, freq, theta, emis)
switch nargin
    case 3
        %%---Emissivity computation at 55° for frequency from 6.9 to 36.5 GHz, and 40° for 1.4 GHz --------
        [ W1, b1, W2, b2, ti, to ] = get_nn_weight(freq);
        [emis] = apply_nn(inputs, W1, b1, W2, b2, ti, to);
    case 4
        % Emissivity taken from input emis, inputs parameter is not used
        % It can be set to NaN
    otherwise
        error("Check function inputs.")
end
if freq==1.41
    % At 1.4 GHz, polynom (degree 2) derived from SMOS angular dependence and fitted
    % to SMAP. Adjustements on the value at nadir (pv(3)=ph(3)=e0) and then adjustements on pv(2) and ph(2).
    % pv(1) and ph(1) not changed (the curvature).
    [pv_ori, ph_ori] = polynomial_1ghz;
    pv(3,:)=emis(1,:)-pv_ori(2)*40-pv_ori(1)*40^2;
    ph(3,:)=emis(2,:)-ph_ori(2)*40-ph_ori(1)*40^2;
    pv(3,:)=(pv(3,:)+ph(3,:))/2.; % V+H/2 --> on suppose qu'à 40° on est à la moitié
    ph(3,:)=pv(3,:);
    pv(2,:)=(emis(1,:)-pv(3,:)-pv_ori(1)*40^2)/40;
    ph(2,:)=(emis(2,:)-ph(3,:)-ph_ori(1)*40^2)/40;
    pv(1,:)=pv_ori(1);
    ph(1,:)=ph_ori(1);
    for i=1:size(emis,2)
        e_v(:,i)=polyval(pv(:,i),theta);
        e_h(:,i)=polyval(ph(:,i),theta);
    end
elseif freq>1 && freq<18
    % At 6 and 10 GHz, the angular dependence is derived from the 1.4 GHz dependence (above) and from TELSEM2 at 18 GHz.
    % For TELSEM2, class8 is selected after some tests: limited differences between class 7 to 9 and the intermediate one is
    % selected.
    % - e0 is a (frequency proportional) mix between the value calculated from the 1.4GHz angular dependence (corresponding at the 55°) and e0 directly calculated from TELSEM
    % - curvature (theta*theta) is from 1.4 GHz derived from SMOS  (all SIT averaged).
    % - the other coeff (theta) to fit the estimated values at 55°
    [pv_ori, ph_ori] = polynomial_1ghz;
    [a0_k0, a0_k1, a0_k2, a0_eveh, a1_eveh, a2_eveh, a3_eveh, b0_eveh, b1_eveh, b2_eveh, b3_eveh] = telsem2_coef;
    theta55 = 55;
    e0_1v=emis(1,:)-pv_ori(2)*theta55-pv_ori(1)*theta55*theta55;
    e0_1h=emis(2,:)-ph_ori(2)*theta55-ph_ori(1)*theta55*theta55;
    e0_1=(e0_1v+e0_1h)/2.;
    e0_18 = a0_k0(1,8)+a0_k1(1,8)*emis(1,:)+a0_k2(1,8)*emis(2,:);
    pv(3,:)=(freq-1.4)*(e0_18-e0_1)/(18-1.4)+e0_1;
    ph(3,:)=pv(3,:);
    % Coef (in theta) to fit the value at emis value at 55°.
    pv(2,:)=(emis(1,:)-pv(3,:)-pv_ori(1)*theta55*theta55)/theta55;
    ph(2,:)=(emis(2,:)-ph(3,:)-ph_ori(1)*theta55*theta55)/theta55;
    % Curvature (theta*theta) is from 1.4 GHz derived from SMOS.
    pv(1,:)=pv_ori(1);
    ph(1,:)=ph_ori(1);
    e_v=polyval(pv(:,1),theta);
    e_h=polyval(ph(:,1),theta);
else
    % At 18, 36 and 89 GHz, TELSEM2 with a modulation on the choice of the
    % class depending on emissivity obtained at 55°
    [a0_k0, a0_k1, a0_k2, a0_eveh, a1_eveh, a2_eveh, a3_eveh, b0_eveh, b1_eveh, b2_eveh, b3_eveh] = telsem2_coef;
    class1=6:9;  % class1=6 to 9 from TELSEM2
    %Class center according to frequencies (set empirically)
    mean_emis_threshClassfreq = [[0.92 0.85 0.75 0.65]
        [0.82 0.81 0.8 0.75]
        [0.8 0.75 0.7 0.5]];
    telsem_idx = dictionary([18.7, 36.5, 89],[1,2,3]);
    j=telsem_idx(freq); % 18 GHz is the first channel in TELSEM2

    mean_emis = mean(emis, 1);

    % Compute the distance to the class
    coef = zeros(4,size(emis,2));
    % class 6 very low emis
    coef(1,:) = abs(mean_emis_threshClassfreq(j,4)-mean_emis);
    % class 7
    coef(2,:) = abs(mean_emis_threshClassfreq(j,3)-mean_emis);
    % class 8
    coef(3,:) = abs(mean_emis_threshClassfreq(j,2)-mean_emis);
    % class 9
    coef(4,:) = abs(mean_emis_threshClassfreq(j,1)-mean_emis);

    % if mean_emis is larger or below the class 6 or 9 only chose the
    % corresponding one
    % Class 9
    coef(1:3, mean_emis > mean_emis_threshClassfreq(j, 1)) = 1;
    % Class6
    coef(2:4, mean_emis < mean_emis_threshClassfreq(j, 4)) = 1;

    % Choose the two nearest classes
    [~, idx] = sort(coef, 1);
    coef(sub2ind(size(coef), idx(3, :), 1:size(coef, 2))) = 1;
    coef(sub2ind(size(coef), idx(4, :), 1:size(coef, 2))) = 1;

    % Normalize coefficients
    coef = 1 - coef;
    coef = coef ./ sum(coef, 1);

    for c=1:numel(class1)
        cl = class1(c);
        e0 = a0_k0(j,cl)+a0_k1(j,cl)*emis(1,:)+a0_k2(j,cl)*emis(2,:);
        a0 = a0_eveh(j,cl);
        a1 = a1_eveh(j,cl);
        a2 = a2_eveh(j,cl);
        a3 = a3_eveh(j,cl);
        b0 = b0_eveh(j,cl);
        b1 = b1_eveh(j,cl);
        b2 = b2_eveh(j,cl);
        b3 = b3_eveh(j,cl);
        theta0 = 0.;
        theta55 = 55.; % coefficient calculated for inputs at 53° in TELSEM2, but little difference expected with 55° (given the error bar).

        %%% Vertical polarization
        S1_v = ((theta-theta55)./(theta0-theta55)) * ((e0-a0)/a0);
        em55_v = a3*(theta55^3) + a2*(theta55^2) + a1*theta55 + a0;
        S2_v =((theta-theta0)./(theta55-theta0))*((emis(1,:)-em55_v)/em55_v);
        S_v = 1 + S1_v + S2_v;
        emtheta_v = a3*(theta.^3) + a2*(theta.^2) + a1*theta + a0;
        emisV_c(c,:) = S_v .* emtheta_v;
        %%% Horizontal polarization
        S1_h = ((theta-theta55)./(theta0-theta55)) * ((e0-b0)/b0);
        em55_h = b3*(theta55^3) + b2*(theta55^2) + b1*theta55 + b0;
        S2_h =((theta-theta0)./(theta55-theta0))*((emis(2,:)-em55_h)/em55_h);
        S_h = 1 + S1_h + S2_h;
        emtheta_h = b3*(theta.^3) + b2*(theta.^2) + b1*theta + b0;
        emisH_c(c, :) = S_h .* emtheta_h;
    end
    e_v = emisV_c'*coef./sum(coef);
    e_h = emisH_c'*coef./sum(coef);
end

return

%------------------------------------------------------------------------

function [doy_clim, grid, sat, group, orbit, idx_col, idx_row, idx_v, idx_h] = which_climato(mm, dd, lt, lat, lon, freq)
% IN     mm [-]             Month
%        dd [-]             Day
%        lt [-]             Decimal Local time
%        lat [°]            Latitude between -90 and 90°
%        lon [°]            Longitude between -180 and 180°
%        freq [GHz]         Frequency, it should be one of the CIMR frequency : 1.41 6.925 10.65 18.7 or 36.5
doy_clim = tge_doy_10day_composite(1995, tge_dayofyear(1995, mm, dd));
if lat >= 0
    % EASE2_N12p5KM
    grid = 'EASE2_N12p5KM';
else
    % EASE2_S12p5KM
    grid = 'EASE2_S12p5KM';
end

[idx_row, idx_col] = latlon_to_ease2(strrep(grid, 'p', '.'), lat, lon);

if freq == 1.41
    %SMAP
    sat = 'SMAP';
    group = 'SMAP';
    idx_v = 1;
    idx_h = 2;
    if lt<12
        % 6.00 AM
        orbit = 'DESC';
    else
        % 6.00 PM
        orbit = 'ASC';
    end
else
    %AMSR2
    sat = 'AMSR2';
    group = 'GCOM-W1';
    if lt<7.5 || lt>=19.5
        % 1.30 AM
        orbit = 'DESC';
    else
        % 1.30 PM
        orbit = 'ASC';
    end
    amsr2_freq_to_idx_v = dictionary([6.925, 10.65, 18.7, 36.5],[1, 3, 5, 8]);
    idx_v = amsr2_freq_to_idx_v(freq);
    idx_h = amsr2_freq_to_idx_v(freq) + 1;
end

return


%------------------------------------------------------------------------
%
% NAME:    [ Y ] = apply_nn( X, W1, b1, W2, b2, ti, to)
%
%          Calculates a Jacobian matrix from analytically deriving the
%          deriative of the outputs wrt the inputs. If [ni] is the number
%	   of net inputs, [no] the number of net outputs, [nh] the number
%	   of nodes in the hidden layer, and [ns] the number of samples
%	   where to calculate the net output and Jacobian (J), J is a matrix
%	   [ns * no * ni], where the element J[x,o,i] gives the derivative
%	   of output o wrt input i evaluated for sample s.
%
%          net has to be a MLP with mapminmax applied to input and
%	   outputs, a first (hidden) layer with tansig activation functions
%	   and a second (output) layer with purelin functions. See the code
%	   below to see how the NN, mapminmax, and J, are implemented as
%	   matrix multiplications.
%      Note :  the function does not return jacobians as it is not of
%      interest for this forward model.
%
% OUT:	Y	net output for X, matrix [no * ns]
%
% IN:	X   net input, matrix [ni * ns]
%	W1	weight matrix of first layer, matrix [nh * ni]
%	b1	bias vector of first layer, vector [nh * 1]
%	W2	weight matrix of second layer, matrix [no * nh]
%	b2	bias vector of first layer, vector [no * 1]
%   ti	structure having the settings of the input
%		mapminmax transformation
%   to	structure having the settings of the output
%		mapminmax transformation
%
%------------------------------------------------------------------------

% HISTORY: 20220308  Created by carlos.jimenez(at)estellus.fr

function [Y] = apply_nn( X, W1, b1, W2, b2, ti, to )


%=== NN output with matrices

% number of samples
nx = size( X, 2 );
% number of inputs nodes
ni = size( W1, 2 );
% number of outputs nodes
no = size( W2, 1 );
% number of nodes in hidden layer
nh = size( W1, 1);


%= input transformed with mapminmax
%      y = (ymax-ymin)*(x-xmin)/(xmax-xmin) + ymin;
%  but coded with a ymin, offset and gain as
%      y = iymi + (x - offset).* gain;

igai = ti.gain;
igai = repmat( igai, 1, nx);
ioff = ti.xoffset;
ioff = repmat( ioff, 1, nx);
iymi = ti.ymin;
iymi = repmat( iymi, 1, nx);

aux = iymi + (X - ioff).* igai;

%= propagation through first with layer W1 b1 and tansig

x1  = W1 * aux + b1;
aux = 2 ./ (1+exp(-2 * x1 ) ) -1;


%= propagation through output layer with W2 and b2 and purelin

aux = W2 * aux + b2;

%= output transformed back with mapminmax

ogai = to.gain;
ogai = repmat( ogai, 1, nx);

ooff = to.xoffset;
ooff = repmat( ooff, 1, nx);

oymi = to.ymin;
oymi = repmat( oymi, 1, nx);

Y = ((aux - oymi ) ./ ogai ) + ooff;

return

%--------------------------------------------------------------

function [ W1, b1, W2, b2, ti, to ] = get_nn_weight( freq )

if freq==1.41
    [ W1, b1, W2, b2, ti, to ] = nn_1GHz;
elseif freq==6.925
    [ W1, b1, W2, b2, ti, to ] = nn_6GHz;
elseif freq==10.65
    [ W1, b1, W2, b2, ti, to ] = nn_10GHz;
elseif freq==18.7
    [ W1, b1, W2, b2, ti, to ] = nn_18GHz;
elseif freq==36.5
    [ W1, b1, W2, b2, ti, to ] = nn_36GHz;
else
    error("NN is not available of the following frequency : %f.", freq)
end

return


%------------------------------------------------------------------------
%         Hard Coding of NNs weights (One NN per frequency)
%------------------------------------------------------------------------
function [ W1, b1, W2, b2, ti, to ]   = nn_1GHz


ti.gain = [2.23681390713175;0.00224835344544316;0.00592135020084487;5.56242940128978;0.00578469335326059;0.0603493854712451;0.0265540636497697;0.0358927707854227;0.0379556556596848;381.669038172038;3.47732194317967;2.93556160418144];

ti.xoffset = [7.62939453125000e-06;0.286437988281250;0;0.519999742507935;101.144927978516;0.0200005844235420;-43.5700988769531;217.487304687500;227.105056762695;-2.32830643653870e-10;0.530708312988281;0.403301596641541];

ti.ymin  = -1.0000000e+00;

to.gain = [3.231205917439197;2.810952083618423];

to.xoffset = [0.499767780303955;0.383983105421066];

to.ymin  = -1.0000000e+00;


W1 = [0.965524693820893,-1.11586396673938,0.188415109032025,-0.180649414044249,-0.296211206371305,-0.236382514768432,-0.715749858252247,-0.0140373630162679,-0.0994251760080192,0.259966566727509,-1.37053409570767,-0.319921436797090;0.0653660948530833,-0.0977245242596171,0.168524149023028,-0.217228681980222,0.105381938574446,0.140771009431585,0.0925519764724498,-0.0554622088528149,0.103679456163968,0.0193097614427350,0.116546439510736,0.114187778268425;-0.348408221110119,0.848424242336654,-0.0204030826729300,-0.0117987373463156,-0.146461822715321,-0.244660850114419,-0.219836857295417,-0.313847100557629,-0.873617292218878,0.216652510082708,1.37470688046133,1.44900621449207;0.117549678953447,0.177516466899407,-0.255061803190362,0.0844781435510945,-0.247177649586194,-0.463906959768987,-0.206979723432727,-1.04145857197902,0.195522332599969,0.355072745727289,-0.328333604393528,0.319454887859760;0.0953718422564140,-0.0208652369898233,-0.0554957059133861,0.371404470694082,-0.796222091130545,1.02055540339257,0.0670657490378346,0.224428098704957,-0.528608990313826,-0.144679453012948,0.520747756005406,0.148645370560308;0.0221139353943182,0.110646171628136,-0.245564749317848,-0.0845119748354057,-0.0477691242831042,0.313416565286261,0.230051695814575,-0.200352050052452,-0.158986515523558,0.0218372793660211,0.295112479895986,0.241921763831430;-0.0525611768309498,-0.214418612844544,0.230572464852586,-0.0308400432326152,0.295672012375322,-0.445127093087154,0.119608721729002,-0.193224871791297,-0.266643731935846,-0.467934261295901,0.00415624939665890,0.138399862201949;-0.00555519597393269,-0.494165046812970,0.887817141151077,0.0867288366309973,-0.419854877947200,-1.29583717737604,0.0407100867929391,0.820906928787809,0.133056241009044,0.0940645490048624,0.977785096605630,-1.08927375488647;-0.278321583246763,1.67891885755353,-0.00296088832098247,-0.118828030784686,-0.108506353228679,-0.100073414718063,-0.343564858506602,0.480544573185194,-0.149957508521244,-0.0284244733563945,-1.02415509873978,1.59057665992245;0.0366265610707355,-0.0469476541284881,0.0655933565056296,-0.141671947755266,0.248932848394767,-0.482921099433639,-0.123395252010168,-0.0205115243417998,0.0147096308170709,0.873030402393290,-0.112866319183703,-0.151759248121547;0.284647068077448,-0.261110654344001,0.258833923511789,0.541659346379042,0.535336145600254,0.138944493979078,-0.679814973867712,-0.0463521502341879,0.119601162144182,0.134587120956064,0.523749725888736,-0.719247183776870;0.0132066734926130,0.244856229684700,0.0708489954529437,-0.548881987573209,-0.846734854496444,0.0149587891025951,-0.647020102331977,0.380126410798934,-1.12967273801462,-0.451787154709674,1.07123432716015,-0.134250536910403;-0.460288217422207,-0.632146586891172,-0.332593719270739,-0.0513568691348368,-0.0298635562563478,0.254121144133306,0.572395449270294,0.640409785382539,0.538349541318094,0.0809832254840536,-0.133473810045895,-0.285316625956622;-0.109644926977655,0.0996167455268921,-0.174501519487552,0.274320502571272,-0.107865417314716,-0.194559299650497,-0.136064839676712,0.0360840379058540,-0.225719514209606,-0.0707353732502252,-0.161142068129241,-0.133470094832261;0.319897871512618,0.226021520748132,0.167553278640455,0.115656155881358,-0.213066173554875,0.0735347736375097,-0.119381682258988,-0.294929875541522,-0.0132507257809987,0.0353984286211386,-0.173683904569791,0.00698029227305931;-0.171910458772671,-0.825223402771280,0.788653829175405,0.0746370653511705,-0.686354299017070,-1.84436328065098,0.274222269417393,0.553683912949067,0.633316019740090,0.0773234823769263,0.891364662812824,-1.22333056533763;-0.0261391386433171,0.429708540101203,0.101341151525105,-0.449968037157189,-0.219042250127090,0.165633655148765,0.646755153506047,0.266928004473473,0.277325120577842,-0.334417274794195,0.669079300685630,-0.0558964179766791;-0.604402909010188,-0.712358911346828,-0.598367384504074,-0.120887219215368,0.310476777324290,0.201065384375194,-0.731698865032615,0.196841456456709,-1.03743254378452,-0.213355831293279,-1.79628923820754,-0.731125120490849;-0.0151592549508896,0.0673385536323225,-0.226957193880657,-0.00815181318904842,-0.201648980006006,0.575201060147137,0.291677976470368,0.446404150331963,0.483148678357082,0.326903350924069,0.318363286077220,-0.623949681661674;0.558791137829126,-2.05224518718954,0.124693144265022,0.719112077617956,-1.51956811332187,1.03014282864532,-0.741666309118768,0.673670683965005,-0.188044651596942,0.0305137794890878,-0.223726560809573,0.215078266851072;0.158300965166184,-0.716147947795252,0.296287389721957,0.0788108322984573,-0.188902072910942,0.336796604995627,0.104692365862030,0.0535232988830103,0.579755171367540,0.261495059174148,-0.0806255899842336,0.0675168919403982;-0.373298799018554,-0.307020945077744,0.499759810877272,0.220978517282524,-0.830939400671433,-1.13493217500707,-0.200772705410072,0.411918049300321,0.425488434666366,0.506068708844131,0.0994044235739751,-0.407996917591264;0.406386737922421,1.03504324732525,0.199174418266423,0.538864061675678,-0.429512407737835,-0.0913337457439226,0.311689246918743,0.513663668883252,0.421288165027111,0.425744215774758,0.692985715929407,0.156473328395931;0.643399053557273,-0.0262702179195163,0.302975817556746,-0.214691877034396,-0.201127878337930,-0.165327539293881,-0.403198903691087,-0.814336844222251,-1.09619625581485,-0.757095922941083,-0.181774523608822,1.14829612041959;0.0853106899194011,0.120082055581727,-0.0738102233214614,0.0183591779676807,-0.109493214263496,-0.127150451504496,-0.102115756889986,-0.233027252027247,-0.238696440441952,-0.0416763764171903,0.0654491408474886,0.185998929005005;-0.114357779170288,-0.573946263416811,-0.769967538449896,0.0949404780123119,0.180236188706711,0.0849125790804945,0.475249475084722,0.876932697159728,1.23076898375690,0.256320315941834,-0.760528311718766,-0.929850013022582;-0.512783119289366,-0.390586324049636,-0.361975075639408,0.253369091712492,0.321390972740674,0.109150738166650,0.616027936117265,0.991814534305267,1.76832413995650,0.369686659136325,0.719185376335385,-1.92000979087107;-0.706656061477177,-0.425146972869303,0.113411848037939,-0.255210613314157,-0.315315923017320,-0.274784125686780,0.120673966169628,0.114779910114169,0.166435681150047,0.213338890160755,-0.0302068308204625,0.216422799393769;0.928295464722280,-0.254190279089193,0.129218393290328,0.0498267954971235,-0.304599419784509,0.0587387135048356,-0.148753049810947,-0.105075859600420,-0.237812160292084,0.379088734370509,-0.741892151899780,-0.554724315428649;0.922726930718701,-0.824449483729856,0.291551121360720,-0.0393598598817591,-0.116057502582331,-0.0367173486046726,0.329672930610421,-0.172110584554599,-0.0504332944648284,-0.501214183056099,0.425717506892993,1.40382940650892;0.106647697622225,-0.110562800004106,0.111959758788092,-0.178291145853303,0.143743984364237,0.174666374073170,0.148361836346396,0.0806649937309200,0.311790609120008,0.129377224468152,0.0792976302684040,-0.00569286655511679;0.691139877556893,-1.07703047900400,-0.0577781589387638,-0.312338237457780,-0.260864842402796,-1.59243958943158,-0.702479573356674,-0.238755336603287,0.0401181517544578,-0.0629907094229129,-2.96464268172741,0.690260675613136;-0.0381936694167963,-0.376361460563157,0.535952043768433,0.182368144727930,-0.0723502651397994,0.264873857863333,0.254120240961559,0.0278337687688772,1.09298176153345,0.0894677889888399,1.52726632882657,-0.484270850789656;0.0479337305606372,-0.0475069571399105,-0.149860572513754,-0.179179612334433,0.213648715092253,-3.85738315682939,-0.250961977824438,-0.178439130132508,-0.549197542461103,1.32092792680968,0.163138700100649,-1.01871526321960;-1.37995020971155,-0.836178990111510,0.897278207439133,0.470215888140832,-0.00790373307090842,1.08075576361852,0.0768790149067133,0.286021638010585,0.498045891716934,0.215116183634710,-0.181838098676592,-0.539672917865195;0.275367815504972,0.433193933144342,-0.287425746149740,-0.0607350502595600,0.0331383477342719,-0.0838043443661409,-0.101821981959482,0.127879543050834,-0.156766631739232,-0.125389040748628,-0.0351654792097912,-0.0808084443025734;-0.902404189056527,-0.758844655958247,0.978224538951487,0.508750899088103,0.456880509746135,-1.14965353102211,0.668597846267021,-0.0741745059129929,0.522963922565937,0.927321946778998,-1.30366384717730,1.01162596235449;0.355580476722031,-0.653075839731077,-0.291599393033340,-0.612304068655333,0.392532779684709,-0.297617003537327,0.668295605745645,-0.405595552017165,0.474802781024110,-0.399943549408057,-0.708395259287209,0.780676847176284;0.817795047548549,-0.0505736382111137,0.719213425470973,-0.331981971498858,0.169550552124850,0.491611287532216,0.224097749158230,-0.718230796072190,0.0400550198917802,-0.151257463863479,0.828606366176056,-0.354251312642848;0.603185973796043,0.457325428128126,1.13330083161604,0.0835611020838190,0.507765595480572,-0.360912693381173,0.233256543183018,0.206838528392257,0.166211941980196,-0.538362502329489,-0.180822749720100,0.711306206021600;0.670992160104029,-0.622453008518028,0.313275869421154,-0.295254029886421,0.918264737281106,-0.579258029543222,-0.0426001687921294,0.826729561448235,-0.00476451653955206,0.372012580313259,-0.764190844188430,-1.12394500950891;-0.792978570281972,-1.73373031731528,0.903242666463791,0.0770294088345450,-0.284749359756322,0.335813540576813,0.213221320477938,-0.0962005419884282,-0.231312950665306,-0.504340181376235,-1.01369719669025,0.349625084224371;-0.522998678354741,0.505671013063340,0.523763636846166,0.0814049246393252,-0.186655396466471,-0.0150441724964412,-0.773910822029058,0.234173070890518,-0.251958561359248,0.494645808636669,0.346469324810080,1.05110954129998;0.758203847812698,-0.550126939401106,0.0345285231935942,0.341514115431129,-1.80385571315806,0.176019707726764,-0.113629225810967,-0.00604402885002141,-0.758121614277282,0.458588911883957,0.777936853634403,0.648590068426131;0.260314082491134,0.249484136182408,0.295484435934605,0.247094750552537,-1.09570337678267,-0.442483533981931,0.523684173050810,-0.00129632049215802,0.960650110611988,-0.235625559718386,0.396979568010019,0.190940055403847;-0.348830827524558,0.182052752120767,-0.231057674902000,0.197524738834483,-0.403248809180526,-0.447305626397982,0.478299563340577,-0.716428770580765,0.275551556368221,-0.729225654858301,-0.0257235882425493,0.609377962376685;-0.301956898209733,-0.817473050047063,0.0501726714324962,-0.644035565874807,0.883110915336698,0.104444156941867,-0.245842121958858,-0.841719005133676,-0.0542648994983347,-0.603128887605675,0.599453857498666,-1.41012055644129;-0.00711319693478412,0.619245838051097,-0.474047036035035,-0.0666421197487532,-0.409111506984857,-0.0148224152198155,0.103033176547031,-0.542574793225302,0.174845107294918,0.281553384106359,-0.0273036687521167,0.906772528560643;-0.548715582237643,-0.0142902905976278,-1.08040983403178,0.143912005660553,0.195632272984803,-0.294656007920083,-0.387384166223857,0.742149317708586,-0.228321875226102,0.233536124983171,1.42923505802689,-2.44864437427282;0.108155011196878,-0.193033558976472,0.247231529775320,0.137919777553691,0.111382359709746,0.317320372206938,-0.139589297185216,0.0270142519015034,0.623212479046401,0.301044646920532,-0.0227884940347804,0.208299788889406;0.119112849923140,-0.359402660641419,0.464286324203809,0.286177322895412,-0.543688075560781,0.667361123496349,-0.0355007570856633,0.240772082361125,0.590015340213668,0.275649968384753,-0.410940492504559,0.100190033360760;0.0596531348423197,0.289008526030686,-0.00875485697470159,0.181006886546905,-0.0846612480684167,-0.00831286545305335,-0.100744351038012,0.0356618070066742,-0.117937391462655,-0.325658309983187,-0.160960892956667,0.0500393949809641;-0.342782566431766,0.290659094869693,0.173613125039792,0.168547504180083,-0.562072724506256,1.15080983810843,0.314056806317988,0.359855586025209,0.515460680625599,0.0173678154448598,0.101505055089530,-0.629962440518383;0.566349711234431,2.30364890913314,0.148795339118818,0.0786572106238032,-0.0612369900943309,-0.162119311341889,-0.856061563764432,0.658852326872775,-1.16546856793940,0.278934916674034,1.02833137114091,-0.364811252803945;-0.482291573518108,-0.306459476846105,0.904104813778648,0.372639473493584,0.430950990086550,1.11757121598589,0.237054737587337,0.399826085275832,0.0449925571571473,0.545272228602859,0.0903601981377168,-0.447461704190373;-0.00778275645521855,-1.45155470236258,0.404445175637104,-0.0122406578861963,-0.246936153590077,0.297576290438602,0.385358437646157,0.0562566459089243,0.298492814844885,-0.100109690560397,0.0192753937507154,-0.508980405497692;0.0708483823543000,0.679431564019634,-0.0374911071033784,0.558075034620341,0.650288811724109,-0.308287464937888,0.198298669576973,-0.242282223245438,0.0712943810416038,-0.0853717522882110,0.630678841288624,0.578323342875029;-0.293746996927999,-0.539197010504972,-0.0176454379752622,0.546675047275798,0.878092958953537,-0.180140314197800,-0.206611694300950,0.405518208230379,0.951286783931328,0.333743462046161,0.443424151389575,-1.93655138902933;-0.323868938947916,1.52091385794797,0.00359993362639841,-0.562260297511603,0.510711012791487,-0.193936032904075,0.246460386414325,0.168937404082586,0.179148281038309,-0.0576660925469337,-0.0670018577894798,-0.0555185220574919;0.0472707989282000,0.526461961316754,-0.291989387956511,-0.435451097400704,-0.113875293485224,-0.342381707663237,-0.501885606414713,0.368919821747091,-1.35035049574879,0.0432685561168679,0.244349240937246,0.816994981887580;0.0455809171228895,0.0286132518259361,0.101242257262245,-0.0620086799844930,-0.168021893934467,0.703183589793034,0.478926867224013,-0.113903570800595,0.632686369241444,1.04836676336370,-1.40422766481435,1.23411238293908;0.111746749848661,0.823912253635385,-0.751424739613707,-0.0150356106797851,0.316774482973420,-0.351186566041594,0.0613898060920142,0.281026077169675,-0.127852158210449,0.403221255493799,1.83611276407487,-1.91891288174048;0.216669950538401,-0.0947180951224102,0.00348858396239879,-0.0833164085943701,0.338480452242025,0.185389121098766,0.166908956919497,0.321287124293016,0.410119907557327,0.227347158063463,0.0102138160736046,-0.0974078871092995;1.59456562086750,1.14766126629748,-0.985176400941389,-0.623003901152899,-0.348534138221335,1.39299761233904,-0.212649334168902,0.0514310871255896,-0.470428566133334,-0.952904215198964,2.47579256823430,-1.54403010185367;0.0682330624425569,-0.250908748339231,0.0211948561107892,0.463974254627568,-1.08489006255149,0.507837254174027,1.02828782578702,-1.02608059013806,0.339888564061764,0.465580943063334,0.0805895448278148,0.817422560545542;0.246739128524446,-0.324841161886299,-0.218371612886608,-0.133599210129130,0.217064812608969,0.196694697534593,0.557044988004166,-0.738038994317216,0.975675159988182,0.249060332581483,0.0154740377787163,0.641374859528963;0.000783172010221786,0.200901345395985,0.138661157007417,0.223096462648402,-0.0877119197270163,-0.476416229859359,-0.388572710395637,-0.511310697204375,-0.287221984090907,-0.247020441039043,0.320410799426045,-0.0377759638735545;0.251474211082944,0.172870084255955,0.293786907993792,-0.349548881198063,0.0368082149577537,-1.52121943383583,-0.0999959808370252,-0.445004892015179,0.0354621784403199,1.51302731411057,-0.895611012666316,-0.208875400519469;0.0223387141203014,0.293305744852053,-0.322214873849688,0.0321541264665773,-0.236840779649877,0.154216106589076,-0.149683292319339,-0.0235989420800076,0.0693304276593291,0.205506443495841,0.802694189605131,-0.756200255992454;0.585744965823784,-0.339376614321744,0.390835541003623,-0.0872417501866315,0.887564205255054,0.0197013539288597,0.0347569550395793,0.259961652805342,-0.223836246816068,-0.335666253177488,-1.55829486074545,0.0199677242283292;-0.114188718873077,0.233875107529580,0.164800949207085,-0.568830521235087,-1.11895764302576,0.200377008646898,0.339801496062381,-0.102287515445019,0.00555844347683015,-0.149601704490867,0.773618974213796,-0.0136521951008524;-0.247773985367982,0.476765039832705,0.0668154552695562,0.722751287058438,-1.10942235005200,1.21651726492194,0.00284269285239985,0.770836201023744,-0.0937799288288730,0.780224012373207,0.208219767468308,0.257429162524789;0.155760406614981,0.148286675790511,0.362249065842381,0.177968408122438,-0.0468640641897116,3.87682593405952,0.447466909092386,0.103550671241554,0.497805473531824,-1.24180107021438,-0.728271410724847,1.45745920308709;0.0303837096993298,-0.221557367245562,-0.213353429805119,-0.717265404971832,-0.0370991462144671,-0.430106404511425,0.661589268989635,-0.123328254895358,0.760852718076459,0.334825543948265,-0.0783662212201119,-0.373675142859366;0.180884870960393,0.424216912221615,-0.0648211394686395,-0.399274968139489,-0.0563519638140482,-0.0943155367756875,0.363443451477812,-0.0401044643630883,-0.558970325266298,0.108269229382033,0.265584072215568,-0.0436989383749410;0.134055845044821,-0.299274773696715,0.839822544856279,-0.0371465168548265,-0.263189345376780,0.148962951351133,-0.108384310474892,-0.344339325232287,-0.0367203560365032,-0.495977306225734,-2.06501651333183,2.24002768418190;-0.0143240022900006,0.202169260250181,0.0953022167953949,0.533366009252235,0.657085161172287,0.669707615690063,-0.462072389327760,0.424329260203489,-0.937075343664642,-0.0875428403978803,-0.218524723436556,1.02988717678306;0.0601059073827099,0.0157581347625972,0.000982354722160783,-0.0468819071741546,-0.0735858447256298,-0.0324465561383037,-0.0166905104242377,-0.0461240705493237,-0.0797498772959357,-0.0233408700318798,0.0503502717856867,0.0543434238177404;0.151329166450202,0.177932376961270,0.392961242688410,-0.450031299463117,0.0952623032088320,0.120015891361673,0.290244681785456,-0.198479205572808,0.0589489771480043,-0.235776371659747,0.318634844797314,0.232734017541910;0.169100335711517,-0.00427781994723624,0.0685699236705124,0.159621276275626,0.243765164066876,0.225093846394755,0.672470789593186,0.461483395945085,0.873403708239980,0.00761337356864847,-0.267899331446941,0.210645502164261];

b1 = [0.351151024591289;-0.0982819177677183;-0.638417037787345;-0.464516310330091;0.600215540062348;-0.0488923604485325;0.123498928088698;-1.02490776965940;0.770525188490601;0.178394489359277;-0.142794766884029;-0.818247818455707;-0.214499682226726;0.151429525846773;0.0280229831579833;-1.85468164908251;-0.295463297920975;0.642785615396938;-0.202857811021693;-1.00030003674120;-0.174147007144586;-0.423391182582716;-0.148474892151319;1.15774420541288;-0.0426444802945735;-0.507140124658506;-2.04246022157947;-0.401388639764450;-0.328958790316062;-0.411120977849070;-0.108324995504607;-1.31218948837264;-0.601890059627060;-1.91859384638938;0.109687732951948;0.0567093629774832;-1.03820420097758;-0.458950978976912;0.445838688381597;0.0980018438366172;0.0378862792129096;-0.413074374709113;-0.779956393527917;-0.221649832627283;-1.04107759214285;-0.723004725307409;0.594822365837245;-0.591821643080575;-0.732251305715972;0.0336774388508363;0.416550584488481;-0.0317576405143338;0.783033323699929;2.73223439185437;1.55425260870868;0.0509662881541038;-0.0108191576188274;-0.663390859111087;0.868139495336997;-0.185040812170129;0.965464745302414;-0.0356516896164261;-0.0532943093386902;1.69712126609779;-0.141758958146178;-0.292156989579841;0.0693979983741757;0.557773345637896;-0.234193164547966;-0.532979122056770;-0.469915765058255;0.480810729255207;2.40615207558524;-0.00393592585726609;-0.373662566059675;0.216727059166597;0.609302822867072;-0.0430573091789415;-0.0353451864770696;-0.221415954779223];

W2 = [0.759624974120253,-0.347554300114251,0.257059889152475,0.716773022458365,-0.812700752485297,-0.421826541470497,0.458668652360333,1.18818038148915,0.641379348506539,-0.300298513171394,0.716987718040251,0.922205291198515,0.862531328050250,0.454408013454080,-0.168104630482520,-0.899403058692250,0.904129260396135,0.545478189286644,-0.480031613057542,-0.169813193227539,-0.845757209023379,-0.641483527689968,0.646302614931512,1.07364406215485,0.260617807304343,0.306591119276715,0.853880044873074,0.602431395414328,0.468268452283741,-0.450687741902564,-0.383756865570479,-0.629494904102675,0.743085446148503,1.78140979797515,-0.638907441660075,0.317312960323567,-0.946115872074272,0.894305687330815,0.825621799257678,-0.315541656846077,-0.267208353068439,-0.366720892099467,0.612405512432426,0.193350565713716,0.652210652357626,-1.05747881871624,0.569515294739746,0.846438076087544,-0.857386347040950,-0.592756276816664,-0.771299336658460,-0.349934519622143,-0.910449652929047,-0.563843485915301,0.992242221360494,0.855588676178796,0.359611666026200,-0.597402213532316,-0.389394252437949,-0.914475896329588,1.48287958928104,-1.27689300814301,-0.335955089743830,-0.681444169776636,-0.359991033998039,-1.26114672729866,0.247291637460755,-0.819692644631012,0.309576014821907,-0.828319960170149,-0.906618236028593,0.488580011266871,1.46609868205525,-0.883977037004226,-0.499084092409138,-1.68146123087171,-0.935587264476298,-0.0180854114084687,-0.530629675143897,0.162294215827411;0.720211766446495,-0.212004533127008,0.258958903650564,0.731655984462033,-0.745210290462153,-0.00850798081054877,-0.273202301539329,1.13287162043917,0.679636340806351,0.235985801093618,0.756030234279590,0.933857107556840,0.868032068026779,0.248591584477440,-0.421507747879470,-0.853111336957196,0.872362434656134,0.528453510181329,-1.10036010133200,-0.174796841587673,-0.671141460423429,-0.613011777224312,0.667252974436740,1.01816420555067,0.417298478775488,0.310156322661921,0.853740003532538,0.696749656277630,0.529158825235922,-0.493305999816591,-0.367929765473501,-0.565227226384214,0.704450991203772,1.64975863262696,-0.589153859392535,0.602031790084040,-0.928043565856389,0.957702441926143,0.820458200054964,-0.305159505035623,-0.267709685963402,-0.381217120573198,0.630277049010878,0.186830415252075,0.648231759399301,-0.938427387618187,0.629798335437555,0.908297238075383,-0.834554234702591,-0.594391470083225,-0.774959075525064,0.403825152969787,-0.934539576252013,-0.572041420402612,0.964156796209998,0.855051460908001,0.362749096992536,-0.581875707344630,-0.402243502187007,-0.986727387569540,1.38575900763069,-1.29727860103548,-0.594846223631821,-0.638648766745994,-0.332278551799596,-1.30834875068233,-0.127541877289321,-0.800297538016563,-0.439581662995194,-0.848285019981329,-0.865387107365151,0.541824056819044,1.39957217517484,-0.892495288560173,-0.648755943299586,-1.66471497421865,-0.938606783206182,0.0748028372906397,-0.568357780046914,0.334502214199711];

b2 = [-0.0598920371290885;-0.0329032508119668];

return


%------------------------------------------------------------------------------

function [ W1, b1, W2, b2, ti, to ]   = nn_6GHz


ti.gain = [2.16592611029879;0.00211822121778521;0.00585936801508902;5.57009709446998;0.00575559508682464;0.0603549007618636;0.0273979188601703;0.0356645790313325;0.0416423621568934;581.975479994452;5.10112142957262;3.30615911264754];

ti.xoffset = [7.62939453125000e-06;0.344390869140625;0;0.519999742507935;102.512039184570;0.0200001131743193;-41.5389709472656;217.145431518555;227.137084960938;-1.47336226552497e-10;0.671546816825867;0.424632817506790];

ti.ymin  = -1.0000000e+00;

to.gain = [4.28229936376372;2.82166288357301];

to.xoffset = [0.625062048435211;0.332074016332626];

to.ymin  = -1.0000000e+00;


W1 = [-0.00748418194958163,-0.0122889940168742,-0.0573015690826989,0.0359212373832173,-0.0196013008679231,0.172745872812107,0.117271084640953,0.0563121041628669,0.0293476816343474,0.00662309608704419,-0.0448179462830128,-0.0258512556086603;0.0440762101149495,0.00738650105507893,-0.00197524420708961,0.117574170449546,-0.0295969588615848,0.225605532817540,0.176058755536309,0.0389656667904867,0.0125109074394475,-0.0649877717746503,-0.0658922714105790,-0.0204811115511894;-0.00158539257238627,0.0142483849757931,0.0537395209870630,-0.0604117262746759,0.0411441182773901,-0.198590620842057,-0.130526380420360,-0.0648032627509076,-0.0299880443571919,0.0105054491191717,0.0486893048069432,0.0383466339017280;-0.0565356158403608,-0.381984322754830,-0.140349697626256,-0.100099673127798,-0.106318781513376,0.0464153939177422,0.216595450847258,0.595101768898774,-0.0686886333438555,0.00921771366135762,0.544775088083000,-0.0385929404296133;-0.00602860296890235,-0.370479291479695,-0.154602806480802,-0.0676374308373209,0.0374907689686089,0.0759935879372123,0.205109332466282,0.0668475009699627,0.0745137142926220,-0.303842668106552,-0.241252024583662,-0.208427944002889;1.79644238030400,-1.10436235233903,-0.811814709329312,-0.866668961299448,-0.0412469934486035,-0.125781761759466,-0.0455179741805998,-0.205151615135600,-0.418754046225832,-0.969640038013590,-0.710836273122973,0.662890543661295;1.27652162241770,0.807934512140170,-0.312648168053851,-0.594107042184300,0.508309752529516,-0.330187628109367,0.316331569368092,0.467630637837109,0.322040022147130,0.0749994821625872,1.03573306919592,0.765511797157630;0.0203583153059176,-0.0331163874826731,-0.0895646734687560,0.0347195302988850,-0.0809123609108026,0.200651127790238,0.131781739031822,0.119505205137345,0.0250682111365397,-0.0183556578924887,-0.0557927675978531,-0.0999793685990736;-0.122133225946372,-0.923643128408784,-0.401064631241204,0.642615044308203,0.932792693758242,0.416782059515342,0.389504893050476,-0.855597517476993,0.544744904874498,-0.0514344534131833,1.85339896664036,-1.38459999855839;0.0519700284332436,-0.106030884796580,-0.0822631684274468,0.157650180102375,-0.198686362412821,0.312095483739172,0.126721123387453,0.128994747642272,0.0638215329575267,-0.0971531272977078,-0.0397886000735347,-0.145128572289066;0.445236556916952,0.214352415896446,0.711716876544682,-0.444216793646240,-1.14828295864696,-1.54226820822260,0.489325876853252,-0.0261162978413669,0.142729500638123,0.0125980922381160,-0.478643301595219,-0.159383729716086;0.281977038377343,-0.0348758868367972,-0.0234275836038523,-0.267822405424285,0.118473919412867,-0.716590508297200,-0.105969456823858,-0.518253767521903,-0.199992824399016,0.451656817959347,0.588425096872529,0.190263955292165;-0.227254416711449,-0.690143995416526,0.229434375673629,0.0537122191470543,0.653380398988707,0.258826323541514,0.0719391976108545,0.253232552012014,-0.198426660417232,0.0992735277224522,-0.0837194131404986,0.293361226116650;0.00157810766768200,0.00466180862920295,0.0514579220223312,-0.0452753158140800,0.0360921681117444,-0.187725993460217,-0.131773818786348,-0.0636706975456847,-0.0222709745882940,0.00842400749964647,0.0497561726982904,0.0448133441367838;-0.0165632308555647,0.137502348604296,0.265707716670102,0.00857923612201662,-0.124074066933502,-0.144717367576546,-0.129731052028622,-0.200950927728583,-0.135998063190638,-0.0141555746080843,0.144198304398755,0.0966105293146560;-0.0103766520337534,-0.219107518526037,-0.441292755599365,-0.249551746332649,-0.427165386662743,0.0956220389902855,-0.608028864074209,-0.132519179648117,0.178908965159914,0.209576224197851,-0.0928288778684100,-0.464144769050251;-0.0435377778935408,-0.615229549892746,0.159251644896017,-0.0847942477888633,0.273142812970258,0.676677085349340,0.211881583152138,-0.799761079361001,0.0222693455260790,-0.583875973343898,-0.206907141103404,0.343543310186779;0.0474789661181813,0.0854614835901648,0.117688174487104,-0.0721372952913180,0.188236584664873,0.00480928333008310,0.216513097187799,0.0764697499860946,-0.0319376080656685,-0.0813623579418968,-0.118596287944203,0.125686703777870;0.151437276528926,0.0144351879755797,0.144695462750661,-0.0354569669006034,-0.0169639346356299,0.0243068327359130,0.0594358693261829,-0.128168355548114,-0.0154873439055506,0.0823987242116253,-0.0654672275804083,0.101331986265641;-0.465968520285545,-0.157090919966337,0.0248139487229466,0.186918827543612,-0.506164267306708,0.233301928333772,0.194126297856987,0.468792727202935,0.172130029177731,0.674112584944127,-0.193340063892501,-0.0280961972814805;0.249958172467391,-1.41088104974133,-0.00959231999341658,0.492867020198617,0.568044837478214,-0.121819093609528,-0.579682908930405,-0.0445846108662332,0.271673332541226,-1.03823174004629,-0.494062176190926,-1.33315375834079;-0.136418333750174,0.912201827714173,-0.119650403539871,-0.222091926350311,-0.382438193517203,-0.311586176269777,1.03113927061379,-0.627636538402228,-0.268462462667269,0.995647375359666,0.555402425867575,1.68202340804214;0.145948567574408,-0.00635747192841746,0.449040416275712,-0.119878735163558,-0.732836298758370,0.382174746920216,0.499850308082596,-0.239396044488613,0.355907423303678,0.0729183892866797,0.482785826523893,0.0983101862899200;-0.00303546804716815,0.0478566575840722,0.0704335149916217,-0.0979856767580242,0.100967129931331,-0.257892161743153,-0.119325322273428,-0.0893977914907757,-0.0505473719314769,0.0362566858048304,0.0521057146827894,0.0694541295695684;0.221479878104170,-0.238843368006186,0.203294828422257,-0.326491551272527,-0.328039717883192,-2.44855658287994,-0.292366836182205,-0.0956557433113928,-0.0902539295276490,0.200682650312599,-0.134987322947171,0.337044138283272;0.141655018470191,0.137537247478105,-0.291256875244404,-0.0324443376273801,0.636808675169325,-0.279237068092131,-0.606362746652521,0.338110990405368,-0.407587591631851,0.0346155676415042,-1.07084936130881,0.279274285295051;0.000819751448772691,0.0578300698117720,0.0410652788282698,-0.140960330921277,0.0927204946190112,-0.332616759931560,-0.106929553519740,-0.0713679724027927,-0.0557073170964690,0.0587581044453765,0.0475293120038070,0.0232637678700580;0.0114584130672068,0.0324386010256101,0.106105294413816,-0.0168777054465297,0.0323962757425437,-0.190031458542556,-0.102458746233290,-0.0935959461876787,-0.0472772673519758,-0.0123070018944771,0.0449121682820584,0.0487179583062589;0.403213218191377,0.0270716039281458,-0.126789730207197,0.780043573716731,-0.181317229049752,0.148977957887803,0.0566642128859707,0.192948760056852,-0.381952034867739,0.169094142317364,0.0160434617971326,0.0428646397273578;-0.353803239764637,-0.140195920734476,-0.677853506340123,0.0245271598437540,0.611319549794500,0.278424330339103,-0.410162660472336,0.579025313868003,-0.112390944223473,-0.526183962717524,0.447483260483835,-0.247203176650885;0.114846639175873,-0.311630620147522,0.305649327210190,-0.108617770844363,-0.116773564112534,0.219737750698754,0.632200526431747,0.527435189703375,0.388578654374889,0.671188536112058,0.551971213166979,-0.125708209870903;0.0340480354870143,1.30442015571255,0.389682869574644,-0.160123474856901,0.248348722668967,-0.487339166322515,-1.17229402152326,-0.833815795048590,-0.544237648592147,-0.0619661753483811,-1.76684599514306,1.12441463979002;-0.524352948422244,-0.153506696668955,0.497730195056462,0.287782498683558,-0.0396344631222767,0.365843153084320,0.415784810757662,-0.119309559593211,0.522900977129655,0.267218712791242,0.182667756487177,0.532471031600094;-0.208052882422496,0.549116893445157,0.804011218802970,0.408998158895705,0.172241589357008,-0.0535196979466398,-0.0747690585968361,0.467411893600353,-0.0490081966357887,0.190978023788948,-0.458743223884645,-1.15885675452314;-0.0486273279914336,-0.00707653681035470,0.601840499325398,-0.0164030095623402,0.341748135517810,0.134254800645327,-0.288845733522318,-0.191370233941293,-0.139296797036100,-0.154465015260539,-0.118334173257577,-0.380752021919975;0.162409405646025,-0.325383356795862,0.0857208802760898,0.410014661034816,-0.508755543003470,0.0881807943213022,0.413046560142099,0.192316533605592,0.445519861311082,-0.315204166698200,-0.613274321008151,-0.198349724555396;-1.27217921274789,0.233017806862264,0.442832777328926,-0.152487454881121,-0.632382311637776,-0.257040014495223,0.403861480042810,0.173106159121111,-0.581256180326087,-0.359924785298555,0.475565495002977,0.644899839880969;-0.0694416082341732,0.532020836018047,-0.972990492067465,-0.393321941657416,-0.779421780729190,2.36938071910273,0.210551866060505,0.482932197470723,-0.339871683968528,-0.697459963644663,0.577426536483979,-0.522388754908064;0.0362361884568312,0.399128356926115,-0.995509872320970,-0.358541626743267,-0.629416122893147,1.16948872663439,0.125265048132648,0.292045897341869,-0.446560266253055,-0.895187763966801,0.712489708288948,-0.889384266505466;1.44629097350536,0.796646253901832,-0.275955904976873,-0.641628352816329,1.10007354192271,-0.591783187404760,0.0577542422708670,0.815645636706012,0.0954489985440024,-0.855761453391330,0.852493640976411,0.625898585785858;0.365633269737328,1.06231814431123,0.0775609217305410,-0.218898357001556,-0.0123054459300737,-0.208231848673909,0.757811868438041,0.383966387031480,0.507803893416061,0.467661437333121,2.01109666440575,1.01457053372523;-0.366828422656580,-1.21861644252080,-0.251294848267764,0.0784468244390063,-0.298301156330895,0.511733877603465,-0.650690576116383,-0.127510446665874,-0.288011914639882,0.00436363776990252,-1.58141000637069,-0.568998570121937;0.0704595303155760,0.0246952254009303,-0.00943083399714493,0.0875431085567886,-0.0795487967024369,0.186160252473068,0.187907659350217,0.0567291210058084,-0.0245295661672796,-0.0876768551107165,-0.0955748433085708,-0.101775339916400;-0.364715684838000,0.740051867655426,-0.172643711620897,-0.146450658307248,-0.161170004624096,-0.447252773344587,1.12786299944459,-0.897193731132057,-0.144148149395743,0.705012243794163,0.968533376271107,1.72154796223180;-0.682374305440085,-1.30243904945232,0.358581958295103,0.136243829578566,-1.22691227232841,0.844958353871653,0.424956923723204,-0.192193622200180,0.374552901396404,-0.822372280322044,-0.318046454323647,-0.752413970135932;0.164060898096519,-0.197924422949590,-0.230716296194127,0.393783725223567,-0.383692915447854,0.183815529297276,-0.227310305822740,-0.346276359954655,-0.372432486937174,-0.649295588661277,-0.0613695652994797,-0.127264306318861;-0.336769785306213,-0.880373074641213,0.324749000801268,0.389215992804102,-1.48273521148810,0.742090144662917,1.32506601160678,-0.423455486010222,0.687125738969002,-0.499880727083328,-0.417247595337103,-0.541765915375200;-0.0730106177387742,0.00655412711222418,0.157191662310131,0.242266059785143,-0.0539212627897957,0.0947637531726748,0.337085856282511,0.218832068837492,-0.290659623736718,-0.530755356037340,0.218741786041431,0.210710279112068;-0.640957132903007,-0.106411818267187,0.613227495733601,0.455308490245906,0.107261528006772,0.136213887712009,0.697873614847224,-0.282288368052983,0.550549927151765,0.575109653380208,-0.170289226762430,1.15562464679296;-0.423305785589975,0.0722633560107776,-0.589461816044871,-0.140449227064026,-0.449097179332946,0.219384178906785,-0.593420487890086,-0.617021705851044,0.621913396008537,-0.340765854300848,-0.344521602586424,-0.748274601466279;0.0892386000434068,-0.0738938896355948,-0.225692917377982,-0.155673427231722,0.178821038952342,-0.130229665537791,-0.337356171940115,0.198181330169608,0.395568064879579,0.275047687170962,0.245093802072399,-0.349936005023781;0.130375470215294,-0.248969928743927,0.267486032879767,0.0262056190298326,0.186824391855175,-0.0461360988699617,0.671917325570916,0.129888372579858,0.371503615095105,-0.181314527763816,1.22994576064375,-1.06983893177647;0.0279850897741101,0.0273475060121780,0.00432154909593053,0.0344768891736719,-0.00504487120473311,0.117251809002767,0.123793985881313,0.0237878988212742,-0.0103958309221032,-0.00995492441081868,-0.0745356617195199,-0.0390892436881404;0.358853177154602,-0.189468125936462,-0.0130036833903142,0.203138660647750,-0.277134439787101,0.181893234847686,0.265595117966750,0.0883071819864578,-0.138648719379775,-0.0213526653660220,0.00672531667204086,-0.230300504155421;0.0320018266892522,-0.0455073779899978,-0.260790444944222,0.274855199875759,-0.394398356606407,0.238377021412419,-0.0370597407883401,-0.568647191184840,-0.0379284897082036,-0.340594580015959,0.257046232919673,-0.209203039371011;-1.13636196948822,0.350030234576420,0.656009913899495,-0.191237292110841,-0.786446876242795,-0.103001917842226,0.844093107275813,-0.0599769695458104,-0.537816326015799,-0.345474354457921,0.965816868473230,0.170547286427314;0.715160025076927,0.167733565972115,-0.355303488275126,-0.0385228179062321,0.341916526468342,0.780721929649138,0.133558988640124,-0.673858234678481,-0.123154536309332,-0.600127045709281,0.250667410487334,-0.486923179143445;-0.128710285577484,-0.512781264402552,-0.639125380316058,0.746328534139271,1.12525287591957,0.428716528625771,0.104795040797569,-0.671016496323974,0.436013335274501,0.137351270502543,1.36244985047397,-1.45236382803710;-0.0654680871868213,0.0821919949230842,1.21059039637280,-0.0882961954800848,0.174954011941995,-0.0823069015066024,0.104445791791512,0.0148538559904796,0.0954592264823080,-0.360646976006154,-0.335377683102424,0.367183109915703;0.300747176517836,-0.563204677513230,0.131910678644902,-0.593331594027794,-0.527247524323436,0.0180282448857419,0.0743921725519339,-0.364773432409395,0.390044475954123,0.0699990146799211,-0.802184265787031,-0.185291301398768;-0.0118664830684181,0.337552741834836,-0.0950739037176729,-0.668624261798803,0.234848259731247,-0.444414582533776,-0.269235989281581,-0.0523324981656687,-0.358660325041685,0.174323420720700,-0.00593116819282353,0.697614789923860;0.129222360559874,-0.361298909189180,0.577802373682580,-0.477486534228260,0.138919851751223,0.0831060703113749,0.374436000663026,0.241097977937428,0.0496352963587169,-0.291413032543767,-0.342100475646902,0.391437724923827;0.0439945506639984,-0.0909731877035138,-0.558802340616674,0.375381190778241,0.896920004827469,-0.173352795859621,-0.405621392534401,-0.251721991423685,0.00539896327927944,0.409767003875469,-0.555925415028548,1.07399786540574;-0.0195141137475495,-0.239384585672568,0.398778171656269,0.170799988107522,-0.00632365672931945,-0.630418253402896,0.0927143384091549,-0.234497110261246,-0.0977881169120816,-0.132536500510739,-0.0853267004394215,-0.203531267735433;0.199009937493331,-1.58706417323339,-0.448115471106877,0.506384659106775,0.792047492913999,0.704666482722129,-0.206105809557764,-0.707685228165874,0.230314817678050,0.248441918432818,0.289372290355280,0.279422759575958;0.0444925880561716,0.0291232139515630,0.0805919383427586,-0.00826080484819908,0.0113393197067062,-0.156027602437458,-0.0769445990447117,-0.0638559602225060,-0.0339506087684128,-0.0249972770621121,0.0368445407057878,0.0214221211361281;0.00457186510203193,0.00601806055200381,-0.0770013476569489,0.0520552292276531,-0.0392440851024606,0.142051133497896,0.0799086272896774,0.0418553164201327,0.0472237397630586,0.0158951891381606,-0.0798402240825489,0.00612414163353300;-0.507446155065232,0.185337971366657,-0.642252286054362,-0.253181183488152,0.494530634970832,-0.484379374896827,-1.33266326718296,0.119914896877426,-0.402800186579688,-0.0699340269321402,0.413981287336249,-0.963781161504405;0.435185991786430,-0.545455084940596,0.243493826299695,-0.0290113481706347,-0.651881242187510,0.254514924569500,0.476768031824949,-0.483171034890765,0.228513891962898,0.145712780162915,-0.815737373811669,2.08230188658208;-0.00253719937203154,0.194236650755505,-0.371479187455433,-0.432473681539015,-0.135182179875796,-0.829719720495862,-0.226198922110073,0.158063875271432,-0.361503120869780,-0.203415503563628,0.302321399901450,0.0704186301865985;-0.442075763541886,-0.702209788514092,0.0818384298166035,-0.0420592284694066,0.301560009759737,-0.350121716405043,0.523216823584425,-0.119189321394151,0.504450389985902,-0.305156524513346,0.613264932215600,-0.739590490843901;0.00774402504675278,0.376374774214121,0.0654537748992363,-0.380951495458593,0.432345123902419,-0.0727832435460474,0.0382061169197539,0.0445030879832720,0.136397650318096,-0.460567065336822,-0.0401786174847528,-0.0671196559002678;-0.549081288928107,-0.156498246068764,-0.713876885237030,-0.134512649385808,0.247898942282269,-0.606546067329613,-0.829859083364077,-0.0977472849826671,-0.177975495507334,-0.119163538673821,-0.00414232672278338,-0.0363721855343699;-0.000836110581809252,0.0431782818624396,0.0485523395667369,-0.113915056458637,0.0950276454853385,-0.265176633585460,-0.126707880245883,-0.0707080072824243,-0.0457910312707612,0.0377627298460085,0.0507482473563878,0.0555477135594086;0.0241749754870944,-0.598785295360565,0.522282630960679,-0.297323557448275,0.222148084105122,-0.145210776638974,0.305908378059443,-0.0955724979588649,0.0442728213534451,-0.600412574480812,-0.327751422325622,0.274676933589366;-0.189694964971675,-0.395719505229817,-0.201613163085941,-0.00393209494314911,-0.301277208366603,0.315802747162595,0.293507135388624,0.719505240942428,-0.0315035646840634,0.163549299788589,0.225493668825551,-0.0938764260575545;0.0960471671834172,0.483135286864049,-0.511478297044275,0.349172841940136,-0.0689971869378897,0.138040504537254,-0.806884215051471,0.00635119995885461,-0.396496397176212,-0.131843682342907,-0.449072103470075,0.267403509829019;-0.0350492672293819,-0.00993529450118840,-0.112753104131331,0.00421555402970759,-0.0608292508035396,-0.0411793704953027,-0.131064533572670,-0.00381746444453154,0.0454405390319201,-0.00912578781773400,0.0116831306229096,0.000289357735041703;-0.0147398060629551,-0.121623259440455,-0.0119303029949572,0.232365848732582,-0.167735117048717,0.362231937052671,0.0723072674129739,0.0622113220732517,0.0923991511255337,-0.0569511676600073,-0.0540188460115756,-0.0537708465881842;0.266987121863963,-0.115250727742771,-0.196545474113703,0.382948791328146,0.204899525371221,-0.116470979834795,0.304590947006779,0.152876359701797,-0.437643488322341,0.551144912077138,-0.0623190414475209,1.24810702387010];

b1 = [0.00749698553331446;0.0336776920592205;-0.0112645801717237;0.165272356935840;-0.0492126007298417;-0.492992431345173;0.649006713070275;0.00174264556495425;-1.07858249116783;0.0145175284360341;-0.612879729940965;-0.292116666506120;-0.120785014247103;-0.0105081718770426;-0.0908149006065783;-0.0113168013093966;0.280452454099749;0.0437861522021422;-0.0295375347569565;-0.271656887313428;-0.540184118079917;0.208492791358284;0.0256680482021213;-0.0191953734973674;-1.74864517937162;0.570216830438584;-0.0706639817871553;-0.00891654991730196;-0.0572099419236944;-1.17026978540539;-0.111700923715879;2.65106216375627;-0.248072530775888;0.428653599516093;-0.0716474585827338;-0.189629250877701;-1.52222617854136;1.10830021353438;0.339161970164969;0.176450881768113;-0.134455208880994;-0.0285162440188893;0.00114833500929072;-0.568736087047081;-1.80960251006182;-0.0193488109565949;-1.54248151977345;-0.0766061898337799;-0.336005488353878;-0.257937449873990;-0.279949118054733;-0.541626766544303;0.00163113309115775;-0.111060986641424;-0.0761792320213938;-1.00626435192402;0.471071785572651;-0.639483589631727;0.489402823580647;-0.0733293425045488;0.523553413999672;0.235629150650772;-0.515245747398421;0.185944119013041;-0.981532692546515;-0.000862954023966287;-0.00720956679306578;0.235640258145830;-0.817712189601703;-0.565342520715108;-1.12402976474815;0.0639849546462601;-0.574538863736907;-0.0213604177007284;0.144943210220840;-0.342709375858265;0.192552111459633;-0.00989362894957349;0.0362239145871840;-0.146268733674436];

W2 = [-0.226216452999212,-0.355085404481171,0.261354287994357,0.517159549393181,0.263491496019804,0.140157292231535,-0.752049569845347,-0.269370869262420,0.797092716426971,-0.332687361087034,0.567203965871403,0.628316155470220,-0.401008267524573,0.252914149218832,0.298163652510797,0.713663865835516,-0.659269004306999,-0.457961515015303,-0.288014116882520,-0.137971353968642,-0.347543531568676,-0.850418294412639,-0.940266284241124,0.294410325610709,-0.765917488477996,-0.554173410518682,0.343343045084249,0.220957140334486,0.341395235507898,0.197886447440896,-0.0154635092847950,-0.331596230715255,-0.658007614080986,-0.648322322225180,0.452409681094508,0.494246951355326,-0.811721088375397,-0.621940497480861,0.665355351553816,0.464836366183527,0.509635924026149,0.714400182773480,-0.342559116774867,0.551126236673630,-0.493242029610630,-0.00247397181948832,0.496576698567952,0.698178175007076,0.637260657593379,-0.617467350158809,0.0247637044473324,-1.09136435116082,-0.233676632146940,-0.346645611919541,0.184184207179759,0.699413528355503,0.519908121606742,-0.745441826561447,0.394277750193403,-0.590174186392112,0.834594203903474,-0.321521678628942,0.639587491270929,0.0147124134320422,-0.409388438167005,0.154650548251770,-0.198752057333195,0.849618388511020,0.619239905832753,-1.27341440240017,0.719982517589253,0.412742942833701,-1.48092307410671,0.306737306836622,-0.411808506418107,0.903930986989978,-0.214885513581397,0.192986712989950,-0.330890848480658,-0.724917559623244;-0.144403875207255,-0.100711582273976,0.159943162950600,0.523002024630542,-0.325294560855382,0.109401273462619,-0.774415474192382,-0.202123198653698,0.968462595400520,-0.342042870915492,0.535691964534309,0.472474505955250,-0.414216460297161,0.145970025259749,0.289355553937481,0.270255561441919,-0.568072992374355,0.119589093546561,0.106677273137421,-0.696262701833365,-0.346504582557674,-0.902775284109860,0.554278745332610,0.235110993590960,-0.589656041455599,0.833200704312391,0.240422371031959,0.206567030015649,0.144517434253788,0.622743661783432,0.542497149859634,-0.299608387538482,0.109055026312850,-0.437543842965981,-0.185253908264901,0.443752686504423,-0.741913080237263,-0.451807431015822,0.519251660581812,0.455930704931630,0.462302975622994,0.613260753640660,-0.0699195780251167,0.577396584287477,-0.487214349735179,-0.455848855612211,0.431544080804898,0.0175534941006059,0.219592332670233,-0.456071021886118,-0.619615980666307,-1.34011400233838,-0.0297409245875398,-0.428942188830045,0.767126468494801,0.676918145471251,0.398694752290294,-0.907502307877207,0.0549221150415880,-0.458973942411847,0.498451060084504,0.623987349497179,0.776136680754444,0.664872755846278,-0.449154604802458,0.166106609695395,-0.157887495236988,0.595811981060683,0.559197926517654,-0.979785210784296,1.20969410126877,-0.123756313172087,-0.837878788775827,0.213977809852417,-0.841936520611238,0.935827169159119,0.649751342132071,-0.103490044534155,-0.313130702328991,-0.769017268522533];

b2 =  [0.292384627891248;0.291143032520085];

return


%------------------------------------------------------------------------

function [ W1, b1, W2, b2, ti, to ]   = nn_10GHz


ti.gain = [2.16592611029879;0.00211822121778521;0.00585936801508902;5.57009709446998;0.00575559508682464;0.0603549007618636;0.0273979188601703;0.0356645790313325;0.0416423621568934;581.975479994452;5.46565141145039;3.47287193983234];

ti.xoffset = [7.62939453125000e-06;0.344390869140625;0;0.519999742507935;102.512039184570;0.0200001131743193;-41.5389709472656;217.145431518555;227.137084960938;-1.47336226552497e-10;0.697082877159119;0.446913689374924];

ti.ymin  = -1.0000000e+00;

to.gain = [4.727096193774356;3.093549344883394];

to.xoffset = [0.665690064430237;0.381629198789597];

to.ymin  = -1.0000000e+00;


W1 = [0.0431847048815179,-0.000677421754823157,0.0474188675442228,0.0193889209811275,-0.0906685412105355,0.144633304479408,0.0930785567568879,0.164545390134599,-0.0153878017240546,-0.100076049578532,-0.0110186176305710,-0.165222460733685;0.0312366203865967,0.0292103519212177,0.0244142462833136,0.0119785477487235,1.14976105967016e-05,0.0866498894008397,0.00534457456532965,0.0402083981766036,-0.00903398707929249,-0.0265770838050115,0.0186808434262724,-0.0575635044090157;0.00376998475408470,-0.0289330456174543,0.00630004319621914,-0.00222627436471691,0.0189914287854469,0.0326917428980144,0.0239473551514276,-0.00277221233571390,-0.0186723030861316,-0.00541764838551586,0.0228170287366871,-0.0344660857633590;-0.0418396373930056,0.000761613355338035,-0.0425541705444578,-0.00938030849547507,0.0134034889543233,-0.156647980764833,-0.0430575686975547,-0.105213149715997,0.000450608773843420,0.0920318860835478,-0.0889591615989646,0.125268989074728;-0.0358995753404967,-0.0204205796085891,-0.0255155754221466,-0.00926852686947299,-6.33085149923706e-05,-0.114691133178631,-0.0151216984491363,-0.0573176931332837,0.0102827580757821,0.0491139985037468,-0.0378469624246001,0.0882738764357516;0.377535052786668,0.674504828364397,-0.0302851936277041,-0.427312457725518,-0.0673515113960211,-0.809274318305039,0.828154183773614,0.157314594663881,-0.440495181160694,0.233263794276720,0.732679739620035,0.848513124948303;0.0315563161017454,-0.00223623052552035,0.0290375278621027,0.00860093487230151,0.00857137441785510,0.145360762658251,0.0329208576866088,0.0658584076627066,-0.0133292223171512,-0.0600267801485131,0.0851444776424383,-0.104619645454599;0.0341827079430431,0.0229182314050763,0.0235728295572883,0.0101757329003157,0.00342936860283611,0.0918137450564574,0.00782233530964108,0.0405312875902397,-0.00961872678014577,-0.0326584643682037,0.0228368520405263,-0.0657856926469337;0.0420085968875881,0.0421716271980656,0.0402626733945834,0.0145535415160365,-0.0148656949758839,0.0913343970988301,-0.00266721726688564,0.0577300710338308,-0.00139399798597295,-0.0319482785341553,0.0157311298974470,-0.0513525268971791;0.0302515763216688,0.0202179899130911,0.0192851700353832,0.00871277501805156,0.00543171950416408,0.0771162784387971,0.00597684772764138,0.0306986488985139,-0.00876380497529941,-0.0249651472005678,0.0159027497751183,-0.0555850135229631;0.101274859728715,-0.115826132237494,0.134277046593832,-0.0748427370615023,-0.206363184087961,-0.0670421711814439,0.00320054953519248,-0.106746383809893,0.108450191493816,-0.219286700234825,0.114517338972656,0.0398240378137634;-0.0583409207180924,0.0314315062965219,0.0255977494122515,-0.0115138429045848,-0.0662161197510381,0.00227141831395829,-0.0110063165206455,0.0599532466438399,0.0267955657266280,0.0204514283361409,-0.0517430613716436,0.0799127294175501;-0.0838720539638613,-0.351182610997262,-0.0331377911639667,-0.0337480014038950,0.607741371950245,-0.483777318314165,0.0221808348851406,0.456875876441576,0.960723546469035,0.869061880905852,0.193107578931175,-0.903680332008600;-0.0316468689558152,-0.0497045125496718,-0.0403328973108061,-0.0130493956363830,0.00753812041382268,-0.0755369836135023,0.00762069340352935,-0.0456042997701455,0.00431015197171147,0.00328402168110984,-0.0220388911920534,0.0211660329764265;-0.143709011849661,0.295934580397536,-0.268130779949000,0.331128475000519,-0.372072240814850,-0.594896542828097,-0.0861104375557086,-0.123283966994228,-0.379989813934028,-0.0241184661946156,-0.101719400617690,-0.654544133656738;0.491734186226385,-0.234474113717821,0.158587461114948,-0.132996316243513,-0.0434883275101884,-0.208159118175486,0.511529709156491,0.241119292172267,0.117536889002143,0.0598324697054819,-0.464716331931522,0.0969840278981980;-0.0442731791779010,-0.0877162535815648,-0.0434710645247287,-0.0159299748485233,0.0401115107684488,-0.0678437235756296,0.0552389436265024,-0.0550079260699640,-0.0230659468268969,0.0116552817664019,0.00594004775101338,0.00917763667179591;0.0327514755556038,0.0283552249903452,0.0239821486284341,0.0105741846864638,0.000535108458398660,0.0780084924985057,0.00211240745840126,0.0354810769877557,-0.00691768907689379,-0.0240867196010552,0.0140001696483399,-0.0507247664149506;0.0455822197326704,0.442401275909240,-0.705664160656979,0.0202214648799205,0.119988081030266,-0.590754855442967,0.0366684107802030,0.401108312181642,0.225771181769897,-0.0710746072544459,-0.00879465893191285,0.605787648240908;-0.0387378408944123,-0.0268600497104384,-0.0292643478302791,-0.00858027612216798,0.0155875520767061,-0.120463824683398,-0.0169228230087952,-0.0745916082223861,0.00772439459288358,0.0593344693861683,-0.0345968747781259,0.0980540056908666;0.00394242439385328,0.00922448051053359,-0.00193017837336537,-0.00236511271668528,-0.00740444854255725,-0.0179010470103204,-0.0137579324140028,-0.000926746938015395,0.00985447832097696,-0.000114283918129747,-0.0102386024796885,0.0147555591278236;-0.0375759239385781,-0.0131792934396293,-0.0332335461726812,-0.00698129251941273,0.0272775598088651,-0.140796959796695,-0.0360589156706854,-0.0995461159556208,0.00663061737285593,0.0793938122947922,-0.0548065345944911,0.126969574138957;-0.339621576637660,0.137988148568172,-0.175427811760092,0.0142341870159563,-0.390182593425542,-0.0455047786532379,0.453921574328772,-0.335368701551368,-0.447229331626418,0.543446132044456,0.154440899332805,0.825742797716766;-0.0224995737446846,-0.110966094490869,-0.0373360867460388,-0.0157820469347588,0.0293499988746132,-0.0376895734223826,0.0697948134737709,-0.0475645173729519,-0.0278673182835384,-0.0183114607597129,0.00593245360038667,-0.0386136731158939;0.162779952675555,0.127640605862281,0.0795441337720342,-0.0596067667847352,0.00172696565019380,0.0196012103669461,0.0988270212097192,-0.251508818335828,-0.759778218880031,0.0716006577606561,0.105238422880652,0.156381302438053;-0.00451913056649470,0.0568749494611868,-0.00869229715355252,0.0428651375396484,-0.0628736689073440,-0.151476373401549,-0.0508970112639459,-0.0333423165563274,0.0139733951059977,0.0615872538073391,-0.103846563631939,0.153950176025615;-0.0429953480935178,-0.00976299765334823,-0.0572451863266090,-0.0212964151458081,0.0433639351758200,-0.138096614806023,-0.0496181684862317,-0.120521834005994,0.00923272434727864,0.0818602403396511,-0.0548004850711287,0.112062900989272;-0.926392076556928,1.02848468793818,0.270672225812692,0.815018944110022,0.427900178038561,-0.326629962260442,-0.0414349722685355,0.370411963929291,0.226104843090665,-0.0467813560728240,1.51603801515016,-1.40878380972167;-0.122051287259338,0.470987779371461,-0.112317566849715,-0.113107985643676,0.216392242204569,-0.433464561062564,-0.741935270451818,-0.412180446982486,-0.196235884472610,0.470297366803455,-0.878046679693402,1.57900422985400;-0.0995257117606834,-0.123045807185451,0.0206710277865470,0.371950939387624,0.0298233322395261,-0.340282774213507,-0.0659249923773448,-0.0569884725755906,0.222184349402052,-0.298999178608269,-0.554659785941257,0.145250291961448;-0.118351824570563,1.40310916202800,0.168290757865634,0.0531423797623494,0.182373802418954,0.175093257668483,-0.345197307220990,0.377321935018161,-0.0371905095874881,-0.142045257234039,-0.201998681609939,0.598693982570734;-0.333041110855161,-0.742676794766719,-0.101425931889464,0.0410164638831049,0.788667157343342,-0.734554906792268,0.0902299138570059,0.450271017770163,1.23399801735229,0.699412683970354,0.702483715709275,-1.29691628316747;-0.133457536690729,0.0618159957170285,-0.139189406410811,-0.221742166817244,-0.0363691933669162,0.182265970740452,0.127960572943398,0.477158493726374,0.0862879091225723,0.105668200960081,0.359302066701225,-0.867510603453678;-0.0395732054575287,0.0181708382392010,-0.0361539930071259,0.0140515296718913,-0.00530903825428844,-0.222862225612349,-0.0604247954665886,-0.120840275825465,-0.0217247287611826,0.135790909619858,-0.178690718695544,0.175812258508906;-0.208655213844052,-0.471093000778408,0.341218345835100,0.301219725576611,0.163932721012135,-0.211876532818038,0.0471115300375105,-0.315541874303501,0.268317562906974,0.334345561587726,-0.190400857826755,0.283973100019477;-0.116906883011047,0.0819322796487274,-0.161156303784339,0.0691918891330354,-0.0234866911881134,-0.310449283638361,1.16965514477876,-0.692461167255106,0.298601068254984,0.983403611586341,2.03852521921848,0.0432658800008663;-0.206519464166826,0.00410898465177199,0.338214265399082,0.129685804641316,-0.375729366106118,0.0331556015684458,0.217363975626063,0.439145164014210,0.293812454728498,-0.00403821423644205,0.563433473035642,-0.243839809302676;0.853969434126323,-1.35724980006023,-0.373871655510047,-0.613107142096196,-0.361469667207611,0.256039042721426,0.224558133156466,-0.576860987272065,-0.104270399696937,-0.133455351442008,-1.29454811311709,1.30179744821107;-0.00421384736387175,-0.0476185452122701,0.579861574813094,0.149204543289521,0.705750163873917,-0.585302122957586,-0.550909843629953,-0.240437122279883,4.33499759878017e-05,0.289516765744160,-0.356211202833022,-0.398535681771771;-0.271902387947064,-0.197135041696345,0.246382939412045,0.110125985861497,-0.0571722634581751,0.0470379817939249,-0.943774462973462,0.519658101896090,-0.132498363120031,-0.655654467885345,-1.21259650596867,-0.866615695792963;0.0302603995412305,-0.137647605440241,-0.0209445051100441,0.145097201875796,-0.342037419008388,0.169580710072091,-0.102531482952643,-0.537269113968370,-0.901977363656519,-1.32040193673500,0.378517332475702,0.561546150163866;0.200615583194296,1.59283983673857,0.214689325368500,-0.0862911316146701,0.0533092302559819,0.365130378682894,-0.265283831134285,0.415712017525173,-0.0866076820218213,-0.235724200685392,-0.537873997342694,1.20789413778683;-0.0524400201511989,0.859653066501163,0.111621822040274,-0.353326624404981,-0.845969557467433,-2.04814713358625,0.372028513045331,-0.321289951885378,-0.0862089197446175,0.949631716942793,0.365385709331662,-0.0302429388909697;0.331962367015591,0.590809393395449,-0.559668552225226,-0.252531430241581,1.40687391299917,-0.488239474793527,-0.294774396160360,-0.262254599003152,0.0887178400850221,0.616109551017070,-0.433193094359275,1.35316793569963;0.120366228308716,0.116964136522116,-0.137990581227248,0.0314755057219260,-0.480225804173314,0.377704505099175,0.462298447600011,1.07173564927510,-0.0591682485297378,-0.546581471604433,0.934966319716984,-0.455439635897426;-0.221075530854677,0.206948432988141,0.162285846385649,0.0256816643070798,0.0948786612965728,-0.136136813354157,0.0750336947668184,-0.0561567772002065,-0.730696148327461,-0.201429531971550,-0.223699112875459,0.425898413070122;0.247781354403032,-0.266314467398668,0.239841333320131,0.115092997466568,-0.444234778625422,0.0601859944238419,0.204906786663345,-0.0407145776467918,-0.332170001422067,-0.151210279702167,-0.296060142056553,-0.488650348001049;0.0440658160305934,-0.0123171937147812,0.0353047276156276,-0.00204470404471473,-0.0626580426741278,0.161490988322888,0.0849338409237360,0.166976154806491,0.0195771778881870,-0.121084144455320,0.102513485666621,-0.227132715253738;0.916383623249440,0.609819758726622,-0.276506514111068,-0.0953967790592182,-0.308922997085808,0.317580089821968,0.935800715698584,0.659379435884294,-0.292856385493975,0.173929499933135,0.688574331799503,0.332852601067449;-0.611878230614175,-0.947848609700028,0.704016583415035,0.240768062207731,-1.40735954083080,0.526754164727453,0.165374950973580,0.239326479316813,0.152774375515309,-0.940246424696203,0.852201864641983,-1.54482465242214;0.0529526881945020,0.0136697249564016,0.0371877733873926,0.00428197776213371,-0.0959292266143472,0.139489852782799,0.0788746628534991,0.176551448635823,0.00346505145643496,-0.0987240615065693,-0.0137498756246835,-0.162084815227320;0.0382587801441302,0.0234800991984918,0.0253775845870732,0.00833754221288423,-0.00296819450709364,0.112887691766777,0.0129704037273697,0.0592616325990641,-0.00865939547155465,-0.0515531291890999,0.0333093791779271,-0.0887796698689118;-0.250435856651500,0.420133089811746,0.286703706671407,-0.0128019084402439,0.306581590674354,-0.0645890788742195,-1.19643624108624,0.0386571051960982,-0.228235531723227,0.130347883664303,0.0812705638846010,0.432125670882383;0.755154103646010,-0.120859836655667,-0.839655456366366,0.0407842507372847,0.920444565492506,-0.284025834357068,0.311185585685188,0.0217601840932721,0.516422433342951,-0.257306667636079,0.548133152513982,0.715864321223122;-0.196438594102249,-0.0627910875432695,-0.676417516631139,0.0110833897810323,-0.105811855932603,1.29865505628239,0.147797791882766,-0.238893858670580,0.127276308029088,-0.204559985597182,-0.213202564782683,0.289040646933013;0.0378995264997219,0.0182392309009913,0.0301947052268520,0.0111553950933635,-0.000738121946511195,0.118726943559561,0.0173863437115622,0.0615676727826874,-0.0103765422464631,-0.0514514553443203,0.0439082726137779,-0.0872723778780784;-0.00175940725317983,0.101997952127086,0.0159030986642427,0.0619138939929737,-0.0972019688366601,-0.149940903156842,-0.0759033770815607,-0.0103001056515503,0.0231637140111497,0.0692187473921881,-0.101020104292789,0.180866276233295;-0.0551709796968172,-0.346746072209368,0.145935477408508,-0.171770595884913,0.0712522291297466,0.449600069197196,0.0525495194279948,0.105281229896132,0.0340139481568730,-0.0788394578404228,0.580005101683808,-0.252834460021771;0.113768092607224,0.0170523072217647,0.155561746646424,0.0933719663606600,-0.100275355451132,0.101156651702738,0.0983486469489405,0.0645749753157818,-0.138666666651917,-0.118527079841720,0.00926899428025889,-0.0326569710218917;0.0358773843813773,0.0158646831397851,0.0179112652924643,0.00286139525181964,0.00452265215038973,0.123350780350701,0.0176146144323393,0.0596173319775326,-0.00848447192882660,-0.0599514196436158,0.0438678581612665,-0.106204775918588;0.0357436365142921,0.0527217671695268,0.0226462747548059,0.00937990233229298,-0.0153597291919714,0.0503606329187799,-0.0238838119786479,0.0335709993068204,0.00981885771756076,-0.0158261816068908,-0.00897851556237932,-0.0214822681527062;-0.120644779944819,-0.569901891386972,-0.668223804696825,-0.226890964563360,0.339266652509957,0.914606351071877,-0.560986637330568,0.781676172992125,-0.235250182438883,-1.47017518326377,0.152768257761226,-0.420753462016386;0.0374058518828675,0.0734516773627352,0.0833121970104316,0.0275123882517601,-0.0890155757777840,0.109189549111664,-0.00418566310351888,0.129843811949918,-0.00860317990707899,-0.0374452722137512,0.0327651368899946,-0.0725725776862501;0.0363745933588167,0.0122422134426439,0.0384636282215976,0.0193529030976409,-0.00484983432608472,0.130895333750579,0.0292408471295261,0.0748702769267285,-0.0130141981804666,-0.0494889401492691,0.0688149722306887,-0.0846527855295812;-0.177140258492542,-0.0437350671540060,-0.229065055054713,0.0326416784873502,0.370435426985597,-0.230400539716691,0.0891524141673784,-0.109597927727087,-0.0413352117434922,0.0623302042305417,-0.102739538044107,0.120149492827315;-0.614407116689405,-0.184613242835745,-0.330076557426141,0.274728883627022,0.487358003654451,-0.855670458510418,0.137001799316712,-0.602660942849993,0.293162826618539,0.323157831006550,0.675835164439791,0.297576256214441;-1.00928251999661,-0.118723053722703,0.377707877758600,0.0601637169128450,-0.519360133585773,-0.403710716020889,-0.674474027339821,0.335307579507400,-0.413001112826916,-0.511354454166388,-0.616644860250986,-0.325798001275247;-0.0424441744930118,-0.00919564839711624,-0.0449462613459415,-0.0128950448864794,0.0326942796389893,-0.140822721683967,-0.0422028549402007,-0.108798862971622,0.00631666215174246,0.0834623532649821,-0.0577179218717766,0.119538237024557;-0.0258979292085552,-0.0181604315365856,-0.0251160653614119,-0.0107600594364620,-0.00678752529993948,-0.0690475219387902,-0.00801490733374737,-0.0260679013532297,0.0123549206276953,0.0105211414627183,-0.0201007836378824,0.0383425966505142;0.464781509553903,1.20576777753719,0.322471718519324,-0.825446049912426,1.56401109503745,-0.968467153795881,-0.104677109485900,0.424024883270850,0.351385881068208,-1.41859798660444,-0.783702118796419,-0.0306850390197010;0.0237480090773859,0.862220012087148,0.899893857756264,-0.876284303781763,-1.40389440259321,-1.03212171395613,-0.235644653366082,1.07126672999873,-0.535231138891890,0.688930809977586,-0.820686209577205,1.63674119401584;0.337617295152756,-1.53865069628288,0.0202874526703154,-0.126670737187128,-0.159144017113996,0.393794262899855,0.702141798941201,-0.861626249806470,0.373831122310115,-0.0633506671266208,-0.171095204907234,0.175672327390834;-0.339144838696641,-0.0362428259036869,-0.0414396517938258,0.195062011535065,-0.756214577287738,-0.116081647237828,0.382209016662222,0.543033626507620,-0.0312385334671004,-0.0631025724281087,-0.154886157666763,0.416384963782385;-0.0376420519219509,-0.0252302512147897,-0.0316829914686249,-0.0126711942748611,0.00641064265967567,-0.113410276158630,-0.0140930714355840,-0.0624594319338620,0.00986596592024442,0.0474286434437957,-0.0361026106521359,0.0812473454526261;0.0890831594052184,0.0266478876255263,-0.301132080857063,-0.0538307380339512,-0.398046310167968,0.344414052441748,0.247594090131152,0.138353454150650,-0.173467326616749,-0.237789746906557,-0.344373526805916,-0.779680886290548;-0.205691558487762,0.0667422423898506,-0.271126161875356,-0.178537209380580,0.505635467922069,-0.0152877139233855,-0.146459121192162,0.0730958393600201,0.288793494183883,0.0793432196393333,0.393552552777003,0.358495653865534;-0.0605801653130906,0.883903791652781,-0.285580479233787,-0.0685954097722885,0.708036404828111,-0.436476619835432,-0.677144572158600,0.289160819069127,-0.316494161938225,-0.0412552554728675,-0.00864944034984066,-0.941331076438061;0.105990068864061,-0.504230097311724,0.0341473459545791,0.341725770656636,0.703900724075032,2.02373624384997,-0.346837591302776,0.00259640061290681,0.0747664897845861,-0.496240745499125,-0.763856976550699,-0.105937822381321;0.0850058116980980,0.738302642988793,0.894890161157336,-0.717757050883658,-1.18496352879049,-1.10828057323364,0.141112505768024,0.739458681010082,-0.255025170227437,0.587235933160700,0.0157189182933031,1.05548993199953;0.488134601681549,0.251714471258239,0.547617800941202,0.256488518987590,-0.0597355125396501,-0.237619171523464,0.165620664113919,0.197114287020284,0.276707674710804,0.129507330793451,-0.248214767110308,-0.544755060834846];

b1 = [0.0667747225480878;0.00292535855010206;-0.00790860464710153;-0.0328956039111549;-0.00930847105020394;-0.908672838172771;0.0132924445921913;0.00349796273688939;0.0112900116054364;0.000601445586478964;0.0197610783747647;0.0420498775486266;-0.301417675991635;-0.00408461624171350;-0.560265913804831;-0.0308957577331113;-0.0128761312549477;0.00241751266203109;-0.614232861230049;-0.0161571275918472;0.00429537393040562;-0.0255442138455647;0.00972758043714327;-0.00845453262465045;0.0300527076041611;0.00175848461411772;-0.0370981830393383;0.223662426856578;-0.0335647469952117;-0.0697517972375411;1.05080492953207;-1.21277585448166;-0.192678948806813;-0.0457678134150187;0.190566494107998;-0.339896160148029;0.232748024363620;-0.765867406787384;0.0179918982320627;0.355387167896057;-0.562936838927574;1.32010382965012;-0.573381893360795;1.09864511828929;-1.04393769629058;-0.106165971985803;0.0885659941880663;0.0408162613565113;0.540264711614004;-1.93160843345785;0.0664823499657186;0.0105184902605048;0.483181305437991;0.165748078331059;0.731224904235969;0.0117271158850824;0.0129153080243112;0.0992060984225980;0.0404884519961467;0.00990324210789132;0.00710053230383213;-1.51383077677795;0.0141460485422837;0.0160611594618026;-0.0448927513013909;-1.16788406852787;-0.540772595405116;-0.0316359533031746;0.00189021828629238;0.756315270558065;-0.00779563239713761;-0.691754756403250;-0.0443385937904401;-0.0115122002008384;0.178763009463690;-0.173089337997729;1.33140798065533;1.60521203951495;-0.183085931902114;0.374780434582744];

W2 = [-0.282369328207818,-0.0969748828710172,0.0130764386186138,0.203159064295952,0.125843335706285,0.391937414920437,-0.135093114990998,-0.0969185996101616,-0.138183605104639,-0.0762821300609409,0.249333917147191,-0.0511166321993168,-0.795388226122514,0.115720316179219,0.705197713521525,-0.661013863747642,0.168348384965467,-0.0894137766025967,-0.916905699651428,0.157185417545123,-0.00469258036612569,0.194208937494392,-0.0736878675854884,0.153721700339925,0.200554324928047,0.0323989147253361,0.231130181736840,0.821786098849415,0.474936792024647,-0.0495903320258512,-0.0571631673151938,0.461229394494019,-0.347298564484689,0.207072162154510,0.0632444957303062,0.755802211192096,0.235039870738271,1.02332582927034,-0.793328665813228,1.10375099666937,-0.261058258827386,-0.334569444429448,0.630702965217667,-0.609101553986495,0.493456688452778,-0.295621077710801,-0.512806120054685,-0.288744276670805,-0.370882493609227,-0.506836754213019,-0.288845174236169,-0.130411397570468,-0.740555450337073,0.629177288574987,-0.987887508129748,-0.134440401641866,-0.0271034496667869,-0.0437212407387268,-0.289250951782342,-0.125349975680669,-0.0988917107222905,0.289315094799513,-0.246869724441291,-0.155095006838320,0.459443605819438,-0.930942284331721,-0.653094835486670,0.211816596597222,0.0666998571956411,0.112087343776976,-0.582054015950671,-0.957251141666864,-0.0847613286397413,0.137881283112363,-0.658581964277278,0.383847593042844,-0.490724073772657,1.31609527027508,0.910514945365830,-0.272234732881620;-0.100244919537991,-0.0904690995383462,-0.0666508185827742,0.193088193987610,0.133858990012803,0.265655654620725,-0.186887719798197,-0.103124471210892,-0.0800204572912739,-0.0866185478477597,-0.148211447959816,0.0901935634505579,-1.65042382585486,0.0520668202231154,0.484059412601023,0.0281811264144743,0.0100569214680300,-0.0799454355112761,-0.829773550175322,0.131539144236178,0.0304551122838601,0.161928777117498,-0.670472186496133,-0.0421324450442463,-0.495195393640347,0.245825770878108,0.141432141371630,0.708470245715294,0.517971028321376,-0.891861161331656,0.874742392727201,0.836297975447823,-0.694674652474428,0.297127005132217,0.671407978378282,0.741398246432277,-0.472205976001724,0.941860428489809,-0.583587039974013,0.860392395401406,-0.703889927343701,-0.880174628840319,0.706540100586434,-0.794349001648483,0.532711141102383,0.394487416246562,-0.347012403725734,-0.216557616583506,-0.269324685311041,-0.678577213753539,-0.0868709742602998,-0.129062803544839,-0.441219589616042,0.586458760467958,-0.659726517146223,-0.138765784337412,0.284602001983808,-0.595725812066028,-0.175349904983471,-0.152595121867270,-0.0229418491296047,0.485333151205089,-0.0567052339751895,-0.151133882927367,-0.0894259265072299,-0.551145687937146,-0.565547426139004,0.156899776034272,0.0749709142677418,0.128302177235303,-0.686585435763046,-0.674737620712230,0.569458848112016,0.124987871113069,-0.232660633636303,0.441714206973140,-0.659330931569627,1.30895708684731,1.10600101223805,-0.656452103236267];

b2 = [0.0382934926562045;0.289253635064078];

return

%------------------------------------------------------------------------

function [ W1, b1, W2, b2, ti, to ]   = nn_18GHz


ti.gain = [2.16592611029879;0.00211822121778521;0.00585936801508902;5.57009709446998;0.00575559508682464;0.0603549007618636;0.0273979188601703;0.0356645790313325;0.0416423621568934;581.975479994452;5.36249120419891;3.13250047833172];

ti.xoffset = [7.62939453125000e-06;0.344390869140625;0;0.519999742507935;102.512039184570;0.0200001131743193;-41.5389709472656;217.145431518555;227.137084960938;-1.47336226552497e-10;0.684799492359161;0.373584717512131];

ti.ymin  = -1.0000000e+00;

to.gain = [5.277736148676918;3.101175065335577];

to.xoffset = [0.686223149299622;0.370642840862274];

to.ymin  = -1.0000000e+00;


W1 = [0.465202605808206,-0.179263871819386,-0.197089306220122,0.595620502098521,-0.0182938143464825,0.112990029034460,-0.204773460346492,-0.334619039649824,-0.624394791193648,-0.530578049593535,-0.173424472110254,-0.0674799763091879;-0.0455139851948763,0.112145826006256,0.138080925037269,0.117679462431914,0.00368893943440964,0.00580436921391576,0.119563810430580,0.00480273722331336,-0.0541916654554612,0.0256581973835312,-0.0247392211839892,0.0639808644486933;-0.0426018314816369,-0.0374000180386443,0.0147639380311336,0.0230489559273279,-0.000896026863280345,-0.0818785113705559,0.0773568428312642,0.0144620816567797,-0.0702091811924354,0.00740188375546681,-0.0474092089036403,0.0909590267895560;0.0554142675371651,-0.134179732965201,-0.157639855953447,-0.137500479191175,-0.00185932068940567,-0.00943374919320020,-0.134230105115011,-0.0124166081895612,0.0521135005308707,-0.0184861614541780,0.0198609115518945,-0.0517876493274318;0.0963610436341655,0.259398240542748,0.361685896187875,0.516513038113779,-0.562961816535681,0.109907885922098,0.166838254552554,-0.429603470747588,0.178104991559825,-0.499383648233831,-0.0306221711883112,0.553699338988777;0.176187545792497,0.247652040810703,0.892943572396414,-0.0646981136110020,-0.548202551804835,-0.192372294917996,0.875215482421972,-0.0146911579018460,0.619666226326665,0.719259663065058,0.228554913474152,1.25409025076463;-0.0935289351764745,-0.161990726020218,-0.0155094878760800,-0.0111317492944862,-0.0330602976636986,-0.200681167795641,-0.0628859194375009,-0.0393776988954808,-0.00635965419940839,0.0138199020561336,-0.0616374509369702,0.109750416864759;0.306651862327877,-0.367805409426562,-0.273098672640674,0.0933741501890680,-0.167227165243901,-0.343536330628012,-0.479279950157847,-0.314712240672585,0.0563283826756613,-0.149169257889686,-0.0281706001995424,0.170272867564851;-0.0486219811182846,-0.243371966496981,-0.0972203202165165,0.0411686814113557,-0.0244427954277366,-0.0713502297485836,-0.0979461958068935,-0.0331871530744189,-0.0193198386638558,-0.0336620456972827,-0.0959620273794193,0.0794379069220632;0.0130999386801546,0.184172531175012,0.119156124807202,0.0529120181512209,0.0801356631976392,0.0360988775400016,0.0555073049491477,-0.0114836202574653,-0.00193281157543540,0.0460818311409772,-0.000430395823654677,-0.0545742253029416;0.325447095407491,-1.31335509393208,0.331476674657777,0.585493625292187,1.04694115348084,1.22364280209049,-0.717079220774767,0.0291353053534231,-0.0987911652091347,-1.30729301244171,-0.607516436606070,-0.278580983306686;-0.0337065413327431,-0.0398286855038727,-0.0855827068625482,-0.0856539556904805,-0.00596477749156250,-0.0447881705288045,-0.00595583512179450,-0.0204967483705430,-0.0725842462860924,-0.0210921461424147,0.0244316438516713,0.0375005652789644;0.711245854310267,0.322060087143275,0.0793035917677531,-0.114340532620893,-0.138072091474540,-1.07246370369641,-0.0443792431315297,-0.133582419682042,0.00523862117981740,0.467342736710847,-0.130659487371695,0.479127742819198;0.426410863365328,-0.439465685280927,-0.265955322233069,0.144828995652347,1.16051139396175,-0.206806585197223,-0.157598788000067,-0.0975809705881190,1.07423415369583,-0.378714878722773,0.560444286806877,0.725429439606678;0.220349035521814,-0.547935962707836,-1.42277348723403,0.151390779241517,-0.273885405752503,-0.193637497492167,0.211407920775810,0.0433592971390566,0.131340089419863,0.154108676678515,0.355009513828215,-1.05136756301311;0.0103611496678992,-0.109677876409879,-0.138790032800115,-0.101939300821674,0.00473023689468708,-0.0232885886641114,-0.108954395931930,-0.00864066000939526,0.0166518702628322,-0.0106059556397030,0.0312430799661104,-0.00309552799061706;-0.00854267922150557,0.467541900736269,0.136672116982854,-0.148323299180395,-0.00856330076687902,0.786823911328141,-0.127995770457281,0.167878390417487,-0.124675789052726,0.198511005448895,0.402474904419500,0.393746850195638;0.621902546756434,-0.585966533628808,-0.769798957106831,0.147617717668769,-0.697572286074965,-0.0117190859070877,0.413577947127505,0.466882184364876,0.143749342991474,0.622920856451211,-0.254978401748009,-0.118062583807829;0.283539066259381,0.192302099234163,-1.02852385036803,-0.151881756528019,0.427654016980595,-0.0716603823722538,-0.335192480071363,-0.268820656635138,-0.397698151409252,0.0310961950173315,-0.876073128031931,0.563468206286978;0.00532670465555091,-0.0126388582326491,-0.00246874805492838,-0.0112710969807518,0.0253672559017686,0.0688887157932519,0.0528174238791431,0.0746747212788990,0.0507032633517989,-0.0640974092199415,0.0185044915026385,-0.0123063249654706;0.323373807036555,-0.305842547478531,-0.278618684513264,-0.488085620182819,0.710548424777512,-0.356134300417916,0.0417762626034435,0.390408615059316,-0.314688314891934,0.522781058905329,-0.781039965372746,0.977226889494696;0.0304442831303498,-0.304007369263399,-0.0424067634855949,0.587467979465454,-0.563869762240588,0.0809575139828885,-0.322868390134706,-0.299127664009366,-0.383144250662886,-1.19364773356390,0.467789532969256,-1.26605583734739;-0.0682959646172107,-0.0942406849859657,-0.0818482816349119,-0.00646575573201756,-0.00707842686757634,-0.0573346512709008,-0.00158880967843965,0.00435184317341736,-0.0523033046326309,-0.0144094205890611,-0.0174843712664423,0.0706311261242733;-0.398325079572653,-0.255504327839603,-0.450613364118142,-0.267401924123777,0.0271570160202553,0.435399176860035,-0.0346424147905746,0.342098090184364,0.0485642755089703,-0.0379582398568006,0.279134960119190,-0.119707613669574;-0.0501902680329274,-0.105584443293766,-0.106153609832047,0.0328141517356888,0.0181708908516841,-0.0991941984390903,-0.00473571204651500,0.00371102668198270,-0.0566062894807219,-0.000957988337925057,-0.0432774545594705,0.0836056952427558;0.0561621658698312,0.0784960301860298,0.0676507297842996,-0.00251714357399085,0.00410726815202155,0.114898963400329,-0.00144257524884526,-0.00245064310646495,0.0503050350249692,-0.0248604175686125,0.00356975564513463,-0.0928597178997123;-0.0514477883193518,0.215181965644438,0.151735670404981,-0.00266836506275130,0.0884105133247145,0.129379904568420,0.204639754353921,0.0641979900948334,0.00963214379078208,0.0151060235565421,0.0543484867353320,-0.133939601715746;-0.212436497661644,-0.458195095241918,-0.607288574035415,-0.0183100417467952,-0.143856129226353,0.901843872324078,-0.0855585998685389,-0.195177625377632,-0.327260675893099,0.777276365205674,0.682948518518383,0.896207604702400;0.547437157391456,0.270096262878474,0.314457533050648,-0.224035465725296,0.545675545899119,1.16100236704478,0.293422373652973,0.235651069017007,0.367111388128583,-0.407174987203051,-0.136420222493097,0.347328028091769;0.0403272638189849,0.0390573083204059,0.0989843996368538,-0.00753703898750485,0.0153181153849440,0.0512395059178609,0.0340307365963107,-0.0340810317352779,-0.0345080425640619,0.0158688294687727,-0.0448875299351361,-0.00824392441228596;0.237183519782034,-0.286608407702768,0.339683431786307,-0.237480555267737,0.0308932160190170,-0.0633599955072989,0.0256608422588638,0.0802602379886501,0.226193180911172,0.196113643839428,-0.264212867627111,0.118308084616839;0.0286422801607381,0.0277037923736068,0.0113548765183477,-0.0566409160099397,-0.00804073282521865,0.0662672662546272,0.0230095547277921,0.0272176427025960,0.00658093334612420,-0.0293685104496076,0.0185492395302116,-0.0682206252293410;-0.315529365878757,-0.833572524364328,-0.476704701734351,-0.282767524032922,-0.965500391603884,-1.02033664700161,0.141678617333879,-0.633188592572731,-0.184404623136007,0.357249994364555,0.919843797235542,-0.347101979850221;-0.0254029803540068,-0.195622679356677,-0.180345214513747,0.0316623830789579,-0.0653578780346511,-0.221445355314637,-0.0685014510654991,-0.144313019966350,-0.0527374647544808,-0.0693111575282059,0.0552656459484668,0.0173200652413695;0.0668245980885408,-0.495463804121462,0.900496600199289,-0.143165850250824,0.167630690119415,0.622540667441362,0.348634559882505,-0.402114781547972,0.299933280908677,-0.0449125946153908,-0.232373182158765,-0.182865959140787;-0.890953189970353,0.725919517976560,0.170870009432839,0.890270640588019,-0.954079201238932,-0.155694794130714,0.426256385462669,0.526471189089271,-0.0582412948481786,0.877680003235645,0.491931160237414,0.197915506612241;-0.0980744505824677,-1.33040469785921,0.263147923277867,0.286793256989223,-1.72181230161899,0.555843013685347,-0.0973207725374095,-0.0974982040810272,-0.260666276042930,-0.127585877966031,0.661327108863687,-1.04333567991916;0.125554985645869,-0.407416379832366,0.114351211270531,-0.363737997493519,0.528246636515170,-0.288289421291705,0.148125338524063,-0.0653495361222172,-0.205317965509507,-0.165535881692467,0.148356145994998,-0.390111727307110;0.375964076869784,-0.0514649517008095,0.245322308598380,-0.251875895400319,0.604080601584547,0.448111146458342,0.328451315303392,1.15523980798736,0.380344857981057,0.887333863617267,-0.0351298311057359,0.735332991364392;-0.263622775320036,1.07811033654753,0.104324451526228,-0.175482638621340,1.76028963023690,-0.437710736466575,-0.143253334370154,0.783646017606164,-0.345601162223914,-0.0608094075976192,0.118795212286410,-0.374095496721952;-0.647458240620276,0.618957098663295,0.277127579339918,-0.0924624999347999,-0.656079268812115,-0.209292641773329,-0.880667862529906,-1.05735346940467,-0.153997649205690,-0.517237120872473,-0.475354693593969,-0.542660390455256;-0.769058008344274,0.229312704181536,0.0467820330392563,-0.0778199642135289,-0.592320963776889,0.0371685523915500,-0.836916294089288,-0.187569761045381,-0.0601291086952947,-0.120375028925883,0.249776705926443,-0.533603165680755;-0.352676838311660,0.630330240934003,-1.08015085460176,-0.190878008910780,0.253383643236106,-0.292465150702570,-0.192986613513695,-0.786326852485944,0.240954906206698,0.0644845118367556,-0.0629822121353433,0.199087524433840;-0.0109104892446625,0.0253708317953385,0.358883456887448,-1.01999690704033,0.272753571998427,0.104981824800304,0.0407920967128296,0.931714456995452,-0.0258898420380617,0.288529680311582,-1.37922913496722,1.71800258549574;0.312953014757435,0.450018065157653,0.994829845827382,0.126344815438642,0.225264826822837,-0.0629581535405485,0.632283976473009,-0.421414530900106,0.139329586751447,0.481890268538407,1.26266306610823,0.983541989450983;-0.131683031205477,1.14770997246753,1.08441415756041,-0.0757606545571717,0.148979037091733,-0.165520474195869,-0.368843210657501,0.324990692509174,0.0527288663252997,-0.103304184895661,-1.54177160281923,1.77435321732306;-0.381858215251915,1.34432604422707,-0.358184065684147,-0.0230700416882932,0.112244277480940,-0.569736650946875,-0.765684213013066,1.39241166365186,-0.866327370673117,-0.406266912370920,0.423512117753603,-0.210187424669040;0.0894231820071959,-0.727773553617130,0.635104171027452,0.308461137701656,0.273756062857795,-0.676032896674799,0.161073992392274,0.257055804943479,-0.259234350004919,0.535047883458416,0.719209136836766,-0.952310571601603;0.364527310525629,0.829778767331043,0.578473871718836,-0.657563218478461,-0.544703511247456,-0.290078318569370,0.367280985013611,0.742346071902827,-0.483971733288852,-0.0835326734712894,0.347636750798229,0.408187181907342;-0.0690269010973622,-0.591381268033946,0.268080835127369,-0.0209987770030889,0.00710599438645432,-0.989395005086737,-0.159778734072010,-0.182979408284638,-0.351900360003695,0.775107487836257,-0.00351125811442375,0.279418888629041;0.103269664031403,-0.398417191476857,-0.600799070166300,-0.209691875303503,-0.288432909247495,0.161153567068222,-0.458755209093657,0.419337808843577,-0.156916975228272,-0.128434170173400,-2.00228670780014,0.230361232027016;-0.0493010283717610,0.358442357516756,-0.717832941814744,-0.125853678628039,-0.0445299320277020,0.121898435540695,-0.455989999191673,0.553598452065606,-0.677498094556304,-0.303096045965794,0.171780636446475,0.0850524610353134;0.139484774337293,-0.793705067109051,0.949866185607618,0.166561042621031,0.0653684600962576,0.0824831534616277,0.138807084208572,0.740221213567871,-0.158588561525790,-0.463921770041682,0.454449110406883,-1.29724138240159;0.0893582087065249,-0.0738540921184329,-0.409927623418046,0.267749776672487,-0.00569925992193736,-0.306577590519804,-0.295410051972662,-0.723988907968705,0.200118680483988,0.719194042080710,0.763693964298216,0.979170744822640;-0.236447387673852,-0.552689203427870,0.0327383994590351,0.345009131179963,1.14043758335043,0.0792400857460859,-0.00275978644085117,0.126014088335629,0.908916757489616,-0.587311571298327,1.05300709457244,-1.24349106732828;-0.153476047869587,0.115292671141092,0.482196284526557,0.398512472457637,-0.0324035472255427,0.184116785824088,-0.165840865706093,0.322126128842799,-0.0411284774296671,0.238511802890959,0.417878560820245,-1.03291767877214;0.682547918616834,0.0226384885174102,0.102083479514221,-0.223943751923040,1.06903378126440,-0.273814507057440,0.581500392967353,0.357696732350144,0.0169467615068280,-0.840993360066742,0.146458264494684,-0.205589965471523;0.238335780511535,0.246626866068822,0.0203228442487750,0.150105284814650,0.0884939183445870,0.205424635090366,0.331666265188948,0.222901059080487,-0.136660837055528,-0.0261403452845183,-0.331660815946963,-0.0494147381768875;0.721324305092198,-0.293907956462299,0.472509008362135,0.336117883114951,-1.16223542079125,0.199374255014752,-0.261067068501696,0.582588842867885,-0.666549709238367,0.311605299986559,-0.600652306818813,0.701715572397175;-0.546486531323792,0.831543612832040,-0.172275114798235,-0.577070451903082,0.257823396310217,-0.159047583900494,0.651782029391681,0.949490420368163,0.174303522143111,0.662075468432468,0.431410820746169,-1.03718801775816;0.254679347132889,0.119733840028315,0.0114882136924224,0.0972371898704375,0.0522959981841669,0.0224211988623504,0.0182103743924304,0.0156472612010684,0.0189188786436302,-0.0106290037169505,-0.0880224795957860,-0.0877838409742175;-0.169453590977352,-0.781000776366769,-0.0412436806193922,0.451360645322829,0.172190630022330,0.526351236924936,-0.0639224363576435,-0.319454625743156,0.838037833289459,-0.232077161794195,-0.278145479297559,-0.608490485440503;0.292527666616788,-0.839593054945591,-0.103076549048780,-0.0751545307303662,1.00099073754651,0.234239552403859,0.515773031274343,1.69469644641681,0.134968245146062,-0.0314864367289242,1.24329899777722,-0.397880617391956;0.0983125300263074,0.0369079642712412,-0.226913164919147,0.118472687093777,0.682812268223774,-0.0325988801696675,-0.462606820746278,-0.225279247602372,-0.0378341784222389,-0.0821874268127738,-1.14973253593840,0.504897280174695;-0.0912346895974599,0.154385876178890,-0.958303495798384,-0.0557122975964277,0.371572145754361,0.376751094545042,-0.704009143640880,-0.306049011991297,-0.469130398922034,-0.637813335097500,-0.413675119992976,-1.44734275681774;-0.0743652621440170,-0.416878567230928,-0.674067051812995,0.877339485682578,-0.102440050150723,0.0519325047160945,-0.265084895127156,-1.02755857810022,0.0923096813840120,0.458287216668121,0.545388797352059,-0.178442194391923;-0.0228341887318494,0.270479552350045,0.349175819475375,0.485008843883772,1.08746763226534,0.755290622735413,-0.211359132430738,0.608889728793047,0.413173021011699,-0.199447455081899,-0.532069844127851,0.0511916886866554;-0.505550942367232,0.226715881425731,0.0398382186425837,0.522492595706367,-0.229566105084474,0.0791996615260349,0.282248781679353,0.346454054183018,0.0927274675642198,0.189637172065165,0.283998005851348,-0.0783278237997171;0.144157767835294,0.467035765755164,0.314060977157233,-0.448914324921143,-0.564183213161031,-0.0698759934327052,0.427032876804403,0.371566103438817,-0.161103276468034,-0.00665706695095163,0.461986362144658,0.130614279227750;0.269046015531136,0.159242199374059,-0.109348811157147,-0.436909980692242,-1.16948503701797,-0.320128442674402,0.161855483962195,-0.570848546335103,-0.484666633557925,0.266678439148850,-0.540573336091805,0.796052002775206;-0.112550803393510,-0.753665507748205,1.10101829799708,0.336220554543688,0.444692840531858,-2.00547752319472,0.000161278562251476,-0.423679407407071,0.772267570671399,-0.442704364035982,-0.321341051487004,-1.41483376159259;0.331291803301185,-0.548660033143883,0.0688589814161463,-0.201862263737188,0.799039859052685,0.322307692363722,0.284854627548794,1.62328406959403,0.288691813540499,0.207575447623793,0.478439381301290,0.0312608163674548;-0.0157323874366356,-0.162353132122088,-0.527537809488518,0.219479921183034,-0.267155629834367,0.195823778792473,-0.184230324589836,-0.0596236170444760,0.0116896430044051,0.0372522590178817,0.0206089636084113,-0.491274293081001;0.291266708443243,-0.851240533174866,0.326693884211148,0.597885728062947,-0.193560163928612,0.391080403966356,-0.532822915860175,-0.835168595921943,-0.200163687133768,-1.39905101484676,0.0463710588567872,-0.490179007019080;0.00614150012577001,0.955181287872197,-0.608786607593072,-0.335755689474491,1.75333045949702,-0.220936355411965,0.370490377977327,-0.630195084938802,0.905132446877496,0.106181335061695,-0.202906498351284,1.21306063041661;0.0505041402132470,0.0961167039967077,0.0899703708368779,0.000879652691165626,-0.00314824574054738,0.0450630354887141,0.0299030259708564,-0.00568361783408164,0.0129116698064695,-0.00137576966079699,-0.0300502457317555,-0.0671017637655246;0.242525087983944,-0.659399604323763,-0.195600053138687,0.0130607385285841,-0.636805093623926,-2.46199407930991,-0.0575656139153990,-0.632966967382260,0.234933681423760,0.692554016987599,-0.427155967668410,1.36624872351364;0.0165195336805786,-0.0390699405042022,-0.0298822409526092,-0.0634174945378281,0.00322085613986451,-0.00449089302979464,-0.0617198064222734,0.0243090632928400,0.0570712659826827,-0.0188442555126156,0.0562009290034885,-0.0451113074801060;-0.295570175272150,0.569084619754760,1.73731241843823,0.0924339789612331,0.515515181082687,-0.262568896932586,-0.507457981362285,0.281994502524292,0.183248899613978,-0.266037349142528,-1.39649717994080,1.58129207784731;0.0197699073264272,-0.0927288413805367,0.673387632088191,-0.112525314153534,0.585643707855163,-0.378361445282227,0.130389534019529,0.0555749783022669,0.0785719946198918,-0.0423055611869791,0.181279604707086,0.434766608543037];

b1 = [0.338455996581804;0.00477786875095028;-0.00726643338490315;-0.0107897094391708;-0.827929253327049;0.530623077428514;-0.0240731920322564;-0.0827240828817696;0.00726744041382703;-0.0168098453594796;-0.324296453441883;0.00462557139911802;-0.435711976758814;-0.126351720819354;-0.970817310446271;-0.0240856835010129;0.395362553966083;-0.626284610598838;0.203028155961825;0.0480403797480533;0.414738469194535;-0.120828497692361;-0.00553608364192215;-0.328760885129803;0.00730596282769042;0.0203956434191444;0.0294378497391451;0.566084072645659;0.742314088424336;-0.00390248144342434;0.198362645455790;0.0211748988243718;-1.32147592054948;0.0829501029534627;0.415994871022211;-0.552525193595931;-1.28140137378220;-0.209515112996828;0.0776962314634874;1.53349872418728;0.610683563776998;-0.0247696749072700;-0.356688080393094;0.250363175720723;0.688779481855009;1.21741381368424;-0.318084341431422;-0.465620410455977;0.231431039283255;-0.942075462135318;-0.0776076881431494;-0.352294169294892;0.221880397284446;-0.194466097378601;-0.988827393692191;0.144962774089325;-0.362882693985438;-0.0358985502675963;-0.0110051176170671;0.299447240817264;-0.0308353293095378;-0.135402938663166;-1.39539375769841;0.250035887646977;-0.186309642055870;-0.437333143176702;0.538562388200743;-0.348917173308412;0.133766885310685;0.381921699829420;-1.20869098580868;-0.930423883666150;0.0827998814084606;0.0186950167403460;0.817189909191036;0.00365186854182985;-2.62297365721710;0.00122332686889238;1.21332010843247;-0.123604784842766];

W2 = [-0.547547760555094,-0.224105534578173,0.00151965204568449,0.253980312014135,0.380472772164987,0.946205406702021,0.192064937308713,0.416202574949445,0.237407227224790,-0.282106906348502,0.293950569651711,0.148407090608536,0.840499205885218,0.429262438849867,-0.822943958130580,0.243540736260035,-0.799938654994642,0.691299735541512,0.784684399319803,-0.00799403968352118,0.790388502961761,0.891902785871678,0.148994971015256,-0.363275264138991,0.168333498428884,-0.142361114330060,-0.376466124845540,0.713018742620020,0.929520630589146,-0.151158542526435,-0.717036294656533,-0.0255380842435430,0.841476917531910,0.00670529401752870,1.13712873915648,-0.201570407471922,-1.07954114886837,-0.194670025126586,1.13615123223106,-0.777113287952257,0.672939721765174,-0.600147971821934,-0.855569362653102,-0.608824632048900,-0.802519547466605,-1.08469391437030,0.707660235206066,-0.770133474767114,-0.287411293906548,0.908913448693872,-0.889037235004729,0.504274899261512,-1.18148082548648,-0.502634418105824,0.896506896493986,0.400621722974284,-0.578408758314844,-0.300974439120754,-0.630914692346226,-0.665035650503525,-0.174057667681127,-0.745882073033419,0.968672139642561,-1.06460088470133,0.865291276113438,-0.726630473582601,1.22602399889695,0.911457187944169,-0.326891409372954,1.52406508962378,0.720950411373636,-1.31515078070277,0.704916189416606,-0.829524723194281,-0.826494021860991,-0.166199976997804,-0.855808755826071,0.0949256868393982,0.684237191576273,0.624883987513719;-0.176495240356910,0.129752342135556,0.172197401669978,-0.130404392893640,0.757519679533219,0.768434498259441,0.267487825145609,0.589640517860331,0.210090210459496,-0.0775576941100431,0.271211950520339,0.0303989344190440,0.777963017564332,0.364070146845904,-0.637727609413150,-0.0559876854768763,-0.824498260650944,0.484970991590411,0.480788140838190,-0.0519071711250481,0.883362694911797,0.601838081058114,0.131107848893524,0.414906441052288,0.200866312067092,-0.164125267396683,-0.253000264535644,0.653983000122814,0.629971361705614,-0.0217775827400821,-0.0452784613385834,-0.104096733543971,0.807973258692115,0.283568984580145,0.866534024754740,-0.254933769432816,-0.841791059090136,-0.586797112049347,1.00293210096396,-0.613861192835930,0.742077250675002,-0.906343067092468,-0.831798043309469,-0.415007949916111,-0.743593744538482,-0.949639520824971,0.457916401231799,-0.840104993245449,-0.606945309114934,0.846940477542489,-0.856748517415267,0.776339212273643,-1.10825972694341,-0.544490221840804,0.890182389667991,0.843471305772031,-0.544159700091810,-0.556513992413051,-0.603803105296141,-0.543765018310632,-0.274569596376252,-0.583950078420257,0.987812071329343,-0.0627391330945571,0.675748310569180,-0.557749791499399,1.22165521296726,0.828390972861255,0.710407250952004,1.52941192776110,0.650750169486695,-1.24204079271816,-0.278615062466138,-0.700692049778187,-0.674395951693083,-0.0958290796146878,-0.730294667527731,-0.102993422726046,0.621432911312111,-0.0192880569217222];

b2 =  [0.114247195128315;0.212090139114435];

return

%------------------------------------------------------------------------

function [ W1, b1, W2, b2, ti, to ]   = nn_36GHz


ti.gain = [2.16592611029879;0.00211822121778521;0.00585936801508902;5.57009709446998;0.00575559508682464;0.0603549007618636;0.0273979188601703;0.0356645790313325;0.0416423621568934;581.975479994452;5.19291201552313;3.20018311594673];

ti.xoffset = [7.62939453125000e-06;0.344390869140625;0;0.519999742507935;102.512039184570;0.0200001131743193;-41.5389709472656;217.145431518555;227.137084960938;-1.47336226552497e-10;0.643690526485443;0.383905887603760];

ti.ymin  = -1.0000000e+00;

to.gain = [4.761602513907577;3.137858153692472];

to.xoffset = [0.617577075958252;0.374222934246063];

to.ymin  = -1.0000000e+00;


W1 = [0.0269471887128286,0.0234334966948785,0.406858715731415,-0.315235454922834,-0.0567707719890437,0.242909549612316,-0.141068072816123,-0.217186103034490,0.254525027551980,0.123537893413463,0.0762244095401828,-0.0175795946421705;-0.530903419703216,-1.33383748726005,-0.141166855833623,0.489789647258411,0.0841188942673084,0.373815256307419,-0.408578047928374,-0.737212342216969,-0.0172431937768387,0.599056517208583,-0.00502452932020081,-1.13564540515724;0.0935284023496395,-0.345338220848921,0.391838144963310,-0.514872703232854,-0.0364158327491969,-0.723684813000490,0.594537521125057,-0.349815163216613,-0.192851037938194,0.0301595057268434,-0.376751156586569,0.290656007163709;0.550366708030547,-0.274757543849979,-0.405758817943743,0.852152134496257,0.213969305187058,0.329637386736728,0.0236541787300287,0.0603867682214910,0.119142589240120,-0.304248724964512,-0.121488141604313,0.543167516611504;-0.237148752409827,0.183740611849855,-0.254425890612870,0.0456500555044321,-1.00440159589061,0.381779310019747,0.965734924810947,-0.231586200709058,0.544565710516910,0.0475961547618404,-0.212006000890299,-0.592086132980785;0.869487888772860,-0.714407965038448,0.276445628168113,-0.0606508350640389,-0.646073040716342,-0.222448502616922,0.464205420495605,0.212072268572236,-0.285571495902958,0.193073773903677,-0.442066145701020,1.27159409343232;-0.241503475273822,0.0889015926607241,-0.202086323780927,0.327882300387076,0.979574765653329,-0.549923265477027,0.683197859392225,-0.502545985443144,0.184376569835362,-0.359940832761508,-0.216379098600797,-0.0436675552029582;-0.254420931268458,0.260928374359432,-0.0729332258548130,0.332598558928023,-0.714933876185779,0.625548329499720,0.861331523086130,0.225652862994534,-0.0491095416814242,0.224236117473984,-0.669349114106491,-0.471372339620071;-0.821858369261960,0.894585719556599,-1.22346234234833,-0.203979200524262,0.708595106119951,-0.0246420681349270,-0.425155668881959,-0.841639119113334,-0.229334080676240,-0.680670668354146,-1.02288756004923,0.402481320905506;-0.578230083986990,-2.32956240284671,-1.22171642981163,0.976683680085012,0.0144772945004425,1.51917651242914,-0.983986301313550,0.187699480926629,-2.22538514152695,0.550985956676262,1.13287509113328,-1.37800704555302;0.122356278168510,-2.01293970939124,-0.393183212371771,0.182374717668263,-0.0796723440940494,0.830187158815658,-1.28507319113227,0.413115091488774,-1.11779977490861,0.161315290318342,-1.28323394674952,-0.113607111786556;0.192489686083449,1.34971777359592,0.622182571858639,-0.0664525092255548,0.526449162862482,-0.256144405899100,-0.298815944295378,0.809688273922054,-2.19122980290566,-0.222992874928226,-0.850709452378550,0.476492629200931;-0.364998161385672,0.238676897106456,0.196824589963368,0.414599312610824,1.61863016163922,-0.202819339963201,0.138792833899479,0.358819698251163,-0.915577730434161,-0.767807600705017,-0.0782178288869375,-0.186775761157952;-0.375722311560921,0.915216853292653,-0.582906788410506,-0.0163012497119145,0.0919854769223111,-0.152394338317089,0.745591738849947,0.372269916269443,-0.150539270057569,1.72457739019768,1.54340023811823,-1.20873253489426;-0.272967531279778,1.53095101032387,-0.585248065495634,0.0918047694022640,-1.21947622141354,2.53638777227989,-0.580945342003043,0.0151337524895724,0.240434912615723,-0.296102179173150,0.286379357572813,0.513216215313273;0.122916126327226,-0.140335516059382,0.791817630971483,0.226680561006895,0.230242001635334,-0.986353638415208,0.0965813735190812,-0.313793076362552,-0.146490684863506,-0.581424740416369,-0.254740966350171,-0.287620663806391;0.0229497435831104,0.115117158894693,-0.412537150910078,0.0195112694871798,0.277212481429062,0.251101735681782,1.09529740969338,0.626525906207961,0.730426160118485,0.325609336604557,-0.729122637581984,1.34230976942690;-0.973927767999274,0.620741767689356,-1.29386602441716,0.00302865634077397,1.16613761419031,-0.246415390281862,-0.0127268586289189,-0.633644048922166,0.0432013006022991,-1.59269404031359,0.311160028088214,-1.47372251067849;-0.129640858735159,0.182469704575523,1.00935001454633,0.541746990139091,-0.165843682687374,0.253751220218682,0.600633474099309,0.646629229256249,0.480779352032534,0.0681834124568309,-0.884746441727761,1.08906512441626;0.903356189525155,-1.68113888477756,1.22110306317723,0.463944334236571,-0.295261449065291,0.506154995387528,-1.02189020915579,0.935035607412693,-1.25821594540690,-2.07488447882801,-1.74490962367836,1.66908329294937;1.32455190826790,2.54223359973639,0.461518141635173,0.0357948836247861,-1.44233464035007,0.0505354569615617,-0.443983117601750,0.298526282850650,-0.255194941199319,-0.343840148797425,-1.81937184949118,2.45744681360185;0.180329405904700,0.428839010110341,-0.173683719713418,0.179327708280171,-1.08084012392718,1.39961524624371,-0.145705465310562,0.0517548442144022,-0.407428270865902,-0.869025294563271,0.361361470172870,-0.540848265378996;0.580118130827540,0.918863877426343,0.677665902582831,0.341713567898363,-0.459900506741443,-1.34622856137577,-0.360097502121139,0.262242927764230,-0.401040354343464,0.112458211833513,0.143680298121604,0.568660288196177;-0.908502903428374,-1.04352474840079,-0.703285598626394,-0.277430988099692,-2.02753789390334,1.35271250611379,0.286685109295330,0.426471494229827,-0.602291411921650,0.0260822585912721,0.868887047321097,-0.575755863714721;0.312637628572895,0.917878836402150,-0.652036146959566,-0.0345446890243715,-2.22650028912544,-0.188839949493434,-0.710319024235236,0.385369570650190,-1.53521017876858,0.943150320519345,-1.16205103628286,1.87196920479198;0.0409070400282568,0.422776632125176,-0.104581081950649,-0.154227772167120,-0.307194193295581,-0.00516925356544496,0.493570811125102,-0.631210851607772,0.239795849133475,0.552706211798516,-0.0634078098628573,1.13284013743503;-0.0176774962362281,0.0593575450649943,0.826737363825631,0.727085761777258,-0.0980523492767233,0.643522420763629,0.697375285174605,1.44360607687329,-0.922779457739177,-1.27725569604938,0.187636159846627,-1.32758455175926;-0.298651271562267,0.255483228907024,-2.63688715804845,0.0540651221665590,0.597174758919836,2.78999468242972,0.203406655536718,0.0154134871999850,-0.418453315949400,-0.606614194062588,1.22742187023768,-1.38856942660087;-0.279266764427777,-0.883054797536915,-0.640887919956985,0.118776322650983,0.173618175577640,0.177412006328768,-0.355315534551136,0.140295127751307,-0.00245282423414949,0.650808249774214,0.314344576283900,0.361548667711165;-0.443482462545475,-0.117473073371291,0.555996078459413,-0.147831515752856,0.313568374008992,-0.332871731689266,0.908605479110307,-0.132066438555526,0.797674406478050,-0.910444972107031,-0.0423023483769479,-0.780971089849902;0.453233091044067,0.0760072651528431,-0.208600650216782,-0.180713958528504,0.538473168921281,-0.216231594072711,-0.387763341921325,0.0357049957458359,-0.313938771986163,0.00908577607289853,-0.430014460621225,1.13595355605073;-0.144377420383476,-0.568774286478683,0.109676554874307,-0.163831533788354,1.59407613691223,0.440545375204449,0.568134908855844,-0.0575826905076709,1.17230693786148,-0.0559141087476358,0.0769392799647779,-0.165848227620919;-0.146786805510923,-0.638516854754136,0.336055976673460,0.557916467816039,0.133634943543057,-0.0271727782356885,0.273324335644049,1.60913535248169,-0.226380179431817,-0.571129342470428,0.0679594170613800,-0.934032910989107;-0.264635414189943,0.381075286821416,-0.00374858458043388,-0.0438802250957884,-0.234130220408106,0.0762954765703724,-0.243993022573980,-0.189747995173463,-0.179481977448195,0.256215154726479,0.375845317242866,0.427665763563920;-0.0811890929030859,0.725477272217607,0.612062494985441,0.0589851519027840,0.411919148888744,0.398841032334047,-1.29022475636208,-0.759183417185307,-0.984192633659176,-0.762692287499924,0.393926592750682,-0.536593063200827;-0.0355188866083718,0.416835285008771,-0.473504701661755,0.314553774658484,0.433254885034162,-1.06755769538399,-0.339317211525577,0.951556606534951,-1.26598635040732,-0.369109532757480,0.129924620084740,-0.831397825717817;-0.995050304086656,-0.136365817393354,-0.542600991018899,0.375488304004600,-1.01952218333546,-0.517918526287364,0.444092211926763,-0.270020186473628,0.127296221917308,-0.0801022470970406,0.281845248107191,-1.20583929176437;-0.258175281610069,1.72150167036105,-0.601404347606784,-0.243149089714219,0.429344380019416,3.28682207781334,-0.796512271871869,0.564018867244739,-0.110404375824701,-0.461386926752972,-0.708451480773838,0.463398005983016;0.496271805328255,-0.762778685486125,0.821277218560931,0.160544975331686,-0.0619358981154643,0.0755448342608462,-0.639377596828229,-0.475895119761709,-0.166561921254402,-0.887536753611905,-1.80502182012913,1.40616867953956;0.557597524607164,-0.144824255027548,-0.311715419762457,-1.37261755846291,1.62179047287985,0.201912610104338,0.327437711836565,-1.45236513037498,0.249463982726115,-1.17667517969845,0.283524332217507,1.53324553421279;0.145777576637540,-0.266451854943127,0.552454301267865,0.212671223679860,1.03205658188052,3.70914557998079,0.259351012157144,-0.374424743504582,0.775371422905782,-0.553198254273544,-0.597330556113601,0.454201792899866;-0.0395826495598159,0.794608180473295,-0.509554201431760,-0.0638778580369269,-1.06734925497276,-0.149045791903373,0.0142542593254982,-1.34184440389658,0.465667540071963,-0.433875322006788,-0.566830715943420,1.38391995958180;0.279680226419775,-0.297214568320515,-0.0507042819418351,0.311952991622856,0.444260913613135,0.384561218930201,-0.592291206174866,0.488620435888234,-0.355756064842269,-0.187160185757251,0.265073807527371,-0.710457857184951;-0.366924183263481,-0.822934982563586,-0.755986063777302,-0.526994171669258,-0.243078122480233,0.220260215201429,-0.332140374558679,0.591512032016074,-1.24012082286813,-0.580259827295520,-0.929141410202601,-0.242129539558371;0.746454167316327,0.871021402037080,-0.240285911121369,0.214586292479682,-1.06556347060049,0.0766372190810913,0.0235383813452673,-0.945341840687605,0.450147165191714,-0.515441039946642,-1.05505046121679,2.10317016394157;-0.170426234216238,0.320176745241259,-1.06917778456750,0.0416963150050889,-0.100279326520500,-0.0822318406407962,0.0244605206618037,0.775972033601021,0.316558784411087,-0.359616484134398,0.0851707694557659,-0.137608853039215;-0.642306672790693,1.74897272832719,-1.00627191768658,-0.352706866270512,0.334119809052189,-0.315769414166586,1.25656181101383,-1.10507554249490,1.29831988293240,2.47903000538705,1.43570549129726,-1.25672146344402;0.679737681123251,-1.11448767783986,1.04454103838956,0.239246259731989,-0.509298859953617,-0.662176432790827,0.529269165718851,0.794480119596201,0.488886376479002,0.314544747220011,1.55320039220300,-1.14572178371697;-0.0501203727648018,-0.188320968599843,-0.252136328648327,-0.261490685160765,0.00775675058399299,0.0238506009652036,-0.462115150806582,0.0778364608753722,0.0202359005032030,0.307754557972916,0.241115735535772,-0.360350481811527;-0.0266174311500488,-0.326609146924980,-0.635592131909652,-0.610308691536370,0.256073435293261,0.158070455160251,-1.21708488908643,0.0641388649147828,-0.305885194685064,0.515839794177052,0.813997698809190,-0.977509573629814;0.541489698179866,-1.27040518241144,0.402018446733247,0.0292936162792301,-0.857569279408385,-3.17141066932721,0.285733013498599,0.0877155534177900,-0.276164420043891,0.706991585735083,-0.565774444626175,1.87616030328961;-0.498924910597735,-0.165628738774864,0.837954474447098,0.140718327439906,-0.302510240261667,-0.0132472849954121,-1.32224719963406,-0.248941053857081,-0.863839592735028,-0.893764359347181,0.349833149691466,-2.04609137446243;-0.189519057226046,-0.436619577647587,-0.875275596970663,-0.648807313895567,-0.409452101774158,-0.247226336027258,-0.894298002524404,1.11264318373131,-1.75645363201477,-0.613767559255548,-0.388693131155444,-0.0183424997525372;0.182887453965427,0.361500902721767,0.488092419571214,-0.401763608012290,0.691624628804035,-0.242921855610910,-1.00340854949843,-0.536028241304994,-0.350384377541152,-0.203807209745364,0.401943128340961,-0.178896285046292;0.273060411920776,0.0337140831176831,0.362276373285456,-0.633613709957718,0.421706496945105,0.495641534908431,0.162200407343530,-0.773833085651738,0.406417494781525,-0.199302931819624,0.572556917195535,0.416076157038741;0.435908726351639,-0.510435411783960,-0.193993131046183,0.150585100724567,-0.390262270009095,-0.239409614397035,0.260967344563739,0.310315485258862,-0.133638693940776,0.406182989883801,0.245946932199440,0.649247758624978;0.167355712937677,0.463385284720784,0.418775686153176,0.240795047826126,1.46781963598382,3.81648087762001,-0.0181397513397231,0.323097777988043,-0.663155394101135,-0.470562566025799,-0.0781745625536741,0.227890604832334;0.402291768755644,0.894668502501775,0.425498583343960,-0.545786794322921,0.342353194409575,-0.124214204823823,0.394604489307575,-0.135861272679988,-0.183482926584127,-0.376291395203481,-1.11164314763977,1.09797510522942;0.168624461811107,-0.263478202472342,0.113941364898184,-0.00808483038390651,-0.302329784612151,0.0310527516146670,-1.47332788779261,1.49982645627054,-1.03747712778723,0.808331535215004,-0.328697257779281,0.555778520236409;-0.823298295567540,-1.14474469317069,-0.556968285216460,-0.264457897448577,-1.61880642316903,1.10429861633074,0.303938776944776,0.00570005386838753,-0.0641330505416714,-0.00827048711656718,0.495851063958084,-0.290532235202761;0.223318794950163,-0.489165521747184,-1.39090754937306,-0.0227364916393833,0.672949441745154,-0.537731976548126,0.0219323315867183,0.0861849925809251,0.284485138619561,0.249339567365873,0.714402476392420,0.955520031471213;0.885554651415358,2.54308255104793,1.29355453221198,-0.976632908957834,0.00275861336150128,-1.47483576845646,1.44701357136535,-0.307874199011490,2.24890738941491,-0.939084272656998,-0.890729545274007,1.10994004700744;-0.362868829546417,1.47373824015433,-0.710155090716077,-0.0320949558169643,-0.945920791779683,1.08916940123164,-0.448304228881222,-0.0519160268703096,0.462231119343746,-0.808719494723007,0.0592288244771173,0.512717848054345;0.696386846678340,0.421973718106253,0.699360504395645,0.360608726679624,-0.369781239246742,-0.760796251044894,-0.573501699705008,0.722463663637365,-1.03351150262029,0.299838157167533,0.566403661280218,-0.502374108842679;0.314938812804007,-0.675226058498672,-0.609371465808034,0.175207226851094,-1.73018321780381,0.205312996057938,0.355656922812924,0.000212976610883361,0.469315963502928,0.420442701189783,1.82177067918166,-0.0730837867079938;-0.857636687778917,0.549071433274771,-1.25895102715043,-0.197012418160738,0.218789482076045,1.96517875517588,0.772480979048187,0.861603175776261,-1.41821609044479,-0.244311782532360,-0.997526859316823,-0.244172033073245;-0.699874272439963,0.755915630024105,-0.771907715493350,0.0602697355410549,1.20810280870047,1.48482132111373,-0.180376525997120,-0.175206138925781,-0.251028159981769,-1.41656717861165,0.493114847574199,-1.90747088944521;0.719939717334688,-0.246017973551313,-0.492489753793464,-0.457229938601703,0.733278564758532,1.39946563339902,-0.217816695467535,0.505722796065550,-0.0119210715793176,-0.0710297797717676,0.442572432237226,-0.551357725733626;1.00458387195175,0.104260783176825,-0.0708801595233415,1.38863241765679,0.0832465100113086,0.195692345777590,0.824937029697335,-0.782355665038964,1.08352272577297,-0.693159737349301,-0.763866617843809,1.39307168015003;0.170327931402047,0.255063323385968,3.80528611420830,0.417319514592592,-0.166807710594263,-3.20784737686567,0.434357320788875,-0.130762282207535,0.967327005519535,-0.130053388580474,-0.717800044013087,0.197218965698674;0.0321428455572624,1.43572363047082,-0.0572579176166008,-0.135855616728365,-0.762532074726472,-0.269607949148927,-0.719088648019140,-1.71914381958264,-0.761233159884611,-0.510858139008043,0.469626474630765,-0.292967555751168;-0.0133352222936882,0.0185558735225696,1.09868431848309,-0.176035485285925,0.638879449573260,-0.0775141924359931,-0.200796920813813,0.339392902384433,0.350809224482948,-0.483041720439354,-1.26179724084826,1.45734191011161;0.586612873226813,0.406800978613436,-0.141277323070724,0.120344655499845,-0.928874887409403,-0.679471530640100,0.542992252578423,0.330396445092413,-0.268686498110121,0.977162809203963,0.755636357751808,-0.228211645454660;-0.0429957876556065,0.258310340137851,0.00252832431574437,-0.0575281295785405,0.0270928051998550,0.00584387750355104,-0.0183073274633619,0.0239520899934563,-0.367213599826621,-0.100421915858969,0.392516116117738,-0.584997877251188;-0.246640137144309,0.846528629761747,0.855500139158558,-0.0625816659201595,-0.0933795925527143,0.324236181565059,0.639483847909282,-0.377197922341469,-0.548879002684864,-0.100300689737693,-0.244200167274850,-0.190545214382038;0.0620090561214462,-0.152966376651125,-0.195098277677763,0.317453996895529,-0.581961185500584,-0.0137325137768569,-0.0289047792517975,-0.00493390768311992,-0.00257905169262092,0.212739998700862,1.82850839295711,-2.01351170899743;-0.0828411872253845,1.58162236446053,-0.464414836920566,-0.0911023846256948,-0.579696590001315,-0.461918083247289,-0.931567476652897,-1.36139971311360,-1.00484224103878,-1.17475425557724,1.28570788655571,-1.20148977242997;0.862136691126262,0.263385919874723,0.243179989364402,-0.504385437410524,0.895417045070251,-1.35712852450352,0.774723591054546,0.0603412586991184,0.467353963434823,-0.270363885989131,-0.199721611486192,-0.384938781597238;-0.233379338226818,-1.30705240899365,0.429872133688501,-0.361300497480816,0.0625838778678057,-0.301004989695902,-0.907521624027677,0.519672652253137,-0.604926959058163,0.565458367135622,-0.509042975118740,-0.247568283964934;-1.17041465498012,0.931809405422770,-0.893412489779592,-0.167032895164843,0.319723029313312,0.0564577389737712,-1.90565475458158,1.30571533055684,-0.523958453782178,0.00311348071612935,-0.326700107549195,0.765851547994701];

b1 = [-0.201304400050797;-0.496296613206391;0.104867994947544;0.0293908321013904;-0.303931716596019;-0.0699965305030768;-0.205357376742383;0.176774925960038;-0.585536269376959;-0.589022070739286;-0.980880904956997;1.26729616418371;0.442199338751318;1.46183561501687;1.67692563658122;-0.0768010291538160;0.229515518081084;-1.57534133433423;0.900151408695689;-1.29446382253186;2.23061742436910;0.531133921099591;0.0679691164415854;-1.52141130119197;-0.291795434602107;-0.00650958811431882;-0.00108234013695434;-0.259635398777786;-0.586396477336302;0.0965145152497833;-0.106079248966113;0.386330724537910;-1.18614337727716;-0.385966099229335;3.38382009556592;-0.470741735810751;-0.145230056300641;2.64014021877879;-0.220207688194925;-0.846706270357850;4.14295530586991;-0.0665796376080675;0.235593140227027;0.402289138441661;0.305130950826525;0.186583887341546;2.12275194846626;-0.546818846675797;-0.251824532614526;-0.229753938087427;-3.11640565121104;1.47689474958800;-1.81127833692826;0.0533220607026907;0.149112639690084;-0.313753028971692;4.76461652082208;0.818577859479877;0.516868324206567;-1.34357166229463;-2.39495644711592;0.808330996023989;-0.0762968426958298;1.02219786337630;-1.30920293730550;0.299478455348455;0.895823292932491;1.17422921342329;0.0466295847828411;0.287471860109281;2.56022188322221;0.393509630845357;-0.593585451781347;-0.153609245717998;0.848766100700713;-0.0219843224285430;1.79748848488629;-0.452407530317236;-0.866157648196079;-1.72301239616324];

W2 = [-0.627082721250988,0.646525630327750,0.635520463266840,1.27934819872721,1.04440848313397,1.24656821145702,1.09514164823908,-0.966456041498784,-1.80995187363562,1.06453575567359,-0.446779553835412,-0.657267489204162,-0.852555015406373,-1.45278226300692,-0.988859871409705,0.804040204604263,0.764224039918492,1.61394206098389,-0.870354426893918,1.30455697827047,-0.459857159065043,-0.792592385538374,1.16669268460257,-0.611280840348738,-0.671842505595909,-1.06544893578467,0.502401274472726,0.729907665654012,0.632024007620324,-1.70060353835766,-0.209458266355009,-1.35654955119256,-0.750672185474042,-0.529850756992541,-0.924818470982428,0.988325748742801,0.738191421011367,-1.07531089279312,-1.34326665717061,-0.217721356153622,-1.70722674186021,-1.59354082031977,-0.893396396013901,-0.592822510348920,0.984290347415470,0.645041931947605,1.45681694325273,-1.33625667212834,0.567003461522704,-1.21070803974467,-1.77136833173898,-0.744497575772020,-0.575996462829086,0.798377009300992,1.32273506419465,0.598040550086965,1.24999320927143,1.01203431218493,1.37869494702591,0.997154887475293,-0.716989157644298,0.942834029753209,1.64006949728509,-1.52111124962653,-0.336170212207612,0.199736225360626,-1.73872862458314,1.09382528866244,-0.450957388223983,0.408390733661515,1.98287410656027,1.04592598908956,-0.994554826335448,-0.355665244221539,1.48971681141485,1.72992601596819,-1.24840939820679,-0.672861082606919,0.655951040726704,0.829375867197400;-0.190244929746121,0.430392956885832,0.455867918949123,0.870166417564464,0.965092859740064,1.00191512827027,0.867201105643266,-0.680987823378779,-1.23534314012850,0.709161933407407,-0.313526475500698,-0.453767509210867,-0.676900490284537,-1.03731049480021,-0.738643389765254,0.628931377868789,0.549124606110550,1.13450177332600,-0.590569314158017,0.938805270899072,-0.317435191473693,-0.521805650416264,0.808031164692202,-0.469618724305810,-0.485491331547627,-0.953453004277544,0.381182235461014,0.508874053602029,0.707504666310558,-1.15632761941311,0.492884239168285,-0.960412405315710,-0.527101321528097,-0.0521549430712277,-0.691393600191179,0.696683002574541,0.628798268842944,-0.721056878407010,-0.958945539645820,-0.155285409803611,-1.27332589871223,-1.15148501574431,-0.803009323581546,-0.415895898141473,0.686834764724988,0.392373363476302,1.03390389123209,-0.883333957810329,-0.226927904851245,-0.714775961560234,-1.25679891969819,-0.575114617118888,-0.411145946320113,0.418176993101547,0.902762184128756,0.137152077007692,0.918650550915753,0.661488931603179,1.05803463372901,0.725427640043243,-0.560736112711560,0.617195369265154,1.23154039682301,-1.08210773961181,-0.224131646775166,0.140098689371752,-1.22206421748069,0.748082775389239,-0.320032926824724,0.299354328566732,1.47816858694421,0.717083654714811,-0.755202235247976,-0.781518203472331,1.13944643884742,1.16869798502657,-0.940359202375975,-0.474679724532673,0.465026024175638,0.580819120067879];

b2 = [0.787206208411435;0.628429418007751];

return


%------------------------------------------------------------------------
%         Polynomial coefficient for angular dependency at 1GHz
%------------------------------------------------------------------------
function [pv_ori, ph_ori] = polynomial_1ghz
% polynomial coefficients (deg 2) for the angular dependence at 1.4GHz
pv_ori = [1.689976329544152e-06,4.398401651997119e-04,0.978962600231171];
ph_ori = [-1.216814962390345e-05,1.427972310921177e-04,0.974004864692688];

return

%------------------------------------------------------------------------
%         TELSEM 2 coefficient for angular dependency
%------------------------------------------------------------------------
function [a0_k0, a0_k1, a0_k2, a0_eveh, a1_eveh, a2_eveh, a3_eveh, b0_eveh, b1_eveh, b2_eveh, b3_eveh] = telsem2_coef
%%% Inputs for the angular dependence from TELSEM2_interp_angle.m
%%% Ancillary coefficients, for the 3 frequencies (18, 36, 85GHz) and the 10 surface classes (indicated by class1 in the TELSEM2 atlases):
%%% a0_k0(3,10),a0_k1(3,10),a0_k2(3,10), a0_eveh(3,10),a1_eveh(3,10),a2_eveh(3,10),a3_eveh(3,10),b0_eveh(3,10),b1_eveh(3,10),b2_eveh(3,10),b3_eveh(3,10)
a0_k0= [ 0.11509  0.091535 0.34796 ;...
    0.10525  0.16627  0.24434 ;...
    0.29217  0.23809  0.28954 ;...
    0.17516  0.19459  0.28697 ;...
    0.10521  0.12126  0.30278 ;...
    0.18212  0.19625  0.14551 ;...
    -0.19202  0.5411   0.03739 ;...
    0.10292  0.5486  -0.058937;...
    -0.022672 0.44492 -0.058448;...
    -0.33894 -0.17621  0.14742];
a0_k1= [ 0.61168 0.59095 0.7918 ;...
    0.60271 0.69213 0.62218;...
    0.32728 0.34334 0.37062;...
    0.51217 0.4491  0.50101;...
    0.48913 0.41932 0.29734;...
    0.64474 0.30637 0.031107;...
    1.0405  0.17538 1.3215 ;...
    0.61819 0.31298 1.7218 ;...
    0.87761 0.47583 1.2583 ;...
    1.0959 0.92842 0.51033];
a0_k2= [ 0.26726 0.32033 -0.14778;...
    0.28547 0.13592 0.13193;...
    0.37178 0.41813 0.33875;...
    0.30203 0.35479 0.20189;...
    0.40663 0.47493 0.40668;...
    0.14811 0.52382 0.86634;...
    0.14286 0.27164 -0.37947;...
    0.2737 0.12001 -0.67315;...
    0.13492 0.065463 -0.19316;...
    0.24905 0.25475 0.34637];
a0_eveh=[0.9592599869E+00 0.9565299749E+00 0.9511899948E+00;...
    0.9560700059E+00 0.9541199803E+00 0.9483199716E+00;...
    0.9461100101E+00 0.9439799786E+00 0.9387800097E+00;...
    0.9317600131E+00 0.9289000034E+00 0.9236800075E+00;...
    0.9208700061E+00 0.9190599918E+00 0.9105200171E+00;...
    0.9162799716E+00 0.8937299848E+00 0.8014699817E+00;...
    0.9570500255E+00 0.9213600159E+00 0.7893999815E+00;...
    0.9639400244E+00 0.9530599713E+00 0.8850200176E+00;...
    0.9685299993E+00 0.9622600079E+00 0.9118800163E+00;...
    0.8997200131E+00 0.9012699723E+00 0.9107499719E+00];
a1_eveh=[0.3627802414E-07 -0.7778328204E-08 0.4396108011E-07;...
    0.2503205394E-06 0.1996262995E-06 0.2929977541E-06;...
    0.4190530660E-06 0.3655744649E-06 0.3519195673E-06;...
    0.5574374313E-06 0.5273076340E-06 0.5376484182E-06;...
    0.1026844529E-05 0.9679998811E-06 0.8616486866E-06;...
    0.3180800832E-06 0.2886778532E-06 0.2310362675E-06;...
    -0.1118036366E-06 -0.1502856577E-06 0.4842232926E-07;...
    -0.8410978580E-08 -0.3478669441E-07 0.2209441590E-06;...
    0.2485776633E-06 0.1800235907E-06 0.2510202251E-06;...
    0.2687000915E-06 0.1740325644E-06 0.3562134339E-06];
a2_eveh=[0.3067140824E-05 0.2520012231E-05 0.4831396382E-05;...
    0.8213598448E-05 0.7378375358E-05 0.1022081960E-04;...
    0.1225889173E-04 0.1165553113E-04 0.1188659007E-04;...
    0.1693615741E-04 0.1648317448E-04 0.1715818144E-04;...
    0.2744720041E-04 0.2642072104E-04 0.2671847506E-04;...
    0.1349592094E-04 0.1261523357E-04 0.5447756394E-05;...
    0.2064244654E-05 0.1919016057E-06 0.5940860319E-06;...
    0.5334760772E-05 0.4130339221E-05 0.4104662821E-05;...
    0.6530796327E-05 0.5727014013E-05 0.7451782039E-05;...
    0.1071246970E-04 0.9539280654E-05 0.1034286015E-04];
a3_eveh=[-0.2004991551E-07 -0.6895366056E-07 -0.2047409282E-06;...
    -0.7322448425E-07 -0.1273002681E-06 -0.2729916844E-06;...
    -0.9421125213E-07 -0.1683332300E-06 -0.2726891637E-06;...
    -0.1317753799E-06 -0.2107972250E-06 -0.3556060904E-06;...
    -0.1889465580E-06 -0.2757958271E-06 -0.4909850304E-06;...
    0.7339644004E-08 -0.4058669560E-06 -0.4146343997E-06;...
    0.6170279931E-07 -0.1998567996E-06 -0.4713119139E-07;...
    -0.1361754887E-07 -0.1765622955E-06 -0.2348146637E-06;...
    -0.3901189061E-07 -0.1305666189E-06 -0.1533838798E-06;...
    -0.2679148992E-07 -0.4441960044E-07 -0.1815613899E-06];
b0_eveh=[ 0.9592599869E+00 0.9565299749E+00 0.9511899948E+00;...
    0.9560700059E+00 0.9541199803E+00 0.9483199716E+00;...
    0.9461100101E+00 0.9439799786E+00 0.9387800097E+00;...
    0.9317600131E+00 0.9289000034E+00 0.9236800075E+00;...
    0.9208700061E+00 0.9190599918E+00 0.9105200171E+00;...
    0.9162799716E+00 0.8937299848E+00 0.8014699817E+00;...
    0.9570500255E+00 0.9213600159E+00 0.7893999815E+00;...
    0.9639400244E+00 0.9530599713E+00 0.8850200176E+00;...
    0.9685299993E+00 0.9622600079E+00 0.9118800163E+00;...
    0.8997200131E+00 0.9012699723E+00 0.9107499719E+00];
b1_eveh=[ 0.3626608347E-07 -0.7786279177E-08 0.4393379172E-07;...
    0.2502746099E-06 0.1995944388E-06 0.2929554341E-06;...
    0.4189516289E-06 0.3655020180E-06 0.3518483140E-06;...
    0.5572838404E-06 0.5271903092E-06 0.5375342766E-06;...
    0.1026605219E-05 0.9677979733E-06 0.8614680951E-06;...
    0.3179358714E-06 0.2884899004E-06 0.2308632219E-06;...
    -0.1118781370E-06 -0.1503948681E-06 0.4834672396E-07;...
    -0.8455684153E-08 -0.3485171618E-07 0.2208606134E-06;...
    0.2485595019E-06 0.1799959364E-06 0.2509846695E-06;...
    0.2686167306E-06 0.1739760478E-06 0.3561317214E-06];
b2_eveh=[ 0.3065537157E-05 0.2518960400E-05 0.4829731552E-05;...
    0.8209894986E-05 0.7375769655E-05 0.1021809931E-04;...
    0.1225203869E-04 0.1165053800E-04 0.1188218721E-04;...
    0.1692612022E-04 0.1647546378E-04 0.1715117833E-04;...
    0.2743142431E-04 0.2640772436E-04 0.2670711910E-04;...
    0.1348545720E-04 0.1260529825E-04 0.5439695997E-05;...
    0.2058213340E-05 0.1860650656E-06 0.5898303925E-06;...
    0.5330772183E-05 0.4126528893E-05 0.4100859314E-05;...
    0.6528573977E-05 0.5725009032E-05 0.7449450095E-05;...
    0.1070590315E-04 0.9534271157E-05 0.1033751869E-04];
b3_eveh=[-0.1370247134E-06 -0.1436897747E-06 -0.2954870411E-06;...
    -0.3118435643E-06 -0.2916583242E-06 -0.4311032171E-06;...
    -0.5048401022E-06 -0.4662823869E-06 -0.5206445053E-06;...
    -0.7210980471E-06 -0.6662896794E-06 -0.7548637200E-06;...
    -0.1110204039E-05 -0.1030801400E-05 -0.1140921199E-05;...
    -0.6330818110E-06 -0.9186441048E-06 -0.7947813856E-06;...
    -0.3242539890E-06 -0.5027602583E-06 -0.2777987334E-06;...
    -0.2747250676E-06 -0.3811997260E-06 -0.4102405455E-06;...
    -0.1994112324E-06 -0.2555484855E-06 -0.2842682534E-06;...
    -0.4413041665E-06 -0.3717419474E-06 -0.4975536854E-06];
%%% Transposition of the data (from fortran)
a0_k0=a0_k0';
a0_k1=a0_k1';
a0_k2=a0_k2';
a0_eveh=a0_eveh';
a1_eveh=a1_eveh';
a2_eveh=a2_eveh';
a3_eveh=a3_eveh';
b0_eveh=b0_eveh';
b1_eveh=b1_eveh';
b2_eveh=b2_eveh';
b3_eveh=b3_eveh';
%%%%%%%%%%%%%%%%%

return

%------------------------------------------------------------------------

function [row,col] = latlon_to_ease2(grid_name, lat, lon)
m1=max(size(lat));
m2=max(size(lon));
if m1==m2
    row=zeros(m1,1);
    col=zeros(m1,1);

    [row,col] = wgs84_convert(grid_name, lat, lon);
else
    disp('!! lat lon have different sizes !!')
    return
end

return

%------------------------------------------------------------------------

% ; NAME:
% ;	wgs84_convert
% ;
% ; PURPOSE:
% ;	Use a WGS84 earth model to convert geographic coordinates to
% ;	azimuthal equal-area or cylindrical equal-area grid coordinates
% ;
% ; CATEGORY:
% ;	Grid coordinate conversion
% ;
% ; CALLING SEQUENCE:
% ;       status = wgs84_convert ( grid, lat, lon, r, s )
% ;
% ; INPUTS:
% ;       grid - EASE-Grid-2.0 grid name
% ;       lat, lon - geo. coords. (decimal degrees)
% ;
% ; OUTPUTS:
% ;	r, s - grid row/column coordinates
% ;
% ; RESULT:
% ;	status = 0 indicates normal successful completion
% ;		-1 indicates error status (point not on grid)
% ;
% ; EXAMPLE:
% ;       status = wgs84_convert( 'EASE2_N25km', 90., 0., r, s )
% ;
% ;       status will be 0, and the returned (r, s) will be (359.5, 359.5)

function [r,s] = wgs84_convert(grid_name, lat, lon)

[grid]=ease2_grid_info( grid_name);

projection=grid.projection;
map_scale_m=grid.map_scale_m;
r0=grid.r0;
s0=grid.s0;


[x,y]=wgs84_convert_xy( projection, lat, lon);

% Convert xy coordinates to row/col for this particular grid
r = round(r0 + ( x./ map_scale_m ) +1) ;
s = round(s0 - ( y./ map_scale_m ) +1) ;
%We add one because Matlab matrix start at 1,1 but in EASE-Grid the first pixel is 0,0

return

%------------------------------------------------------------------------

% ; NAME:
% ;	wgs84_convert_xy
% ;
% ; PURPOSE:
% ;	Use a WGS84 earth model to convert geographic coordinates to
% ;	azimuthal equal-area or cylindrical equal-area projection coordinates
% ;
% ; CATEGORY:
% ;	Projection coordinate conversion
% ;
% ; CALLING SEQUENCE:
% ;       status = wgs84_convert_xy ( projection, lat, lon, x, y)
% ;
% ; INPUTS:
% ;       projection: projection ['N','S', or 'M']
% ;             N = Northern azimuthal equal area
% ;             S = Southern azimuthal equal area
% ;             M = cylindrical equal area
% ;       lat, lon - geo. coords. (decimal degrees)
% ;
% ; OUTPUTS:
% ;	x, y - projection x, y coordinates (meters)
% ;
% ; RESULT:
% ;	status = 0 indicates normal successful completion
% ;		-1 indicates error status (point not on grid)
% ;
% ; EXAMPLE:
% ;       status = wgs84_convert_xy( 'N', 90., 0., x, y )
% ;
% ;       status will be 0, and the returned (x, y) will be (0.0, -0.0)
% ;
% ; REFERENCE:
% ;
% ; Snyder, John P. Map Projections--A Working Manual. USGS Professional
% ; Paper 1395. US Govt Printing Office, Washington, DC. 1987.
%
% ; Lambert Azimuthal equal-area projection formulas for the ellipsoid:
% ; pp. 187-190.
% ;
% ; Cylindrical equal-area projection formulas for the ellipsoid: pp. 81-85.
% ;
function [x,y] = wgs84_convert_xy(projection, lat, lon)



%;; epsilon for test in neighborhood of pole for azimuthal projections
epsilon = 1.0D-6;
map = ease2_map_info(projection);


dlon = lon - map.map_reference_longitude;
dlon = easeconv_normalize_degrees( dlon );
phi = easeconv_deg2rad( lat );
lam = easeconv_deg2rad( dlon );

sin_phi = sin( phi );

q = ( 1.0 - map.e2 ).* ( ( sin_phi./ ( 1.0 - map.e2.* sin_phi.* sin_phi ) )- ( 1.0./ ( 2.0.* map.map_eccentricity ) ) * log( ( 1.0 - map.map_eccentricity.* sin_phi )./ ( 1.0 + map.map_eccentricity.* sin_phi ) ) );

qp = 1.0 - ( ( 1.0 - map.e2 )./ ( 2.0.* map.map_eccentricity )* log( ( 1.0 - map.map_eccentricity )./ ( 1.0 + map.map_eccentricity ) ) );

switch upper(projection)
    case 'N'
        if ( abs( qp - q ) < epsilon )

            rho = 0.0;

        else

            rho = map.map_equatorial_radius_m.* sqrt( qp - q );

        end

        x = rho.* sin( lam );
        y = -1.0.* rho.* cos( lam );

    case 'S'

        if ( abs( qp + q ) < epsilon )

            rho = 0.0;

        else

            rho = map.map_equatorial_radius_m.* sqrt( qp + q );

        end

        x = rho.* sin( lam );
        y = rho.* cos( lam );

    case 'M'

        x = map.map_equatorial_radius_m.* map.kz.* lam;
        y = ( map.map_equatorial_radius_m.* q )./ ( 2.0 * map.kz );

    otherwise
        printf("Programming error")
        return
end

return

%------------------------------------------------------------------------

% ; NAME:
% ;	ease2_grid_info
% ;
% ; PURPOSE:
% ;       Returns grid parameters for the specified EASE-Grid-2.0 grid
% ;
% ; CALLING SEQUENCE:
% ;       status = ease2_grid_info( grid )
% ;
% ; ARGUMENTS:
% ;       grid - EASE-Grid-2.0 grid name, assumes pattern is 'EASE2_P',
% ;              where P=projection='N', 'S' or 'M' e.g. EASE2_M25km
% ;
% ; KEYWORDS:
% ;       MAP_SCALE_M - map scale (meters)
% ;       COLS - grid columns
% ;       ROWS - grid rows
% ;       R0 - column that is mapped to map_reference_latitude
% ;       S0 - row that is mapped to map_reference_latitude
% ;
% ; RESULT:
% ;       0 - successful completion
% ;       -2 - undefined grid name
% ;
% ; EXAMPLE:
% ;       status = ease2_grid_info( 'EASE2_N25km', COLS=cols, ROWS=rows )
% ;
% ;       will return status of 0 and cols = 720, rows = 720
% ;
function [grid]=ease2_grid_info( grid_name)%                           PROJECTION=projection, $
%                           MAP_SCALE_M=map_scale_m, $
%                           COLS=cols, $
%                           ROWS=rows, $
%                           R0=r0, $
%                           S0=s0

grid=struct;
switch upper(grid_name)
    case 'EASE2_N25KM'
        grid.projection='N';
        grid.map_scale_m=25000.0;
        grid.cols=720;
        grid.rows=720;
        grid.r0=0.0;
        grid.s0=0.0;

    case 'EASE2_N12.5KM'
        grid.projection='N';
        grid.map_scale_m=12500.0;
        grid.cols=1440;
        grid.rows=1440;
        grid.r0=0.0;
        grid.s0=0.0;
    case 'EASE2_N6.25KM'
        grid.projection='N';
        grid.map_scale_m=6250.0;
        grid.cols=2880;
        grid.rows=2880;
        grid.r0=0.0;
        grid.s0=0.0;
    case 'EASE2_N36KM'
        grid.projection='N';
        grid.map_scale_m=36000.0;
        grid.cols=500;
        grid.rows=500;
        grid.r0=0.0;
        grid.s0=0.0;
    case 'EASE2_N09KM'
        grid.projection='N';
        grid.map_scale_m=9000.0;
        grid.cols=2000;
        grid.rows=2000;
        grid.r0=0.0;
        grid.s0=0.0;
    case 'EASE2_N03KM'
        grid.projection='N';
        grid.map_scale_m=3000.0;
        grid.cols=6000;
        grid.rows=6000;
        grid.r0=0.0;
        grid.s0=0.0;
    case 'EASE2_N01KM'
        grid.projection='N';
        grid.map_scale_m=1000.0;
        grid.cols=18000;
        grid.rows=18000;
        grid.r0=0.0;
        grid.s0=0.0;
    case 'EASE2_N100KM'
        grid.projection='N';
        grid.map_scale_m=100000.0;
        grid.cols=180;
        grid.rows=180;
        grid.r0=0.0;
        grid.s0=0.0;

    case 'EASE2_S25KM'
        grid.projection='S';
        grid.map_scale_m=25000.0;
        grid.cols=720;
        grid.rows=720;
        grid.r0=0.0;
        grid.s0=0.0;

    case 'EASE2_S12.5KM'
        grid.projection='S';
        grid.map_scale_m=12500.0;
        grid.cols=1440;
        grid.rows=1440;
        grid.r0=0.0;
        grid.s0=0.0;
    case 'EASE2_S6.25KM'
        grid.projection='S';
        grid.map_scale_m=6250.0;
        grid.cols=2880;
        grid.rows=2880;
        grid.r0=0.0;
        grid.s0=0.0;
    case 'EASE2_S36KM'
        grid.projection='S';
        grid.map_scale_m=36000.0;
        grid.cols=500;
        grid.rows=500;
        grid.r0=0.0;
        grid.s0=0.0;
    case 'EASE2_S09KM'
        grid.projection='S';
        grid.map_scale_m=9000.0;
        grid.cols=2000;
        grid.rows=2000;
        grid.r0=0.0;
        grid.s0=0.0;
    case 'EASE2_S03KM'
        grid.projection='S';
        grid.map_scale_m=3000.0;
        grid.cols=6000;
        grid.rows=6000;
        grid.r0=0.0;
        grid.s0=0.0;
    case 'EASE2_S01KM'
        grid.projection='S';
        grid.map_scale_m=1000.0;
        grid.cols=18000;
        grid.rows=18000;
        grid.r0=0.0;
        grid.s0=0.0;
    case 'EASE2_M25KM'
        grid.projection='M';
        grid.map_scale_m=25025.2600081;
        grid.cols=1388;
        grid.rows=584;
        grid.r0=0.0;
        grid.s0=0.0;
    case 'EASE2_M12.5KM'
        grid.projection='M';
        grid.map_scale_m=12512.63000405;
        grid.cols=2776;
        grid.rows=1168;
        grid.r0=0.0;
        grid.s0=0.0;
    case 'EASE2_M6.25KM'
        grid.projection='M';
        grid.map_scale_m=6256.315002025;
        grid.cols=5552;
        grid.rows=2336;
        grid.r0=0.0;
        grid.s0=0.0;
    case 'EASE2_M36KM'
        grid.projection='M';
        grid.map_scale_m=36032.220840584;
        grid.cols=964;
        grid.rows=406;
        grid.r0=0.0;
        grid.s0=0.0;
    case 'EASE2_M09KM'
        grid.projection='M';
        grid.map_scale_m=9008.055210146;
        grid.cols=3856;
        grid.rows=1624;
        grid.r0=0.0;
        grid.s0=0.0;
    case 'EASE2_M03KM'
        grid.projection='M';
        grid.map_scale_m=3002.6850700487;
        grid.cols=11568;
        grid.rows=4872;
        grid.r0=0.0;
        grid.s0=0.0;
    case 'EASE2_M01KM'
        grid.projection='M';
        grid.map_scale_m=1000.89502334956;
        grid.cols=34704;
        grid.rows=14616;
        grid.r0=0.0;
        grid.s0=0.0;
    otherwise
        fprintf('unknown grid projection')
        return
end


grid.r0 = ( grid.cols -1)./ 2.0;
grid.s0 = ( grid.rows -1)./ 2.0;

return

%------------------------------------------------------------------------

% ; NAME:
% ;       ease2_map_info
% ;
% ; PURPOSE:
% ;       Returns map parameters for the specified EASE-Grid-2.0 projection;
% ;       by definition, the EASE-Grid-2.0 projections use the WGS84 ellipsoid
% ;
% ; CALLING SEQUENCE:
% ;       status = ease2_map_info( projection )
% ;
% ; ARGUMENTS:
% ;       projection - EASE-Grid-2.0 projection, one of 'N', 'S' or 'M', where:
% ;                    N - Northern Lambert azimuthal equal-area
% ;                    S - Southern Lambert azimuthal equal-area
% ;                    M - Global cylindrical equal-area
% ;
% ; KEYWORDS:
% ;       MAP_EQUATORIAL_RADIUS_M - WGS84 equatorial radius
% ;       MAP_ECCENTRICITY - WGS84 eccentricity
% ;       E2 - eccentricity^2
% ;       MAP_REFERENCE_LATITUDE - map reference latitude (degrees)
% ;       MAP_REFERENCE_LONGITUDE - map reference longitude (degrees)
% ;       MAP_SECOND_REFERENCE_LATITUDE - phi1 = true latitude (only for
% ;                                       projection='M')
% ;       SIN_PHI1 - sin( MAP_SECOND_REFERENCE_LATITUDE ) - sin(phi1) (only
% ;                  for projection='M')
% ;       COS_PHI1 - cos( MAP_SECOND_REFERENCE_LATITUDE ) - cos(phi1) (only
% ;                  for projection='M')
% ;       KZ - cos(phi1) / sqrt(1 - e^2 sin^2(phi1)) (only for projection 'M')
% ;
% ; RESULT:
% ;       0 - successful completion
% ;       -2 - undefined projection name
% ;
% ; EXAMPLE:
% ;       status = ease2_map_info( 'N', MAP_REFERENCE_LATITUDE=lat0,
% ;       MAP_REFERENCE_LONGITUDE=lon0 )
% ;
% ;       will return status of 0 and lat0=90.0D and lon0=0.0D
% ;
function map=ease2_map_info( projection)

map.map_equatorial_radius_m = 6378137.0; %WGS84
map.map_eccentricity = 0.081819190843;   %WGS84
map.e2 = map.map_eccentricity^2.0;

switch upper(projection)
    case 'N'

        %Parameters specific to this projection
        map.map_reference_latitude = 90.0;
        map.map_reference_longitude = 0.0;

    case 'S'

        %Parameters specific to this projection
        map.map_reference_latitude = -90.0;
        map.map_reference_longitude = 0.0;

    case 'M'

        %Parameters specific to this projection
        map.map_reference_latitude = 0.0;
        map.map_reference_longitude = 0.0;
        map.map_second_reference_latitude = 30.0;
        map.sin_phi1 = sind( map.map_second_reference_latitude );
        map.cos_phi1 = cosd( map.map_second_reference_latitude );
        map.kz = map.cos_phi1./ sqrt( 1.0 - map.e2.* map.sin_phi1.* map.sin_phi1 );

    otherwise
        fprintf("Programming error.")
        return
end

return

%------------------------------------------------------------------------
% ; NAME:
% ;	easeconv_deg2rad
% ;
% ; PURPOSE:
% ;       Converts input value from degrees to radians using double-precision value of pi
% ;
% ; CALLING SEQUENCE:
% ;       rad = easeconv_deg2rad( angle )
% ;
% ; ARGUMENTS:
% ;       angle - angle in degrees
% ;
% ; RESULT: angle in radians
% ;
% ; EXAMPLE:
% ;       angle = 90.
% ;       rad = easeconv_deg2rad( angle )
% ;
% ;       returns 1.5707963D
% ;
function [angled]=easeconv_deg2rad(angle)

angled = angle.* pi./ 180.0;

return



%------------------------------------------------------------------------
% ; NAME:
% ;       easeconv_normalize_degrees
% ;
% ; PURPOSE:
% ;       Normalizes input value to range [-180, 180]
% ;
% ; CALLING SEQUENCE:
% ;       norm_lon = easeconv_normalize_degrees( lon )
% ;
% ; ARGUMENTS:
% ;       lon - longitude (decimal degrees)
% ;
% ; RESULT: normalized longitude
% ;
% ; EXAMPLE:
% ;       lon = 185.
% ;       lon = easeconv_normalize_degrees( lon )
% ;
% ;       Returned value will be -175.

function [lon]=easeconv_normalize_degrees(inlon)

lon = inlon ;

while ( single(lon) < -180. )
    lon = lon + 360.;
end

while ( single(lon) > 180. )
    lon = lon - 360.;
end

return

%------------------------------------------------------------------------

function yd = tge_dayofyear(varargin)
%DAYOFYEAR Ordinal number of day in a year.
%
%   DAYOFYEAR(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND) returns the ordinal
%   day number in the given year plus a fractional part depending on the
%   time of day.
%
%   Any missing MONTH or DAY will be replaced by 1.  HOUR, MINUTE or SECOND
%   will be replaced by zeros.
%
%   If no date is specified, the current date and time is used.  Gregorian
%   calendar is assumed.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 12:52:04 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

nargsin = nargin;
error(nargchk(0, 6, nargsin));
if nargsin
    argv = { 1 1 1 0 0 0 };
    argv(1:nargsin) = varargin;
else
    argv = num2cell(clock);
end
[year, month, day, hour, minute, second] = deal(argv{:});

days_in_prev_months = [0 31 59 90 120 151 181 212 243 273 304 334];


% Day in given month.
yd = days_in_prev_months(month)' ...               % days in prev. months
    + ( tge_isleapyear(year) & ( month > 2 ) ) ...   % leap day
    + day ...                                    % day in month
    + ( second + 60*minute + 3600*hour )/86400;  % part of day

return

%------------------------------------------------------------------------




% tge_doy_10day_composite
%
%   It returns the middle, starting, and ending day of year (DOY)
%   for the 36 10-day time intervals used in the BIOMAP database. These
%   are strictly not 10-days for some intervals, as they are adjusted
%   to facilitate monthly means, so 3 intervals make up a month.
%
%   If no input argument is given, it outputs the middle day of the
%   intervals [docm]. This will be used to label the database files,
%   and it is common for all database years, independent of having
%   a leapyear or not.
%
%   If [iyear] is given, it will outputs as well the start [docs] and
%   end [doce] of the 36 intervals.
%
%   If [doy] is also fiven, in this case [docm, docs, doce] are the
%   middle, starting, anf ending day of the interval corresponding to
%   that [doy].

%   Author:      C. Jimenez
%   Time-stamp:  2023-10-29
%   E-mail:      carlos.jimenez@estellus.fr

function [docm, docs, doce] = tge_doy_10day_composite( iyear, doy )

docs = [1 11 21    32 42 52   60 70 80  91 101 111   121 131 141   152 162 172   182 192 202    213 223 233   244 254 264      274 284 294   305 315 325   335 345 355 ];

docm = docs+4;

if nargin > 0

    isleap = ( ~rem(iyear, 4) & rem(iyear, 100) ) | ~rem(iyear, 400);

    if isleap
        docs = [ docs(1:6) 1+docs(7:end) ];
    end

    doce = docs-1;
    if isleap
        doce = [ doce(2:end) 366];
    else
        doce = [ doce(2:end) 365];
    end


    if nargin == 2

        if doy<1 | doy>366
            error('Doy is incorrect')
        end

        ind=find((docs-doy)<=0);
        docs = docs(ind(end));

        ind=find((doce-doy)>=0);
        doce = doce(ind(1));

        [~,docm] = oge_nearestinvec(docm,doy);

    end

end

return

%------------------------------------------------------------------------

function t = tge_isleapyear(year)
%ISLEAPYEAR True for leap years.
%
%   ISLEAPYEAR(YEAR) returns 1's for the elements of YEAR that are leap
%   years and 0's for those that are not.  If YEAR is omitted, the current
%   year is used.  Gregorian calendar is assumed.
%
%   A year is a leap year if the following returns true
%
%       ( ~rem(year, 4) & rem(year, 100) ) | ~rem(year, 400)
%
%   A year is not a leap year if the following returns true
%
%      rem(year, 4) | ( ~rem(year, 100) & rem(year, 400) )

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 12:51:45 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

error(nargchk(0, 1, nargin));

if nargin == 0               % If no input argument...
    clk = clock;              % ...get current date and time...
    year = clk(1);            % ...and extract year.
end

t = ( ~rem(year, 4) & rem(year, 100) ) | ~rem(year, 400);

return

%-----------------------------------------------------------------


function [ind,zn] = oge_nearestinvec(z,zi)
% [ind,zn] = oge_nearestinvec(z,zi)
%
%          Find the nearest value in a vector z to a given value zi
%
% FORMAT   [ind,zn] = oge_nearestinvec(z,zi)
%
% OUT      ind  the index in vector z of closest value to zi
%          zn   the closest value in z to zi
% IN       z    the vector
%          zi   the value

% 2003-02-10   Created by Carlos Jimenez
[ aux, ind ] = min( abs( z -zi ) );
zn           = z( ind );

return

%------------------------------------------------------------------------


function [utc, utcdaydiff] = tools_loctime2utc(lt,lon)
%TOOLS_UTC2LOCTIME Conversion from local time to UTC
%   Detailed explanation goes here
if isdatetime(lt)
    lt = hour(lt) + minute(lt)/60;
end
utc = lt - lon*12/180;
utcdaydiff = zeros(size(utc));

ind = find( utc >= 24);
utc(ind)      = utc(ind) - 24;
utcdaydiff(ind) = 1;


ind = find( utc < 0);
utc(ind)      = utc(ind) + 24;
utcdaydiff(ind) = -1;

return




%------------------------------------------------------------------------

function [F,TF] = fillmissing2_manual(A, method, missing_data_to_fill)
% Method must be in :
%        "nearest"  - Nearest (in terms of index) non-missing entry.
%        "linear"   - Linear interpolation of non-missing entries.
%        "natural"  - Natural neighbor interpolation of non-missing entries.
%        "cubic"    - Cubic interpolation of non-missing entries.
%        "v4"       - Biharmonic spline interpolation of non-missing entries.
arguments
    A (:,:)
    method char {mustBeMember(method,{'nearest','linear','natural', 'cubic', 'v4'})} = 'linear' 
    missing_data_to_fill (:,:) = ismissing(A); % Fill only these locations
end
if size(A)~=size(missing_data_to_fill)
    error("A and missing_data should have the same size")
end
missing_data = ismissing(A);
allMissingLocations = missing_data;


sp = {1:size(A,1),1:size(A,2)};

F = A;
TF = false(size(A));

% Quick return when there are no missing values or empties
if ~any(missing_data_to_fill,"all") || isempty(A)
    return
end

[X,Y] = find(~allMissingLocations);
if isempty(X)
    % Quick return when all values are missing
    return
end
X = sp{1}(X);
Y = sp{2}(Y);
[xq,yq] = find(missing_data_to_fill);
xq = sp{1}(xq);
yq = sp{2}(yq);

fillVals = griddata(X,Y,double(A(~allMissingLocations)),xq,yq,method);
if ~isempty(fillVals)
    F(missing_data_to_fill) = fillVals;
    TF(missing_data_to_fill) = ~ismissing(fillVals);
end

return



