% y  vector freq1 v-pol h-pol freq 2 v-pol hpol .....
%    output for standard retrievasl using v and h, no plans
%    for the moent to use t3 and t4
%
% yf is the output with harmonics for other applications
%    e.g. rt simulations.
%
% emis to output the emis

%               This function is implemented as a general forward 
%		model for different retrieval aplications in the 
%		fra, and it has
%		been developed to be used by different projects
%		including the ESA SCEPS contract number 4000138129/22/NL/IA
%		and the Study on CIMR synergistic global ocean and atmospheric 
%		products, contract number ITT 23/224353.


function [bts, bts_harm_iso, bts_harm_cs1, bts_harm_cs2, esea, eice, elan, tra, tbup, tbdown, TCWV, TCLRW ] = fom_atmos_ocean_seaice_land_tbs_core( atmos_input, surf_input, sensor_input, dir_input, do_iso );


%=== externally calculated cos and sin

if isfield( surf_input, 'COSOWD' )
  do_ecos = 1;
else
  do_ecos = 0;
end




%=== only one azimuth angle is allowed
%    so the t3 and t4 correspond to the 
%    relative angle owd - aa

%if length( sensor_input.AA ) > 1
%  error('Only one azimuth angle is allowed')
%end


%=== if FM used for retrievals so far
%    only the total tbs Vpol and Hpol 
%    need to be output, the harmonic
%    calculations and t3 and t4 are not
%    output and format is tbs stored as
%    v-nza1-fr1, h-nza1-fr1, v-nza1-fr2, h-nza1-fr2 .... 


if nargout == 1
  do_harm = 0;
else
  do_harm = 1;
end



%=== If MASK is one value mot mixed scenes allowed 
%    apart from sea-ice and ocean when SIC < 1. 
%    Reminder: MASK 1 ocean 2 land 3 sea-ice

if length( surf_input.MASK ) == 1

  pc_sea = 0;
  pc_ice = 0;
  pc_lan = 0;

  % sea

  if surf_input.MASK == 1
    pc_sea = 1;
  % sea-ice
  elseif surf_input.MASK == 2
    pc_lan = 1;
  elseif surf_input.MASK == 3 & surf_input.SIC > 0 
    pc_ice = surf_input.SIC;;
    pc_sea = 1 - pc_ice;
  % lan

  end


else

%=== If MASK is a 3x1 vector it is used as
%    precentages of ocean, sea-ice, and land
%    and the SIC value if exist is ignored.

  pc_sea = surf_input.MASK(1);
  pc_lan = surf_input.MASK(2);
  pc_ice = surf_input.MASK(3);

  
end


%=== overwriting pc_sea if seaice emis
%    already contains ocean contribution.
%    This is a tmp solution, see below
%    in the ice emis calculation



if isfield( surf_input, 'SIMONTH' ) & pc_ice > 0
  pc_sea = 0;
end



nch = length( sensor_input.F );
nza = length( sensor_input.ZA );
naa = length( sensor_input.AA );


if sum( sensor_input.VPOL + sensor_input.HPOL ) == 2 * nch

  do_pol = 0;

else

  do_pol = 1;
  aux  = [];
  for a = 1:nza
    aux    = [aux sensor_input.VPOL sensor_input.HPOL];
  end
  indpol = find( aux == 1);

end



%= atmosphere

if sum( sensor_input.TOA ) > 0 

  [ tra, tbdown, tbup, TCWV, TCLRW, ~, ~, ~, tbc ] = fom_atmos_tbs_rosenkranz( sensor_input.ZA, sensor_input.F, sensor_input.H, atmos_input.P, atmos_input.SP, atmos_input.T, atmos_input.Q, atmos_input.CLRWC, atmos_input.LOWFREQ, atmos_input.TCWV, atmos_input.TCLRW );

  %= calcuation result in a matrix fre x nza
  %  if nza == nch, i.e., each frequency at a
  %  za we take the diagonal to have the calculation 
  %  each fre-za


  if nza == nch

    tra    = diag(tra);
    tbdown = diag(tbdown);
    tbup   = diag(tbup);
    tbc    = diag(tbc);

  end 

  %= if not the atmosphere for all channels
  if sum( sensor_input.TOA ) < length( sensor_input.F )

    aux     = repmat(sensor_input.TOA',1,nza);
    tra     = tra.* aux;
    tbdown  = tbdwon .* aux;
    tbup    = tbup .* aux;

  end


  %= preparing for contribution calculation 

  if ~do_harm

    %= no harmonics 
    %  v-fr1-za1 h-fr1-za1 v-fr2-za1 ..... v-fr1-za2 h-fr1-za2 .....
     

    atra    = [];
    atbdown = [];
    atbup   = [];

    % further check this bug
    % for a = 1:nza
    for a = 1:1

      aux      = repmat(tra(:,a),1,2)'; 
      atra     = [ atra; aux(:)];
      aux      = repmat(tbdown(:,a),1,2)'; 
      atbdown  = [ atbdown; aux(:)];
      aux      = repmat(tbup(:,a),1,2)'; 
      atbup    = [ atbup;aux(:)];

    end

  else

    %= with harmonics format  
    %  nza x (4 (v-h-t3-t4) x nfreqs ) 
    % using original tra, tbdown, 
    % modifying tbup for total contribution

    atbup = zeros( 4, nch, nza );
    for a = 1:2
      atbup(a,:,:) = tbup;
    end      
   
  end


else

  % no atmosphere for all freqs

  if ~do_harm

    atbdown = zeros( nch * nza * 2, 1 );  
    atbup   = zeros( nch * nza * 2, 1 );
    atra    = ones( nch * nza * 2, 1 ); 
    TCWV    = 0;
    TCLRW   = 0;

  else

    atra    = ones(nza,4,nch);
    atbdown = zeros(nza,4,nch);
    atbup   = atbdown;
    TCWV    = 0;
    TCLRW   = 0;    

  end 

   

end 





%============================= land surface

if pc_lan > 0

  do_land = 1;

  if ~isfield( surf_input, 'LEMIS' )

    %===  No external emis so calculations 

    if isfield( surf_input, 'SND' )

      %=== originally using a climato for land, it can still be used by not having a snow
      %    density (SND) field as input. If using the current land one, notice that it is
      %    still only valid for snowed covered land


      ifreq = sensor_input.F;
      ifreq( ifreq == 1.4135 ) = 1.41;

      path_rsrc = [ dir_input, '/EmisSnow' ];

      %= original with climato and previoulsy calculated and daily for SNA and SNTL1 
      %  for SNA, SNTL1. Using a nan value in SNT for cases where we just want
      %  to use the climato interpolated in angle

      if ~isnan( surf_input.SNTH )

        %= dynamical emissivity

        elan = fom_snowland_emis_surfem( ifreq, sensor_input.ZA, surf_input.LAT, surf_input.LON, surf_input.MONTH, surf_input.DAY, surf_input.LTIME, path_rsrc, surf_input.LC , surf_input.SDOR , surf_input.DSNA , surf_input.SND , surf_input.SNTH , surf_input.LST , surf_input.SNT , surf_input.DSNTL1 , surf_input.SNMLT, surf_input.CLEMIS, surf_input.CLAGB );

      else

       %= climato with angular dependence 

        elan = fom_snowland_emis_surfem( ifreq, sensor_input.ZA, surf_input.LAT, surf_input.LON, surf_input.MONTH, surf_input.DAY, surf_input.LTIME, path_rsrc, [] , [] , [] , surf_input.SND , [] , [] , [] , [] , [], surf_input.CLEMIS, surf_input.CLAGB );

      end

    else

      elan = emis_land_angle_interpolate( surf_input.LON, surf_input.LAT, surf_input.MONTH, dir_input, sensor_input.F, sensor_input.ZA  );
 
    end

  else

    %=== using passed value

    elan = surf_input.LEMIS; 

  end


  %= using climato till SCEPS-II parameterizaion is implemented


  if ~do_harm
  
    elan = elan(2,:,:);
    aelan = reshape( elan, 2, nch * nza);
    aelan = aelan(:); 

    ts       = surf_input.LST *ones( 2*nch*nza, 1);  
    bts_lan  = atbdown .* atra .* ( 1-aeland ) + ts .* atra .* aeland;

  else
  
    ts       = surf_input.LST *ones( nch, nza );  

    bts_lan_harm_iso        = zeros( 4, nch, nza );
    bts_lan_harm_iso(1,:,:) = tbdown .* tra .* ( 1 - squeeze(elan(1,:,:)))  + ts .* tra .* squeeze(elan(1,:,:));
    bts_lan_harm_iso(2,:,:) = tbdown .* tra .* ( 1 - squeeze(elan(2,:,:)))  + ts .* tra .* squeeze(elan(2,:,:));

    bts_lan_harm_cs1 = zeros( 4, nch, nza );
    bts_lan_harm_cs2 = zeros( 4, nch, nza );


  end


else

  if ~do_harm

    if nza == 1
      % common angle all freqs
      bts_lan = zeros( nch * nza * 2, 1 );
    else
      % individual angle per freq
      bts_lan = zeros( nza * 2, 1 );
    end

  else

    elan            = nan( 2, nch, nza );
    bts_lan_harm_iso = zeros( 4, nch, nza );
    bts_lan_harm_cs1 = zeros( 4, nch, nza );
    bts_lan_harm_cs2 = zeros( 4, nch, nza );

  end

end



%==========================  sea surface
%			     also needed for sea ice if SIC < 1


if pc_sea > 0


  %=== emis calculation, preparing
  %    matrices for freq x nza 

  %=== common za and aa from all channels

  if nza == 1

    %=== emis calculation, preparing
    %    matrices for freq x nza 

    auc  = ones(nch,nza);

    freq = auc;
    for c = 1:nch
      freq(c,:) = sensor_input.F(c);
    end
    freq = freq(:)';

    angle = auc;
    for c = 1:nza
      angle(:,c) = sensor_input.ZA(c);
    end
    angle = angle(:)'; 
   
    auo = ones(1,nch*nza);
    sst = surf_input.SST * auo;
    sss = surf_input.SSS * auo;
    ows = surf_input.OWS * auo;

    %=   relative to satellite azimuth direction defined with the same convention
    %    CONVENTION: angle between the ground projection of the emission direction 
    %    vector and the east direction, counterclockwise positive

    rowd = wrapTo180(surf_input.OWD - sensor_input.AA) * auo;

    itra = auc;
    if naa == 1
      itra = tra;
    else
      for c = 1:naa
        itra(:,:) = tra;
      end
    end
    itra = itra(:)';


    %= emis esea and refl_sea are output as (4 [v-h-t3-t4], nch * nza) with nch * nza 
    %  as fr1-z1, fr2-za1 ... fr1-za2
  
    %=  emis_sea_harmonics components as (1, nch * nza) with  fr1-z1, fr2-za1 ... fr1-za2)

    [~, emis_sea_harmonics, z_sea ] = fom_ocean_emis_surfem( sst, sss, ows, rowd, freq, angle, itra);


  else

    %=== individual za and aa angles per channel

    if 0   % old before new SCEPS2, will ti work the new for retrievals?

    freq  = sensor_input.F;
    angle = sensor_input.ZA;

    auo  = ones(1,nch);
    sst  = surf_input.SST * auo;
    sss  = surf_input.SSS * auo;
    ows  = surf_input.OWS * auo;
    rowd = wrapTo180(surf_input.OWD - sensor_input.AA) * auo;

    end

    auo  = ones(1,nch * nza);

    sst  = surf_input.SST * auo;
    sss  = surf_input.SSS * auo;
    ows  = surf_input.OWS * auo;
    rowd = wrapTo180(surf_input.OWD - sensor_input.AA) * auo;

    freq = [];
    for f = 1:nch
      freq = [ freq sensor_input.F(f)*ones(1,nza) ];
    end 

    angle = [];
    for a = 1:nch
      angle = [ angle sensor_input.ZA ];
    end 

    itra = tra';
    itra = tra(:)'; 

    %=  emis_sea_harmonics components as (1, nch * nza) with  fr1-za1, fr2-za1 ... fr1-za2)

    [~, emis_sea_harmonics, z_sea ] = fom_ocean_emis_surfem( sst, sss, ows, rowd, freq, angle, itra);


  end


  %=  when no harmonics are computed the bts is always the full tb
  %   with dependence on azimuth angle, i.e., isotropic and
  %   anisotropic components. To only get the isotropic the user
  %   needs to reads the bts_harm_iso. For instance, when using the  
  %   the FM for retrievals of TOA-bts, when there is still no observation
  %   i.e., no observing azimuth angle, the bts_harm_iso needs to be
  %   used. But if we are retrieving instrument bts with a certain eaa
  %   the full bts is more approriate if the rowd is known to some degree.

  if ~do_harm

    %= no harmonics 
    %  v-fr1-za1 h-fr1-za1 v-fr2-za1 ..... v-fr1-za2 h-fr1-za2 .....
   
    %= without anisotropic 

    if do_iso == 1

      aesea     = [ (emis_sea_harmonics.evn + emis_sea_harmonics.ev0  )' (emis_sea_harmonics.ehn + emis_sea_harmonics.eh0  )' ]'; 

    else

      if do_ecos 

        coseaa = cosd(sensor_input.AA');
        sineaa = sind(sensor_input.AA');
        cosowd = surf_input.COSOWD * auo;
        sinowd = surf_input.SINOWD * auo;

        cosrowd    = cosowd .* coseaa + sinowd .* sineaa; 
        sinrowd    = sinowd .* coseaa - cosowd .* sineaa;

        cos2rowd   = cosrowd.^2 - sinrowd.^2; 
        %sin2rowd  = 2 * sinrowd .* cosrowd;

        aesea     = [ (emis_sea_harmonics.evn + emis_sea_harmonics.ev0 + emis_sea_harmonics.ev1 .* cosrowd + emis_sea_harmonics.ev2 .* cos2rowd )' (emis_sea_harmonics.ehn + emis_sea_harmonics.eh0 + emis_sea_harmonics.eh1 .* cosrowd + emis_sea_harmonics.eh2 .* cos2rowd )' ]'; 

      else 

        aesea     = [ (emis_sea_harmonics.evn + emis_sea_harmonics.ev0 + emis_sea_harmonics.ev1 .* cosd(rowd) + emis_sea_harmonics.ev2 .* cosd(2*rowd) )' (emis_sea_harmonics.ehn + emis_sea_harmonics.eh0 + emis_sea_harmonics.eh1 .* cosd(rowd) + emis_sea_harmonics.eh2 .* cosd(2*rowd) )' ]'; 

      end

    end

    arefl_sea = z_sea(1:2,:) .* ( 1 - aesea );

    aesea     = aesea(:);
    arefl_sea = arefl_sea(:);


    %=== corresponding bts

    if nza == 1
      % common angle all freqs
      SST      = surf_input.SST * ones( 2*nch*nza, 1 );
    else
      % individual angle per freqs
      SST      = surf_input.SST * ones( 2*nza, 1 );
    end
    bts_sea  = atbdown .* atra .* arefl_sea + SST .* atra .* aesea;

  else

    %= with harmonics format  
    %  nza x (4 (v-h-t3-t4) x nfreqs ) 

    %  fOR toa v-POL
    %
    %  th = toa_tot_iso + th_toa_c1*c1 + th_toa_c2*c2 = tba_up  + tau_up * zv* (1 - (eh0 + eh_0 + eh_c1 *  c1  + eh_c2 * c2) ) *
    %  * tba_dn   +    tau_up   * sst * ( eh0    +   eh_h0  +  eh_c1 *  c1  + eh_c2 * c2 ) = 
    %  = tba_up      +     tau_up * zv * (1 - eh0 ) * tba_dn     -     tau_up   * zv * eh_h0 * tba_dn      +   
    %  +  tau_up  * sst * ( eh0    +   eh_h0  +   eh_c1 *  c1  + eh_c2 * c2 )
    %
    %  Decomposition:
    %
    %  th_toa_spec =  tau_up * zv* (1 - eh0 ) * tba_dn       +     tau_up  * sst * eh0  +  tba_up
    %  th_toa_eh0  =  tau_up * eh_h0 * (sst - zv* tba_dn )
    %  th_toa_c1   =  tau_up * eh_c1 * (sst - zv* tba_dn )
    %  th_toa_c2   =  tau_up * eh_c2 * (sst -zv* tba_dn )
    %
    %  Total:
    %
    %  th = th_toa_spec   +    th_toa_eh0   +  th_toa_c1 *c1     +  th_toa_c2  * c2
    % 
    %  Fot 3RD STOKES

    %  tU =    tU_toa_s1 * s1   +    tU_toa_s2 * s2     =       - tba_up *  (eU_s1 * s1 + eU_s2 * s2) )  + 
    %  + tau_up   * sst *  (eU_s1 * s1)  +  tau_up  * sst * (eU_s2 * s2) = 
    %
    %  tU_toa_s1  =  tau_up   *eU_s1  * (sst    -   tba_dn)
    %  tU_toa_s2  =  tau_up   *eU_s2  * ( sst   -   tba_dn)


    SST = surf_input.SST * ones( nch, nza );  

    bts_sea_harm_iso = zeros(4,nch,nza);
    bts_sea_harm_cs1 = zeros(4,nch,nza);
    bts_sea_harm_cs2 = zeros(4,nch,nza);

    esea1 = reshape( emis_sea_harmonics.ev0 + emis_sea_harmonics.evn, nch, nza );
    esea2 = reshape( emis_sea_harmonics.eh0 + emis_sea_harmonics.ehn, nch, nza );

    zsea1 = reshape(z_sea(1,:),nch,nza);
    zsea2 = reshape(z_sea(2,:),nch,nza);


    %= isotropic component wind direction independent
    %  including specular and oiughness

    bts_sea_harm_iso(1,:,:)  = (tbdown .* tra .* zsea1 .* ( 1 - esea1 ) ) + SST .* tra .* esea1;
    bts_sea_harm_iso(2,:,:)  = (tbdown .* tra .* zsea2 .* ( 1 - esea2 ) ) + SST .* tra .* esea2;

    %= anysotropic components to be used inside the simulator
    %  with the sin and cos of relative wind direction

    bts_sea_harm_cs1(1,:,:) = tra .* reshape(emis_sea_harmonics.ev1,nch,nza) .*  ( SST - zsea1 .* tbdown );
    bts_sea_harm_cs1(2,:,:) = tra .* reshape(emis_sea_harmonics.eh1,nch,nza) .*  ( SST - zsea2 .* tbdown );
    bts_sea_harm_cs1(3,:,:) = tra .* reshape(emis_sea_harmonics.e31,nch,nza) .*  ( SST - tbdown );
    bts_sea_harm_cs1(4,:,:) = tra .* reshape(emis_sea_harmonics.e41,nch,nza) .*  ( SST - tbdown );

    bts_sea_harm_cs2(1,:,:) = tra .* reshape(emis_sea_harmonics.ev2,nch,nza) .*  ( SST - zsea1 .* tbdown );
    bts_sea_harm_cs2(2,:,:) = tra .* reshape(emis_sea_harmonics.eh2,nch,nza) .*  ( SST - zsea2 .* tbdown );
    bts_sea_harm_cs2(3,:,:) = tra .* reshape(emis_sea_harmonics.e32,nch,nza) .*  ( SST - tbdown );
    bts_sea_harm_cs2(4,:,:) = tra .* reshape(emis_sea_harmonics.e42,nch,nza) .*  ( SST - tbdown );

    %= emis

    esea  = nan(2,nch,nza);
    esea(1,:,:) = esea1;
    esea(2,:,:) = esea2;

  end

else

  if ~do_harm

    if nza == 1
      % common angle all freqs
      bts_sea = zeros( nch * nza * 2, 1 );
    else
      % individual angle per freq
      bts_sea = zeros( nza * 2, 1 );
    end

  else

    esea            = nan( 2, nch, nza );
    bts_sea_harm_iso = zeros(4,nch,nza);
    bts_sea_harm_cs1 = zeros(4,nch,nza);
    bts_sea_harm_cs2 = zeros(4,nch,nza);

  end

end


%============================= sea ice surface



if pc_ice > 0

  do_ice = 1;

  % ice surface, sea+ice

  if ~isfield( surf_input, 'SIEMIS' )

    %===  No external emis so calculations 

    if isfield( surf_input, 'SITH' )

      %= parameterized sea-ice

      ifreq = sensor_input.F;
      ifreq( ifreq == 1.4135 ) = 1.41;

      %= CIMR L-band is at 1.4135, we round the freq to the SMAP one 
      %  to use the derived emis using SMAP obs as the NN selection
      %  searchs for a 1.41 GHz freq in the code

      % original instantaneous fields
      %
      %eice = fom_seaice_emis_surfem(  ifreq, sensor_input.ZA, [], [], [], surf_input.MONTH, surf_input.DAY, surf_input.SITH, surf_input.SISNTH,surf_input.SIA, surf_input.SIR, surf_input.SIT, TCWV);   

      % SIT and TCWV daily fields
      %
      eice = fom_seaice_emis_surfem(  ifreq, sensor_input.ZA, [], [], [], surf_input.MONTH, surf_input.DAY, surf_input.SITH, surf_input.SISNTH,surf_input.SIA, surf_input.SIR, surf_input.DSIT, surf_input.DTCWV);   
 
    else

      if ~isfield( surf_input, 'SIMONTH' )

        ifreq = sensor_input.F;
        ifreq( ifreq == 1.4135 ) = 1.41;

        path_rsrc = [ dir_input, '/EmisSeaIce' ];
        eice = fom_seaice_emis_surfem(  ifreq, sensor_input.ZA, surf_input.LAT, surf_input.LON, path_rsrc, surf_input.MONTH, surf_input.DAY);

      else

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

        eice = emis_land_angle_interpolate( surf_input.LON, surf_input.LAT, surf_input.SIMONTH, dir_input, sensor_input.F, sensor_input.ZA  ); 

      end
  
    end  

  else

    %=== using passed value
   
    eice = surf_input.SIEMIS; 

  end 

  %=== corresponding bts


  if ~do_harm


    aeice = reshape( eice, 2, nch *nza);
    aeice = aeice(:); 

    ts       = surf_input.SIT *ones( 2*nch*nza, 1 );  
    bts_ice  = atbdown .* atra .* ( 1-aeice ) + ts .* atra .* aeice;

  else

    ts       = surf_input.SIT *ones( nch, nza );  

    bts_ice_harm_iso  = zeros( 4, nch, nza );
    bts_ice_harm_iso(1,:,:) = tbdown .* tra .* ( 1 - squeeze(eice(1,:,:)))  + ts .* tra .* squeeze(eice(1,:,:));
    bts_ice_harm_iso(2,:,:) = tbdown .* tra .* ( 1 - squeeze(eice(2,:,:)))  + ts .* tra .* squeeze(eice(2,:,:));

    bts_ice_harm_cs1 = zeros( 4, nch, nza );
    bts_ice_harm_cs2 = zeros( 4, nch, nza );

  end


else

  if ~do_harm

    if nza == 1
      % common angle all freqs
      bts_ice = zeros( nch * nza * 2, 1 );
    else
      % individual angle per freq
      bts_ice = zeros( nza * 2, 1 );
    end

  else

    eice            = nan( 2, nch, nza );
    bts_ice_harm_iso = zeros( 4, nch, nza );
    bts_ice_harm_cs1 = zeros( 4, nch, nza );
    bts_ice_harm_cs2 = zeros( 4, nch, nza );


  end

end   



%===== TB top of atmosphere Vpol and Hpol for  mixed scenes
%      corresponding to the selected incidence angle 
%      and azimuthal angle in the control file

if ~do_harm

  %= outputing total bts V and H-pol
  %  used for retrievals where we did no deal
  %  with t3 and t4 for the moment

  bts          = atbup +  pc_sea * bts_sea + pc_lan * bts_lan + pc_ice * bts_ice;
  if do_pol
    bts = bts(indpol);
  end

else

  %= for RT simulations outputting the harmonics

  bts = [];

  %= TOA harmonic TB components
  %            this outputs the TOA contributions for the ocean flat and
  %            rough and land/sea-ice emission (iso), the first harmonic (cs1) an
  %            the second (cs2) using expected input to OSS


  bts_harm_iso = atbup +  pc_sea * bts_sea_harm_iso + pc_lan * bts_lan_harm_iso + pc_ice * bts_ice_harm_iso;
  bts_harm_cs1 = pc_sea * bts_sea_harm_cs1 + pc_lan * bts_lan_harm_cs1 + pc_ice * bts_ice_harm_cs1;
  bts_harm_cs2 = pc_sea * bts_sea_harm_cs2 + pc_lan * bts_lan_harm_cs2 + pc_ice * bts_ice_harm_cs2;

end


%=== removing cosmic background from tbdown to be output
%    tbdown	= tb1 + tra .* TBC;

tbdown = tbdown - ( tra .* tbc );


return





%==========================================================================
%  SCEPS-I land climato emis before we have the SCEPS-II parameterizarion


function emis  = emis_land_angle_interpolate( lon, lat, month, dir_input, freqs, angles )


nf = length( freqs );
na = length( angles );

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

e_clim  = nan(nf,2);

%= land emis per frequency  
[e_clim(1,:), class1 ] = emis_land_surface_table( lon, lat, month, dir_input, floor(freqs(1))); 
for f = 2:nf
  e_clim(f,:) = emis_land_surface_table( lon, lat, month, dir_input, floor(freqs(f))); 
end 



%= better have a made up emis than a -999
%  value for the few cases where no emis
%  lan exists

e_clim( e_clim(:,1) == -999, 1 ) = 0.95;
e_clim( e_clim(:,2) == -999, 2 ) = 0.90;


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
    
    theta  = angles(a);

    %=  Calculation of the emis at theta=0° with a multilinear regression
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








