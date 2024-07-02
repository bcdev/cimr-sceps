%=================================================================
%
% Mscript to faciliate the use of SURFEM-ocean, hard-coding the 
% MLPs weights and biases and providing the FASTEM correction
% for the reflectivity to account for the multiple reflection of
% the downwelling atmospheric radiation. It can be run as:
%
% (1) Emis calculation:
%    [e ] = fom_emis_sea_surface( sst, sss, ows, owd, freq, angle) 
%
% (2) Emis and reflectivity derivatives calculation:
%    [e, r ] = fom_emis_sea_surface( sst, sss, ows, owd, freq, angle) 
%
% (3) Emis, reflectivity, and  derivatives calculation:
%    [e, r, e_dt,e_ds,e_dw ] = fom_emis_sea_surface( sst, sss, ows, owd, freq, angle, trans) 


%
% FORMAT   [e, r, e_dt,e_ds,e_dw ] = 
%			fom_emis_sea_surface( sst, sss, ows, owd, freq, angle, trans)
%        
% OUT   e [-]			Sea surface emissivity, V-pol, H-pol, S3, and
%				S4 components
%       r [-]                   Sea surface reflectivity,  V-pol, H-pol, S3, and
%				S4 components
%
%	e_dt [-]		Derivative De/DSST
%	e_ds [-]		Derivative De/DSSS
%	e_dw [-]		Derivative De/DOWS	       		
%
% IN    sst [K]			Sea surface temperature
%       sss [psu]		Sea surface salinity 
%       ows [m/s]		Ocean wind speed
%	owd [degrees]		Ocean wind direction
%	freq [GHz]		Observing frequency
%	angle [degrees]		Zenith viewing angle
%	trans [-]		Atmospheric transmission 

function [e, r, e_dt, e_ds, e_dw ] = fom_emis_sea_surface( sst, sss, ows, owd, freq, angle, trans)


if nargout >= 2 & nargin < 7
 error('Calculation of r demanded but no trans given')
end


%=== calculate emis and jacobians

if nargout <= 2
  [ e ] = SURFEM_ocean_K( freq, angle, sst, sss, ows, owd);
else
  [ e, e_dt, e_ds, e_dw] = SURFEM_ocean_K( freq, angle, sst, sss, ows, owd);
end


%=== calculate reflectivity correcting with FASTEM 6 routine

if nargin == 7

  [z_v,z_h] = zreflmodv_tan( freq, trans, angle, ows);
  r = [ z_v .* (1-e(1)); z_h .* (1-e(2)); zeros(size(z_v)); zeros(size(z_v)) ];

end



return


%------------------------------------------------------------------------
%
%%SURface Fast Emissivity Model for Ocean (SURFEM-Ocean) with emissivity Jacobians
%developped by Lise Kilic, Catherine Prigent and Carlos Jimenez
% Code HISTORY: 20220520  Created by lise.kilic(at)obspm.fr
% HISTORY: 20220520  Created by lise.kilic(at)obspm.fr

%----INPUTS----
%freq_vec : frequency in GHz
%theta_vec : incidence angle in degrees
%SST_vec : Sea Surface Temperature in K
%SSS_vec : Sea Surface Salinity in psu
%OWS_vec : Ocean Wind Speed in m/s
%phi_vec : relative wind direction in degrees
%----OUTPUTS----
%e=[ev; eh; e3; e4]

%--
%The size of the inputs is 1 x number of cases
%The size of the outputs is 4 x number of cases (4 is the number of
%polarizations (V, H, S3, S4))
%
%------------------------------------------------------------------------


function [e,e_dt,e_ds,e_dw] = SURFEM_ocean_K(freq_vec,theta_vec,SST_vec, SSS_vec, OWS_vec,phi_vec)

if nargout == 4
  do_k = 1;
else
  do_k = 0;
end


%-declare variables----
len=max(size(freq_vec));
freq_vec=reshape(freq_vec,len,1);
theta_vec=reshape(theta_vec,len,1);
SST_vec=reshape(SST_vec,len,1);
SSS_vec=reshape(SSS_vec,len,1);
OWS_vec=reshape(OWS_vec,len,1);
phi_vec=reshape(phi_vec,len,1);

evn=zeros(1,len);
ehn=zeros(1,len);

if do_k
  evn_dt=zeros(1,len);
  ehn_dt=zeros(1,len);
  evn_ds=zeros(1,len);
  ehn_ds=zeros(1,len);
end


%%---Flat sea calculation (neutral wind)--------


for ii=1:len
    
  freq=freq_vec(ii); 
  SST=SST_vec(ii);
  SSS=SSS_vec(ii);
  theta=theta_vec(ii);


  fhz=freq*1e9; % frequency in Hz
  SSTc=SST-273.15; % Temperature in °C

  %-Estimation e neutral (flat sea no wind)-------
  %Dielectric constant module
  [epsi,epsi_dt,epsi_ds]= epsilon_MW_K(SSTc, SSS, fhz);
  % !       Inputs:
  % !               SST : Sea surface temperature in degrees C
  % !               SSS  : Sea surface salinity in psu
  % !               f    : frequency in Hz
  % ! Convention imag(epsi)>0
                
  cothet = cosd(theta);
  sithet = sind(theta);
  epsibar=conj(epsi);

  %fresnel reflection coefficient
  rv = (epsi*cothet - sqrt(epsi-sithet*sithet))/(epsi*cothet + sqrt(epsi-sithet*sithet));
  rh = (cothet - sqrt(epsi-sithet*sithet))/(cothet + sqrt(epsi-sithet*sithet));
  rvbar = (epsibar*cothet - sqrt(epsibar-sithet*sithet))/(epsibar*cothet + sqrt(epsibar-sithet*sithet));
  rhbar = (cothet - sqrt(epsibar-sithet*sithet))/(cothet + sqrt(epsibar-sithet*sithet));

  Rvv0=real(rv).^2 +imag(rv).^2;
  Rhh0=real(rh).^2 +imag(rh).^2;

  if do_k

    %derivative of R
    %v
    rv_re_depsi=(cothet*(epsi-2*sithet*sithet))/(sqrt(epsi-sithet*sithet)*(sqrt(epsi-sithet*sithet)+epsi*cothet)^2);
    rv_im_depsi=(1i*cothet*(epsi-2*sithet*sithet))/(sqrt(epsi-sithet*sithet)*(sqrt(epsi-sithet*sithet)+epsi*cothet)^2);
    %conjugate
    rv_re_depsibar=(cothet*(epsibar-2*sithet*sithet))/(sqrt(epsibar-sithet*sithet)*(sqrt(epsibar-sithet*sithet)+epsibar*cothet)^2);
    rv_im_depsibar=(-1i*cothet*(epsibar-2*sithet*sithet))/(sqrt(epsibar-sithet*sithet)*(sqrt(epsibar-sithet*sithet)+epsibar*cothet)^2);
  %h
   rh_re_depsi=(-cothet)/(sqrt(epsi-sithet*sithet)*(sqrt(epsi-sithet*sithet)+cothet)^2);
   rh_im_depsi=(-1i*cothet)/(sqrt(epsi-sithet*sithet)*(sqrt(epsi-sithet*sithet)+cothet)^2);
   %conjugate
   rh_re_depsibar=(-cothet)/(sqrt(epsibar-sithet*sithet)*(sqrt(epsibar-sithet*sithet)+cothet)^2);
   rh_im_depsibar=(1i*cothet)/(sqrt(epsibar-sithet*sithet)*(sqrt(epsibar-sithet*sithet)+cothet)^2);

   Rvv0_depsi=(rv_re_depsi*rvbar+rv_re_depsibar*rv) -1i*((rv_im_depsi*rvbar+rv_im_depsibar*rv));
   Rhh0_depsi=(rh_re_depsi*rhbar+rh_re_depsibar*rh) -1i*((rh_im_depsi*rhbar+rh_im_depsibar*rh));

   %derivative of emissivity

   evn_dt(ii)=-Rvv0_depsi*epsi_dt;
   ehn_dt(ii)=-Rhh0_depsi*epsi_dt;

   evn_ds(ii)=-Rvv0_depsi*epsi_ds;
   ehn_ds(ii)=-Rhh0_depsi*epsi_ds;  

  end

  evn(ii)=1-Rvv0;
  ehn(ii)=1-Rhh0;

end

if do_k

  evn_dt=real(evn_dt);
  ehn_dt=real(ehn_dt);
  evn_ds=real(evn_ds);
  ehn_ds=real(ehn_ds); 
 %%---end of flat sea calculations------

end

%%---Rough Sea Calculations------------

inputs=[freq_vec,theta_vec,OWS_vec,SST_vec,SSS_vec,evn',ehn']';

%%--Estimation ev0 and eh0 isotropic signal due to wind and Jacobians
%net = load('./MatNNData/e_isotropic_net_data.mat');
[ W1, b1, W2, b2, ti, to ]   = nn_isotropic;

if do_k
  %[ outputs, J_wind0 ] = nn_analytical_jacobians( inputs, net.W1, net.b1, net.W2, net.b2, net.ti, net.to, net.fcn );
  [ outputs, J_wind0 ] = nn_analytical_jacobians( inputs, W1, b1, W2, b2, ti, to );
else
  [ outputs, J_wind0 ] = nn_analytical_jacobians( inputs, W1, b1, W2, b2, ti, to );
end

ev0=outputs(1,:);
eh0=outputs(2,:);
clear outputs


%%--Estimation of ev1,eh1,eS31,eS41,ev2,eh2,eS32,eS42----
%-for anisotropic wind (windir) and Jacobians
%net = load('./MatNNData/e_anisotropic_net_data.mat');
[ W1, b1, W2, b2, ti, to ]   = nn_anisotropic;
if do_k
  %[ outputs, J_windir ] = nn_analytical_jacobians( inputs, net.W1, net.b1, net.W2, net.b2, net.ti, net.to, net.fcn );
  [ outputs, J_windir ] = nn_analytical_jacobians( inputs, W1, b1, W2, b2, ti, to );
else
   outputs = nn_analytical_jacobians( inputs, W1, b1, W2, b2, ti, to );
end


ev1=outputs(1,:);
eh1=outputs(2,:);
e31=outputs(3,:);
e41=outputs(4,:);
ev2=outputs(5,:);
eh2=outputs(6,:);
e32=outputs(7,:);
e42=outputs(8,:);
clear net outputs


%%---Estimation on the total surface emissivity----

ev=evn + ev0 + ev1.*cosd(phi_vec')+ev2.*cosd(2.*phi_vec');
eh=ehn + eh0 + eh1.*cosd(phi_vec')+eh2.*cosd(2.*phi_vec');
e3=e31.*sind(phi_vec')+e32.*sind(2.*phi_vec');
e4=e41.*sind(phi_vec')+e42.*sind(2.*phi_vec');


if do_k

  %%---Estimation of the jacobians of the emissivity---
  %(de/dt) t=temperature
  ev0_dt=squeeze(J_wind0(:,1,4))';
  eh0_dt=squeeze(J_wind0(:,2,4))';
  ev0_devn=squeeze(J_wind0(:,1,6))';
  ev0_dehn=squeeze(J_wind0(:,1,7))';
  eh0_devn=squeeze(J_wind0(:,2,6))';
  eh0_dehn=squeeze(J_wind0(:,2,7))';

  ev1_dt=squeeze(J_windir(:,1,4))';
  eh1_dt=squeeze(J_windir(:,2,4))';
  e31_dt=squeeze(J_windir(:,3,4))';
  e41_dt=squeeze(J_windir(:,4,4))';
  ev1_devn=squeeze(J_windir(:,1,6))';
  ev1_dehn=squeeze(J_windir(:,1,7))';
  eh1_devn=squeeze(J_windir(:,2,6))';
  eh1_dehn=squeeze(J_windir(:,2,7))';
  e31_devn=squeeze(J_windir(:,3,6))';
  e31_dehn=squeeze(J_windir(:,3,7))';
  e41_devn=squeeze(J_windir(:,4,6))';
  e41_dehn=squeeze(J_windir(:,4,7))';

  ev2_dt=squeeze(J_windir(:,5,4))';
  eh2_dt=squeeze(J_windir(:,6,4))';
  e32_dt=squeeze(J_windir(:,7,4))';
  e42_dt=squeeze(J_windir(:,8,4))';
  ev2_devn=squeeze(J_windir(:,5,6))';
  ev2_dehn=squeeze(J_windir(:,5,7))';
  eh2_devn=squeeze(J_windir(:,6,6))';
  eh2_dehn=squeeze(J_windir(:,6,7))';
  e32_devn=squeeze(J_windir(:,7,6))';
  e32_dehn=squeeze(J_windir(:,7,7))';
  e42_devn=squeeze(J_windir(:,8,6))';
  e42_dehn=squeeze(J_windir(:,8,7))';

  ev_dt=evn_dt + ev0_dt + ev0_devn.*evn_dt + ev0_dehn.*ehn_dt + cosd(phi_vec').*(ev1_dt + ev1_devn.*evn_dt +ev1_dehn.*ehn_dt) + cosd(2.*phi_vec').*(ev2_dt + ev2_devn.*evn_dt +ev2_dehn.*ehn_dt);
  eh_dt=ehn_dt + eh0_dt + eh0_devn.*evn_dt + eh0_dehn.*ehn_dt+ cosd(phi_vec').*(eh1_dt + eh1_devn.*evn_dt +eh1_dehn.*ehn_dt) + cosd(2.*phi_vec').*(eh2_dt + eh2_devn.*evn_dt +eh2_dehn.*ehn_dt);
  e3_dt=sind(phi_vec').*(e31_dt + e31_devn.*evn_dt +e31_dehn.*ehn_dt) +sind(2.*phi_vec').*(e32_dt + e32_devn.*evn_dt +e32_dehn.*ehn_dt);
  e4_dt=sind(phi_vec').*(e41_dt + e41_devn.*evn_dt +e41_dehn.*ehn_dt) +sind(2.*phi_vec').*(e42_dt + e42_devn.*evn_dt +e42_dehn.*ehn_dt);

  %de/ds s=salinity
  ev0_ds=squeeze(J_wind0(:,1,5))';
  eh0_ds=squeeze(J_wind0(:,2,5))';

  ev1_ds=squeeze(J_windir(:,1,5))';
  eh1_ds=squeeze(J_windir(:,2,5))';
  e31_ds=squeeze(J_windir(:,3,5))';
  e41_ds=squeeze(J_windir(:,4,5))';

  ev2_ds=squeeze(J_windir(:,5,5))';
  eh2_ds=squeeze(J_windir(:,6,5))';
  e32_ds=squeeze(J_windir(:,7,5))';
  e42_ds=squeeze(J_windir(:,8,5))';


  ev_ds=evn_ds + ev0_ds + ev0_devn.*evn_ds + ev0_dehn.*ehn_ds + cosd(phi_vec').*(ev1_ds + ev1_devn.*evn_ds +ev1_dehn.*ehn_ds) + cosd(2.*phi_vec').*(ev2_ds + ev2_devn.*evn_ds +ev2_dehn.*ehn_ds);
  eh_ds=ehn_ds + eh0_ds + eh0_devn.*evn_ds + eh0_dehn.*ehn_ds+ cosd(phi_vec').*(eh1_ds + eh1_devn.*evn_ds +eh1_dehn.*ehn_ds) + cosd(2.*phi_vec').*(eh2_ds + eh2_devn.*evn_ds +eh2_dehn.*ehn_ds);
  e3_ds=sind(phi_vec').*(e31_ds + e31_devn.*evn_ds +e31_dehn.*ehn_ds) +sind(2.*phi_vec').*(e32_ds + e32_devn.*evn_ds +e32_dehn.*ehn_ds);
  e4_ds=sind(phi_vec').*(e41_ds + e41_devn.*evn_ds +e41_dehn.*ehn_ds) +sind(2.*phi_vec').*(e42_ds + e42_devn.*evn_ds +e42_dehn.*ehn_ds);

  %de/dw w=wind speed
  ev_dw=squeeze(J_wind0(:,1,3))' +cosd(phi_vec').*squeeze(J_windir(:,1,3))' + cosd(2.*phi_vec').*squeeze(J_windir(:,5,3))';
  eh_dw=squeeze(J_wind0(:,2,3))' +cosd(phi_vec').*squeeze(J_windir(:,2,3))' + cosd(2.*phi_vec').*squeeze(J_windir(:,6,3))';
  e3_dw=sind(phi_vec').*squeeze(J_windir(:,3,3))'+sind(2.*phi_vec').*squeeze(J_windir(:,7,3))';
  e4_dw=sind(phi_vec').*squeeze(J_windir(:,4,3))'+sind(2.*phi_vec').*squeeze(J_windir(:,8,3))';

  %--condition for Wind Speed=0
  indn=find(OWS_vec==0);
  ev(indn)=evn(indn);
  eh(indn)=ehn(indn);
  e3(indn)=0;
  e4(indn)=0;
  ev_dt(indn)=real(evn_dt(indn));
  eh_dt(indn)=real(ehn_dt(indn));
  e3_dt(indn)=0;
  e4_dt(indn)=0;
  ev_ds(indn)=real(evn_ds(indn));
  eh_ds(indn)=real(ehn_ds(indn));
  e3_ds(indn)=0;
  e4_ds(indn)=0;

end

%---Final emissivity vectors and jacobians-----
e=[ev;eh;e3;e4];

if do_k
  e_dt=[ev_dt;eh_dt;e3_dt;e4_dt];
  e_ds=[ev_ds;eh_ds;e3_ds;e4_ds];
  e_dw=[ev_dw;eh_dw;e3_dw;e4_dw];
end

return

%------------------------------------------------------------------------
%
% ! Computes the dielectric constant of sea water according to the model
% ! by Meissner and Wentz (2004) including corrections by Meissner and Wentz
% ! 2012 and Meissner et al. 2014, with the corresponding jacobians (K)
% ! Translated from fortran to matlab by L. Kilic (2021) 
% ! Computation of the jacobians (depsi/dt and depsi/ds) by L. Kilic (2022)
% !
% !       Inputs:
% !               SST : Sea surface temperature in degrees C
% !               SSS  : Sea surface salinity in psu
% !               f    : frequency in Hz
% !
% !       Outputs:
% !               epsi_MW : complex dielectric constant of sea water
% !               epsi_dt : derivative according to SST
% !               epsi_ds : derivative according to SSS
% !
% !
% !
% !       E. Dinnat
% 
% ! Meissner, T., & Wentz, F. J. (2004). The complex dielectric constant of pure
% ! and sea water from microwave satellite observations. IEEE Transactions on
% ! Geoscience and Remote Sensing, 42(9), 18361849.
% ! https://doi.org/10.1109/TGRS.2004.831888
% !
% ! Meissner, T., & Wentz, F. J. (2012). The Emissivity of the Ocean Surface
% ! Between 6 and 90 GHz Over a Large Range of Wind Speeds and Earth Incidence
% ! Angles. IEEE Transactions on Geoscience and Remote Sensing, 50(8), 30043026.
% ! https://doi.org/10.1109/TGRS.2011.2179662
% !
% ! Meissner, T., Wentz, F. J., & Ricciardulli, L. (2014). The emission and
% ! scattering of L-band microwave radiation from rough ocean surfaces and wind
% ! speed measurements from the Aquarius sensor. Journal of Geophysical Research C:
% ! Oceans, 119(9), 64996522. https://doi.org/10.1002/2014JC009837
%
%------------------------------------------------------------------------

function [epsi_MW,epsi_dt,epsi_ds] = epsilon_MW_K(SST,SSS,f)


%       ! Intermediate variables
% !     ai coefficients from Table III in Meissner and Wentz 2004
ai=[ 5.7230D+00, 2.2379D-02, -7.1237D-04, 5.0478D+00, -7.0315D-02, 6.0059D-04, 3.6143D+00, 2.8841D-02, 1.3652D-01,  1.4825D-03,2.4166D-04 ];
%!     bi coefficients from Table VI in Meissner and Wentz 2004
bi=[-3.56417D-03,  4.74868D-06, 1.15574D-05, 2.39357D-03, -3.13530D-05, 2.52477D-07, -6.28908D-03, 1.76032D-04, -9.22144D-05, -1.99723D-02,1.81176D-04, -2.04265D-03,  1.57883D-04];
% !    updated bi (i = 0 -> 2) coefficients from Table VI in Meissner and Wentz
% !    2012
bi_new1=[ -0.33330D-02,  4.74868D-06, 0.0D0];
% !    updated bi (i = 3 -> 5) coefficients from Table VII in Meissner and Wentz
% !    2012 (now uses 5 coefficients instead of 3), includes change of
% !    sign on coef#4 reported in Meissner et al. 2016
bi_new2=[0.23232D-02, -0.79208D-04, 0.36764D-05, -0.35594D-06, 0.89795D-08];
%       REAL*8 :: fGHz, f0, SST2, SST3, SST4, SSS2
%       REAL*8 :: es, e1, einf, nu1,nu2 ! Fresh water parameters
%       REAL*8 :: sigma35, R15, alpha0, alpha1, RTR15, sigma
%       REAL*8 :: es_s, e1_s, einf_s, nu1_s, nu2_s ! Salt water parameters
%       REAL*8 :: c1, c2, c3, c4 ! conversion coef from fresh to salt water
%       COMPLEX*16, PARAMETER :: j = (0.D0, 1.D0)

      if SST<-40
        disp('Error in epsilon_MW. SST is lower than -40 degree C.')
        return
      end
      
      if SST>34
        disp('Error in espilon_MW. SST is larger than +34 degree C.')
        return 
      end
      
      if SSS>40
        disp(' Error in espilon_MW. SSS is larger than +40 psu.')
        return
      end
      
      if (SSS<0)
        disp(' Error in espilon_MW. SSS is negative.')
        return
      end
        

      fGHz = f/1.D09 ;%! convert fequency from Hz to Ghz
      f0 = 17.97510D0; %! (1/(2*pi*epsi_0) term in conductivity term in GHz.m/S

      SST2 = SST * SST;
      SST3 = SST2 * SST;
      SST4 = SST3 * SST;

      SSS2 = SSS * SSS;
      SSS3 = SSS2 * SSS;
%  
%       ! Fresh water parameters
%       ! es:  Static dielectric constant for pure water by Stogryn et al. 1995
      gi=[3.70886D4, -8.2168D1, 4.21854D2];

      es = ( gi(1) + gi(2)*SST ) / ( gi(3) + SST );
      e1 = ai(1) + ai(2)*SST + ai(3)*SST2;
      nu1 = ( 45.00D0 + SST ) / ( ai(4) + ai(5)*SST + ai(6)*SST2 );
      einf = ai(7) + ai(8)*SST;
      nu2 = ( 45.00D0 + SST ) / ( ai(9) + ai(10)*SST + ai(11)*SST2 );
 
%       ! Salt water
%       ! Conductivity [s/m] by Stogryn et al. 1995

      ci=[2.903602D0, 8.60700D-2, 4.738817D-4, -2.9910D-6, 4.3047D-9, 37.5109D0, 5.45216D0, 1.4409D-2, 1004.75D0, 182.283D0, 6.9431D0, 3.2841D0, -9.9486D-2, 84.850D0, 69.024D0, 49.843D0, -0.2276D0, 0.198D-2];
      sigma35 = ci(1)+ ci(2)*SST + ci(3)*SST2 + ci(4)*SST3 + ci(5)*SST4;
      R15 = SSS*( ci(6) + ci(7)*SSS + ci(8)*SSS2 ) / ( ci(9) + ci(10)*SSS + SSS2 );
      alpha0 = ( ci(11) + ci(12)*SSS + ci(13)*SSS2 ) / (ci(14) + ci(15)*SSS + SSS2 );
      alpha1 = ci(16) + ci(17)*SSS + ci(18)*SSS2;
      RTR15 = 1.0D0 + ( SST - 15.0D0 )*alpha0/ ( alpha1 + SST ) ;
      sigma = sigma35 * R15 * RTR15;
      
      
      %Derivative of sigma
      %dt
      sigma35_dt=ci(2) + 2*ci(3)*SST + 3*ci(4)*SST2 + 4*ci(5)*SST3;
      RTR15_dt=alpha0*(alpha1+15)/(alpha1+SST)^2;
      sigma_dt=R15*(sigma35_dt*RTR15+ sigma35*RTR15_dt);
      %ds
      alpha0_ds=( (ci(12)+2*ci(13)*SSS)*(ci(14)+ci(15)*SSS + SSS2) - (ci(11)+ ci(12)*SSS + ci(13)*SSS2)*(ci(15)+2*SSS))/(ci(14)+ci(15)*SSS+SSS2)^2;
      alpha1_ds=ci(17) +2*ci(18)*SSS;
      R15_ds=( (ci(6)+2*ci(7)*SSS + 3*ci(8)*SSS2)*(ci(9)+ci(10)*SSS+SSS2) - (ci(6)*SSS +ci(7)*SSS2 + ci(8)*SSS3)*(ci(10)+2*SSS))/(ci(9) + ci(10)*SSS +SSS2)^2;
      RTR15_ds= ( (SST-15)*(alpha1 +SST)*alpha0_ds + alpha0*(SST-15)*alpha1_ds )/(alpha1+SST)^2;
      sigma_ds=sigma35*(R15_ds*RTR15+R15*RTR15_ds);
      

%       ! Salt water parameters 
      c0 = exp ( bi_new1(1)*SSS + bi_new1(2)*SSS2 + bi_new1(3)*SSS*SST );
      es_s = es * c0 ;

      if ( SST < 30 )

        c1 = 1.0D0 + SSS * ( bi_new2(1) + bi_new2(2)*SST + bi_new2(3)*SST2 + bi_new2(4)*SST3 + bi_new2(5)*SST4 );
        c1_dt = SSS*(bi_new2(2)+2*bi_new2(3)*SST + 3*bi_new2(4)*SST2 + 4*bi_new2(5)*SST3);
        c1_ds = bi_new2(1) + bi_new2(2)*SST + bi_new2(3)*SST2 + bi_new2(4)*SST3 + bi_new2(5)*SST4;
      else

        c1 = 1.0D0 + SSS * ( 9.1873715D-04 + 1.5012396D-04*(SST - 30.D0 ) );
        c1_dt = SSS*1.5012396D-04;
        c1_ds = 9.1873715D-04 + 1.5012396D-04*(SST - 30.D0 );
      end

      nu1_s = nu1*c1;

      c2  = exp ( bi(7)*SSS + bi(8)*SSS2 + bi(9)*SSS*SST ) ;
      e1_s = e1*c2;

% !     c3 = 1.0D0 + SSS*(bi(10) + bi(11)*SST) ! Expression for c3  amended in Meissner et al. 2016
      c3 = 1.0D0 + SSS*(bi(10) + 0.5D0*bi(11) * (SST + 30.D0) );
      nu2_s = nu2*c3;
      c4 = 1.0D0  + SSS * ( bi(12) + bi(13) * SST );
      einf_s = einf*c4  ;
      
      %derivative of es_s nu1_s e1_s nu2_s einf_s
      %es
      es_dt=(gi(2)*gi(3)-gi(1))/(gi(3)+SST)^2;
      c0_dt=bi_new1(3)*SSS*c0;
      es_s_dt=es_dt*c0 + c0_dt*es;
      c0_ds=(bi_new1(1)+2*bi_new1(2)*SSS + bi_new1(3)*SST)*c0;
      es_s_ds=es*c0_ds;
      %nu1
      nu1_dt = ( ai(4) - 45.00D0*ai(5) - ai(6)*(90.00D0*SST +SST2 ) )/ ( ai(4) + ai(5)*SST + ai(6)*SST2 )^2;
      nu1_s_dt = nu1_dt*c1 + nu1*c1_dt;
      nu1_s_ds = c1_ds*nu1;
      %e1
      e1_dt= ai(2) + 2*ai(3)*SST;
      c2_dt=bi(9)*SSS*c2;
      e1_s_dt=e1_dt*c2 + e1*c2_dt;
      c2_ds=(bi(7) + 2*bi(8)*SSS + bi(9)*SST)*c2;
      e1_s_ds=e1*c2_ds;
      %nu2
      nu2_dt= ( ai(9) - 45.00D0*ai(10) - ai(11)*(90.00D0*SST +SST2 ) ) / ( ai(9) + ai(10)*SST + ai(11)*SST2 )^2;
      c3_dt = SSS*0.5D0*bi(11);
      nu2_s_dt=nu2_dt*c3 + nu2*c3_dt;
      c3_ds=bi(10) + 0.5D0*bi(11) * (SST + 30.D0);
      nu2_s_ds=nu2*c3_ds;
      %einf
      einf_dt = ai(8);
      c4_dt = SSS * bi(13);
      einf_s_dt=einf_dt*c4 + einf*c4_dt;
      c4_ds= bi(12) + bi(13) * SST ;
      einf_s_ds = einf*c4_ds;
     
      
      epsi_MW = ( es_s - e1_s  ) / ( 1.0D0 - 1i * ( fGHz/nu1_s ) ) + ( e1_s - einf_s ) / ( 1.0D0 - 1i * ( fGHz/nu2_s ) ) + einf_s + 1i*sigma*f0/fGHz ;
      
      %derivative of epsi_MW
      
      epsi_dt=((es_s_dt-e1_s_dt)*( 1.0D0 - 1i * ( fGHz/nu1_s ) ) - (es_s - e1_s)*1i * ( fGHz/nu1_s^2 )*nu1_s_dt)/( 1.0D0 - 1i * ( fGHz/nu1_s ))^2;
      epsi_dt=epsi_dt + ((e1_s_dt-einf_s_dt)*( 1.0D0 - 1i * ( fGHz/nu2_s ) ) - (e1_s - einf_s)*1i * ( fGHz/nu2_s^2 )*nu2_s_dt)/( 1.0D0 - 1i * ( fGHz/nu2_s ))^2;
      epsi_dt=epsi_dt + einf_s_dt + 1i*sigma_dt*f0/fGHz ;
      

      epsi_ds=((es_s_ds-e1_s_ds)*( 1.0D0 - 1i * ( fGHz/nu1_s ) ) - (es_s - e1_s)*1i * ( fGHz/nu1_s^2 )*nu1_s_ds)/( 1.0D0 - 1i * ( fGHz/nu1_s ))^2;
      epsi_ds=epsi_ds + ((e1_s_ds-einf_s_ds)*( 1.0D0 - 1i * ( fGHz/nu2_s ) ) - (e1_s - einf_s)*1i * ( fGHz/nu2_s^2 )*nu2_s_ds)/( 1.0D0 - 1i * ( fGHz/nu2_s ))^2;
      epsi_ds=epsi_ds + einf_s_ds + 1i*sigma_ds*f0/fGHz ;



return



%------------------------------------------------------------------------
%
% NAME:    [ Y, J ] = nn_analytical_jacobians( X, W1, b1, W2, b2, ti, to)
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
%
% OUT:	Y	net output for X, matrix [no * ns]
%       J       net Jacobian for X, matrix [ns * no * ni]
%
% IN:	X       net input, matrix [ni * ns]
%	W1	weight matrix of first layer, matrix [nh * ni]
%	b1	bias vector of first layer, vector [nh * 1]
%	W2	weight matrix of second layer, matrix [no * nh]
%	b2	bias vector of first layer, vector [no * 1]
%       ti	structure having the settings of the input
%		mapminmax transformation
%       to	structure having the settings of the output
%		mapminmax transformation
%	fcn	strucure conating the net activation and transformation
%		functions.	
% 
%------------------------------------------------------------------------

% HISTORY: 20220308  Created by carlos.jimenez(at)estellus.fr

function [ Y, J ] = nn_analytical_jacobians( X, W1, b1, W2, b2, ti, to )


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



if nargout > 1

%=== analytical Jacobian
%    chain of rule
%
%    z(x) = ((x - oymi ) ./ ogai ) + ooff 
%    f(x) = x
%    k(x) = W2 * x + b2;
%    g(x) = tansig(x)
%    h(x) = W1 * x + b1;
%    u(x) = iymi + (x - ioff).* igai;
%
%    dz/dx = dz/df * df/dk * dk/dg * dg/dh * dh/du * du/dx


% dz/df
dz_df   = (1./ogai(:,1));
dz_df   = repmat( dz_df, 1, nh); 
 
% df/dk
df_dk = 1;

% dk/dg
dk_dg = W2;

% dh/du
dh_du = W1;

% du/dx 
du_dx   = igai(:,1);
du_dx   = repmat( du_dx, 1, no)'; 


J = nan( nx, no, ni );


for n = 1:nx

  % dg/dh
  dg_dh = ( 4* exp(-2*x1(:,n)) ) ./ ( 1 + exp(-2*x1(:,n)) ).^ 2; 
  dg_dh = repmat( dg_dh,1,ni);

  %    dz/dx = dz/df * df/dk * dk/dg * dg/dh * dh/du * du/dx
  J(n,:,:) =  ( dz_df * df_dk .* dk_dg) * (dg_dh .* dh_du) .* du_dx;

end

end


return




%------------------------------------------------------------------------
%
% For the restitution of the sea surface reflectivity with SURFEM that we use the same method than with FASTEM.
%
% Meaning that at the end of the code we have:
%
%    Reflectivity(1)  = zreflmod_v * (ONE-Emissivity(1))
%    Reflectivity(2)  = zreflmod_h * (ONE-Emissivity(2))
%    Reflectivity(3:4)  = ZERO
% and zreflmod_v and zreflmod_h are computed the same way that it is in FASTEM
%
%------------------------------------------------------------------------

function [zreflmod_v,zreflmod_h] = zreflmodv_tan(Frequency,Transmittance,Zenith_Angle,Wind_Speed)

[zreflmod_v,zreflmod_h] = zreflmodv(Frequency,Transmittance,Zenith_Angle,Wind_Speed);

N=length(Frequency);
[zinfv,zinfh] = zreflmodv(Frequency,Transmittance,56.*ones(1,N),Wind_Speed);
[zsupv,zsuph] = zreflmodv(Frequency,Transmittance,57.*ones(1,N),Wind_Speed);

%tangent curve for Zenith angles >56 degrees
ind=find(Zenith_Angle>56);
zreflmod_v(ind)=zinfv(ind)+(zsupv(ind)-zinfv(ind)).*(Zenith_Angle(ind)-56);
zreflmod_h(ind)=zinfh(ind)+(zsuph(ind)-zinfh(ind)).*(Zenith_Angle(ind)-56);

% negative values set to zero
zreflmod_v(zreflmod_v<0)=0;
zreflmod_h(zreflmod_h<0)=0;

return


function [zreflmod_v,zreflmod_h] = zreflmodv(Frequency,Transmittance,Zenith_Angle,Wind_Speed)

 transmittance_limit_lower = 0.00001 ;
 transmittance_limit_upper = 0.9999  ;
 Transmittance(Transmittance<transmittance_limit_lower)=transmittance_limit_lower;
 Transmittance(Transmittance>transmittance_limit_upper)=transmittance_limit_upper;

 t_c = [0.199277E+00, 0.166155E+00, 0.153272E-01, 0.399234E+01,-0.130968E+01, -0.874716E+00   ,-0.169403E+01   ,-0.260998E-01   , 0.540443E+00   ,-0.282483E+00   ,    -0.219994E+00   ,-0.203438E-01   , 0.351731E+00   , 0.208641E+01   ,-0.693299E+00   , 0.867861E-01   , 0.619020E-01   , 0.595251E-02   ,-0.475191E+01   ,-0.430134E-01   , 0.248524E+01   , 0.388242E-01   , 0.194901E+00   ,-0.425093E-01   , 0.607698E+01   ,-0.313861E+01   ,-0.103383E+01   ,-0.377867E+01   , 0.180284E+01   , 0.699556E+00   ,   -0.506455E-01   ,-0.262822E+00   , 0.703056E-01   , 0.362055E+01   ,-0.120318E+01   ,-0.124971E+01   , 0.154014E-01   , 0.759848E-01   ,-0.268604E-01   ,-0.802073E+01   , 0.324658E+01   , 0.304165E+01   , 0.100000E+01   , 0.200000E-01   , 0.300000E+00 ];
  
 t_c=repmat(t_c',1,length(Frequency));


  
  ZERO     = 0  ; 
  ONE      = 1  ; 
  PI = 3.141592653589793238462643383279   ;
  DEGREES_TO_RADIANS = PI/180.0   ;


        cos_z = cos( Zenith_Angle*DEGREES_TO_RADIANS );

%        !Using the Cox and Munk model to compute slope variance
        variance = 0.00512 * Wind_Speed + 0.0030;
        varm     = variance * t_c(43);
        variance = varm .* ( t_c(44) .* Frequency + t_c(45) );
        
        variance(variance >= varm )= varm(variance >= varm );
        variance(variance <= ZERO ) = ZERO;

%        !Compute surface to space optical depth
        opdpsfc = -log(Transmittance ) .* cos_z;

%        !Define nine predictors for the effective angle calculation

        zx = zeros(9,length(Frequency)); 

        ONE = ones(1,length(Frequency));
        zx(1,:) = ONE;
        zx(2,:) = variance;
        zx(4,:) = ONE ./ cos_z;
        zx(3,:) = zx(2,:) .* zx(4,:);
        zx(5,:) = zx(3,:) .* zx(3,:);
        zx(6,:) = zx(4,:) .* zx(4,:);
        zx(7,:) = zx(2,:) .* zx(2,:);
        zx(8,:) = log(opdpsfc);
        zx(9,:) = zx(8,:) .* zx(8,:);

        zrough_v = ONE;
        zrough_h = ONE;

        for i = 1:7
           j = i-1;
%          !Switched h to v Deblonde SSMIS june 7, 2001
          zrough_h = zrough_h + zx(i,:) .*(t_c(1+j*3,:) + zx(8,:).*t_c(2+j*3,:) + zx(9,:).*t_c(3+j*3,:) );
          zrough_v = zrough_v + zx(i,:) .*(t_c(22+j*3,:)+ zx(8,:).*t_c(23+j*3,:)+ zx(9,:).*t_c(24+j*3,:));
        end
        
        zreflmod_v = (ONE-Transmittance.^zrough_v)./(ONE-Transmittance );
        zreflmod_h = (ONE-Transmittance.^zrough_h)./(ONE-Transmittance );

return







%------------------------------------------------------------------------
%
%------------------------------------------------------------------------

function [ W1, b1, W2, b2, ti, to ]   = nn_isotropic;

ti.gain = [    2.8591851e-03
   2.2471910e-02
   4.0000000e-02
   6.2500000e-02
   5.0632911e-02
   2.3269984e+00
   2.2429120e+00];


ti.xoffset = [    5.0000000e-01
   0.0000000e+00
   0.0000000e+00
   2.7115000e+02
   0.0000000e+00
   1.4029873e-01
   3.5592941e-03];

ti.ymin  = -1.0000000e+00;

to.gain = [   1.9594736e+00
   3.2770936e+00];

to.xoffset = [  -2.8205421e-01
  -1.1063987e-05];

to.ymin  = -1.0000000e+00;



W1 = [  -9.2979879e-01   4.4571587e-01  -1.7062251e+00  -8.3148321e-03   8.8258917e-03   6.2399664e-01  -8.2878242e-01
   4.1461704e+00   2.4964055e-01  -4.4311678e-01  -3.9028611e-02   1.9626402e-02  -1.4928297e-01  -1.7351189e-01
   6.7916759e-02  -9.9813995e-02  -4.5443032e-01   2.7365763e-02   3.9750499e-03   2.3763076e-01   8.2098597e-01
  -1.4722619e-02  -1.3479225e+00   4.8354880e+00   8.7738758e-03  -4.2596351e-03   2.2636867e-01  -4.2103878e-01
   8.4057173e-02  -4.7387358e-01   2.7071709e-01  -2.9345094e-02   1.3307743e-02  -3.5555809e-01  -4.4159311e+00
   1.4495187e+00  -7.5977867e-01   1.0110700e-01   2.3144523e-02   8.8411296e-03  -3.3403279e-01   2.4711803e-01
  -5.1312804e-02  -1.1151146e+00   2.4704847e+00  -6.3171809e-03  -1.5071587e-04   4.2148163e-02  -3.5330817e-01
   5.7952974e-01  -5.6962362e-01   1.5617670e-01  -1.7928660e-02  -2.6088910e-03  -1.0404523e-01   3.0656032e-03
   5.0206411e-02   2.8451973e-01   2.0794085e+00  -1.1312285e-02  -6.2077850e-03  -8.3137033e-02   5.5121670e+00
   8.3190051e-02  -6.2685587e+00   1.2066237e+01  -2.9563476e-03  -1.2024903e-03   4.3895227e-01  -1.0780342e+00
   1.3384641e-03  -2.8723239e-01   6.8066845e-01  -7.9047835e-03   4.7113356e-04  -9.2635657e-02  -1.1558706e+00
  -5.2302845e+00   5.3930237e-01   1.9852186e-01  -1.4673135e-02  -1.9877275e-02   7.6745861e-02  -1.7787875e-01
  -3.5355859e-02  -2.4744932e-01   8.4479142e-01  -1.1038183e-02   1.0534527e-03  -7.6277296e-02  -6.5378632e-01
   1.9358369e+00   5.0596129e-02  -5.7432370e-01  -5.3975200e-02  -3.8812749e-02  -9.9739823e-01  -8.0326549e-01
  -2.4383054e+00  -1.5206303e-01  -4.3148399e-01  -9.8161160e-02  -2.7414470e-02  -1.3494500e+00  -2.1084344e+00
  -1.1116440e+00   3.0285878e-01   3.8813813e-01  -9.6929317e-02  -1.3996950e-02  -1.6737284e-01   3.6000308e-01
  -2.2898453e-01   3.9213364e-01  -1.9554059e+00   2.1531371e-02  -1.8965496e-03  -8.9125974e-02  -7.1204184e-01
   3.2346468e+00  -2.3472668e-01  -6.7730977e-01   2.3401243e-02   2.2066186e-02  -4.5424809e-02  -5.7101693e-04
   8.0897244e-02  -8.6784862e-02   3.5464727e-01   2.7974663e-03   4.5356615e-03  -1.4308270e-02  -3.7655189e+00
  -2.5835578e-01  -1.5092543e-01  -8.9982824e-02   6.3071881e-02  -4.3308953e-02   5.0784082e-01  -2.5622071e-01
   1.9499358e-03   1.7483992e-01  -2.5435567e+00  -1.3391246e-02  -1.9420712e-03  -1.7160377e-01   2.6361945e-01
  -3.3691936e-02   8.7087075e-01  -2.7508919e+00   1.5650491e-02   2.3903542e-03  -4.1324623e-01   8.5276068e-01
   1.0516539e+00  -3.3535385e-01   1.0954073e-01   8.0777310e-02   2.5212141e-02   8.4196448e-01   1.5010275e+00
   4.8585223e-01   7.5301514e-01   2.9529383e-01   4.2477442e-02   2.7668819e-02  -1.9328132e-01  -2.5883452e-01
   3.5496171e-02  -2.8820117e-01   7.5778944e-01  -1.9802679e-03  -4.6233189e-04   1.0419740e-01  -1.2325865e+00
   1.6666296e-01   2.1438626e-01  -1.9662614e-01   1.5964113e-01   5.5616525e-02  -4.6228809e-02   1.0292634e+00
   4.6557264e-03   4.4933589e-01  -1.4034206e+00  -4.9456415e-04   1.0745853e-03   3.9598937e-01   5.1946898e-01
   2.1976360e+00   3.9240255e-01  -7.8704339e-01   4.9924572e-02   2.0881443e-02   3.1885411e-02   3.3566105e-01
   2.2376790e+00   1.0192159e+00   5.7181550e-02   1.6352632e-02   7.0687888e-02   3.9149583e+00   9.7118252e+00
   4.8293830e-01  -9.3600977e-01  -6.9094258e-01  -1.0787445e-01   1.8937550e-02  -2.1663415e-01  -1.4295595e+00
  -3.1099939e-01  -1.8654057e-02  -1.3377822e-01  -2.4743606e-02   1.4474119e-03  -4.2031534e-01  -2.3298025e+00
  -2.9826216e+00   3.3107815e-02   6.5132501e-01  -5.1206533e-02  -2.1932144e-02   3.4169608e-02  -1.8734167e-01
  -2.1426634e-02  -3.7918600e-01   4.2002889e-01  -2.4305402e-03   2.9280851e-03  -4.7285640e-03  -3.1390404e+00
  -4.3705507e-01   5.1661650e-01  -6.5704158e-02   9.6857810e-02   9.4517490e-02  -4.5606367e-01   2.3533582e-01
  -2.2897180e-02  -5.6534599e-01   2.6546113e-01   2.7476675e-02   1.4060951e-03   1.1378572e+00  -1.4845594e+00
   1.9457110e+00   4.0806453e-02  -5.7859371e-01  -5.8536964e-02  -3.8473931e-02  -1.0524793e+00  -8.2195268e-01
  -3.0121716e+00  -9.5233990e-02  -9.9346172e-01   1.8606577e-02  -1.4377912e-02  -4.3904801e-02  -2.4947284e-01
  -7.0001092e-02   4.1602810e-01   3.7603031e-01  -1.4220068e-02  -3.0816121e-03  -5.4705729e-01   1.8739568e-01
  -1.4947763e-01  -7.5608312e-01  -4.4073418e-02  -8.1712325e-03  -2.9563576e-03  -3.2785903e-01  -7.1430632e-03
  -8.2299071e-01  -2.5717835e-03   3.2371193e-02  -1.0848985e-01  -2.9422600e-02  -1.0551377e+00  -1.6372524e+00
   3.0196511e-01   6.9227767e-01   2.0394639e-01  -3.1746292e-02  -1.1676353e-02  -2.9033413e+00  -2.5570774e-01
  -5.7096230e+00   4.2305987e-01  -1.5475364e-01   1.6798078e-02  -1.4993857e-02   6.4242757e-02  -1.5786027e-01
   6.4541847e-01  -4.2286948e-01   4.8489840e-02   1.5806858e-02  -1.4042590e-02   3.2794425e-02   4.1871133e-02
  -1.5185166e+00  -2.2120500e-01   1.0565461e+00   4.5187863e-02  -9.6709381e-02  -2.6591743e+00  -5.3377191e+00
   2.9818048e-01   6.5599532e-01  -7.3518215e-03  -8.4895008e-02  -2.5191490e-03  -8.0881528e-01   5.9587761e-01
  -6.4435723e-02   6.2299747e-01  -1.9091091e+00   9.2290035e-04   1.8096769e-03   1.2153726e-01  -9.5780758e-02
  -2.0464116e-01  -2.0007676e-01   5.0562330e-02   4.2543679e-02  -6.5675256e-03  -1.5208592e-01  -8.3026800e-01
  -4.2637154e-02  -1.7761359e-02  -2.9713148e-01  -2.0645428e-02  -5.4515034e-03  -9.9346491e-01   7.9045873e-01
   1.0280711e+00  -1.0229325e+00  -3.9103399e-01   9.2261227e-02   7.4560034e-02   3.8579725e+00   5.5211981e+00
   2.3637511e-02   1.4871169e-01   3.4777176e-01   4.7626590e-03  -3.5521106e-03  -1.3568540e-01   1.1687029e+00
   5.4443009e-01  -2.4545695e-01   1.1245371e-01  -4.5982947e-02  -8.8125449e-02   2.4591063e-01   1.3175241e-01
   2.0773224e-02  -1.8521673e-01  -1.1527595e+00   5.6313631e-03   4.2924693e-03   1.2302489e-01  -3.3533704e+00
   4.7746255e+00  -3.0113040e-02   1.0724819e+00   8.7356767e-03   2.5642492e-02   6.6632662e-02   1.7953296e-01
  -1.1092410e+00  -2.7215496e-02   1.8294548e-01   1.4429474e-03   6.4980061e-03  -5.3948577e-01   1.7014427e-01
  -1.0374110e-01  -4.2921459e-02  -1.1617897e-01   1.0943411e-02   2.4512365e-04   4.2070047e-01  -9.7117565e-01
   1.0182716e-02   5.5804124e-01   4.5140839e-01   3.2470199e-02   9.3234071e-03   9.7838574e-01  -3.9948599e-01
  -1.4606294e-01  -2.1172088e+00  -1.5195805e-02   3.9726073e-03   4.7688248e-04   8.2993992e-02   1.9979563e+00
   5.8788527e-02   1.4300600e+00  -4.2042164e+00  -9.2430923e-04   3.0364616e-03   3.3918820e-02   1.5358186e-01
   5.4441250e+00   1.8000032e-01   5.5105031e-01   8.5757623e-03   1.8400092e-02   3.2333442e-02   2.6270002e-01
  -5.2908915e-02  -3.6401674e-02   6.1747335e-01   3.9771111e-04  -1.2030089e-03  -1.1917454e-01  -5.0649838e-01
   5.1652739e-02  -2.5015322e-01   6.6440742e-01   1.5296889e-03   2.5891434e-03   8.1299763e-02  -1.7084275e+00
   7.9470727e-01   6.0773983e-01  -4.9363612e-01  -1.1723786e-01   4.9607568e-02  -1.2793774e+00  -1.3306503e+00
   1.0390464e-01   2.6731256e-01  -2.3160785e+00  -2.0838671e-02   5.5952614e-03   2.9460316e-01  -1.1575763e-01
  -2.6303193e-01  -1.1050511e+00  -6.2795343e-02   1.1644727e-03  -5.9895860e-03   9.2784891e-01  -1.3624272e+00
   2.9974855e-01   9.3754176e-01  -1.1596077e+00   5.6377298e-03   2.2665416e-03  -3.5853236e-01   5.8272201e-01
  -2.4490651e-02   8.9316297e-01  -2.8773242e+00   3.7788737e-03   1.6113553e-03  -2.8475986e-01   5.6880343e-01
   1.0340627e+00  -2.0797560e-01   5.6435599e-02   9.7616174e-02   3.0595238e-02   8.8429854e-01   1.5972117e+00
  -4.2689664e-03  -7.9457931e-01   1.5695643e+00  -4.8436964e-04  -6.6750715e-04  -1.9851849e-01  -5.1283713e-01
  -2.9819096e-02   5.2947095e-01   5.3530538e-01  -2.4821000e-03   9.8354304e-04  -8.4915489e-02  -7.5824008e-01
  -1.9320861e-02   7.7387820e-01  -1.2523638e+00   5.8457695e-03  -2.1181860e-04   1.8870666e-02   6.7158743e-01
   8.3464990e-02   4.7357035e-01  -2.2917511e-02   1.4636141e-02  -8.9632346e-04   8.6534858e-01  -7.6886546e-01
   2.9149816e+00   1.7365248e-01   9.6564375e-01   3.2821645e-01   8.5883631e-02   9.3572728e-02   4.6459762e-01
   7.4685093e-03   9.8836857e-01  -3.4214659e+00   5.1091374e-03   5.7598939e-03   8.4117868e-01   4.2585451e-02
   1.4738944e+00   6.4405319e-02  -1.9416046e-01   1.5502632e-01   4.2862797e-02   5.4174292e-01   6.6159988e-01
  -4.7496584e+00  -3.3934944e-01   2.9499805e-01  -7.3731487e-03   3.4259282e-02   1.6457787e-01   3.9720553e-04
   2.1871525e-01  -1.8518453e+00   5.6368414e+00  -2.5667487e-02  -2.8239602e-05   1.2701574e+00  -6.7686608e-01
  -1.8317044e-01  -8.9662036e-02  -6.5550811e-02   4.1773595e-02  -2.1952845e-03   1.2168807e+00  -3.7250253e-01
  -2.6904104e-02   2.3497769e-01   1.2329240e+00   3.3855034e-02  -9.3408806e-03  -1.0031460e+00  -2.5263504e-01
  -2.4194940e-01   1.3149530e-01   3.9417321e-01   7.8823933e-02  -2.3309398e-03   1.6480540e-01   2.0104613e-01
  -7.1768815e-02  -3.7008693e+00   3.9778344e-01   4.6561350e-03   2.6308628e-03   6.5948234e-01   6.1144307e-01
   2.2215043e+00   1.0054307e+00   3.4702999e-02   1.6652618e-02   7.2188187e-02   3.8852323e+00   9.6230482e+00
   2.2514366e+00   2.4200557e-01  -1.9853490e-01   1.5928052e-01   2.3456124e-02   3.8848921e-01   3.6236636e-01
  -1.7117008e-01  -2.1801923e+00   3.4307230e-02   7.8118381e-03  -2.3502325e-03   2.3473237e-01   8.3641878e-01
  -4.5338604e-01  -6.0928010e-01  -1.3964309e-01  -2.2895488e-01  -1.3436763e-01   5.4702356e-01  -5.4748473e-01
   2.0483505e+00  -5.7822858e-01  -7.2169964e-03   4.6500021e-02   9.7677687e-03   1.0243914e+00   1.4932022e+00
  -4.2693010e-01  -8.7946060e-01  -1.9584795e-03   3.5110213e-02  -6.0163396e-03   9.0905557e-01  -8.0336854e-01
  -3.0183107e-01  -6.4924766e-01  -1.5097909e-01  -3.1918935e-01  -1.6682236e-01   3.3265715e-01  -7.6006474e-01
  -4.2331599e+00   2.1714556e+00   9.9230030e-02   5.2796620e-02   2.0401657e-03  -2.0898005e-02   6.3001552e-01
   1.5159840e-01  -8.4283895e-01  -2.4776936e-01   3.2819751e-02  -1.2658812e-03   1.0980257e+00   1.6805004e+00
   8.8947965e+00   2.5939988e-01  -7.2041861e-01   7.2769131e-02   1.9150598e-02  -8.9783678e-02   2.9837979e-01
   5.1621308e-01  -1.2434515e-01   2.2147832e-01   2.2527534e-02   1.2796891e-02   9.1235528e-01   1.5002131e+00
  -1.7300941e-02  -5.5034623e-01   1.3570241e+00   2.2887989e-03  -4.9925923e-04  -1.5706310e-02  -7.0882810e-01
  -2.7220385e-01  -9.3059844e-01   1.1897795e+00  -3.5110995e-03  -2.8838096e-03   3.2593281e-01  -5.5435555e-01
   8.7382172e-01   1.4812502e-01  -1.2887409e-01   1.7814444e-02   1.7870056e-02   2.2048542e+00   3.3338001e+00
  -1.5122702e-01   2.9480328e-02  -2.6435499e-01   1.2582424e-02  -5.4278491e-02   3.0179754e-02  -3.3596852e-01
  -7.3901977e-03   1.7742340e+00   1.0546910e-01   1.0318889e-02  -2.3597687e-03  -7.3774279e-01   2.8350842e-02
  -5.1968914e+00  -2.2072615e-01  -6.1835098e-01  -8.8215870e-04  -1.8393217e-02  -6.3641053e-02  -1.9807606e-01
  -1.0350136e-01  -2.7094011e-01   8.2072819e-01  -9.6174697e-03   1.3744206e-03  -9.8715825e-02  -2.7115954e+00
  -9.9020621e-02  -1.0322486e+00  -2.4274821e-02  -6.4753824e-03  -1.5082613e-03  -1.9589470e-01   6.6195547e-01
  -2.0361767e+00   6.2064825e-01   1.0239006e-01  -4.4594669e-02   2.0435296e-02  -5.5920129e-01   8.0013844e-02
  -2.3928981e-01   4.2506111e+00   8.4714030e-02  -4.0794880e-02   3.0417362e-03  -2.3060713e+00  -4.1114497e-01
   8.9762861e-04   1.2687510e-01   6.4210203e-01  -6.5214566e-03  -2.1601007e-03  -1.7280536e-01   1.9609870e+00
   2.7529535e-01  -3.6170388e-01   7.2112786e-03   2.6703889e-02   3.9133304e-03   1.7390417e+00  -3.8636709e-01
   2.1448671e+00   2.4970981e-01   3.8689279e-01   8.3547751e-02   2.6255448e-02   1.1292990e+00   2.1076419e+00
  -1.2585349e-01   1.2234767e-01   7.6365493e-01  -8.1991312e-03   1.2017672e-03  -1.4589281e-01  -4.2092569e-01
   1.0392508e-01  -8.9524945e-02   3.4278838e-01   6.5913294e-03   6.0142409e-03  -4.4300079e-01  -2.3605519e+00
   9.1767061e-01   1.3155031e-01   6.3908338e-02  -6.6384375e-02  -3.3526097e-02   2.5365453e-01   4.1522743e-01
   8.3166180e-02  -1.1459832e-01   4.4557944e-01   1.1702505e-03   4.8450681e-03   2.4562842e-02  -2.7920209e+00
   2.8232487e-02  -2.8778147e-01  -1.3564955e+00   1.5725738e-02   6.8140455e-03   4.9603715e-01  -6.8216688e-01
   4.7907167e-01  -1.7537064e+00   5.6709781e+00  -6.0188577e-02  -3.1196981e-02   2.4282497e+00   9.1421971e-01
   3.2472030e+00   5.7756426e-01  -2.3953666e-01   3.0806863e-02   1.3117870e-02   2.5567034e+00   1.4037983e+00
   1.4117163e-01   3.3931331e-01   3.7644001e+00  -1.6417755e-02  -1.6291982e-02  -1.0367724e-01   9.5359116e+00
   5.9436365e+00   3.2583613e-01  -6.1835991e-01   6.3803117e-02   1.4314761e-02   2.2847757e-02   3.4108264e-01
   1.5992651e-01  -1.7994584e-01   4.3447846e-01  -6.3250609e-02   7.9793826e-03  -1.2341413e+00  -1.4876583e+00
  -2.9363224e-01  -2.5880443e-02   2.7357307e-01   2.0298399e-02   1.7152941e-03   1.1520119e+00   1.6544107e+00
  -6.0397770e+00   1.0279589e+00  -6.8207974e-03   6.0692461e-02  -1.3665172e-02  -2.4633027e-01  -9.9479295e-01
   2.1739924e+00   3.7534486e-01  -7.6896243e-01   3.9219001e-02   1.5120752e-02   5.2904174e-02   2.4239223e-01
  -2.0949715e-01  -1.4851704e-01   1.8768047e-01  -1.3760515e-01  -4.3827466e-02   8.9940635e-02  -1.1470065e+00
   2.4426992e-01  -5.1378440e-01   4.0821245e-01   7.2050442e-03   9.1054881e-03  -2.6628527e-01  -1.4994969e+00
  -1.9804572e-03   5.2425243e-01  -1.5743700e+00  -2.6811979e-03  -3.6739941e-04  -1.8461967e-01   5.1528520e-01
  -2.0373919e-01  -5.0731456e-01  -5.6206359e-01   3.8039428e-03   1.2360624e-02  -3.8030987e-02  -3.7917469e-01
   4.4831709e-02  -2.9485185e-02  -4.3494758e-01  -3.1618972e-02  -2.4238884e-03  -4.3010804e-01   6.9159413e-01
   1.1642752e+00   5.9654968e-02  -3.1574179e-01  -1.4829999e-01   2.0069061e-01   3.6196377e-01   9.7201758e-01
   3.5914657e+00   2.1788338e-01  -6.6906126e-01   2.1157231e-02   7.2160777e-03  -3.4903701e-02   4.0905916e-02
   8.7875708e-01   1.0434406e-01  -1.5991373e-01   2.0636785e-02   1.9244743e-02   2.2149412e+00   3.1417915e+00
  -4.1691029e-01  -1.1992734e-01  -1.8207633e-01   1.0199683e-01  -6.6250469e-02   7.0477219e-01   9.2985345e-02
   5.5769254e-02  -2.4142753e-01   5.8342596e-01  -1.7699390e-03   7.7548130e-04   4.5042599e-02  -2.1235827e+00
   2.1427991e-02  -2.9280294e+00   7.8994428e+00   4.8748413e-04  -3.6254758e-03   3.1067939e-01  -4.5639382e-01
   3.0453009e-01   7.2568442e-01   2.2837062e-01  -2.9408074e-02  -1.3506292e-02  -3.0424528e+00  -3.6166996e-01
   1.4150755e-01   2.2402409e-01   1.2852193e+00  -3.6619906e-02  -1.5522820e-03  -3.7891196e-02   5.4947069e-01
   1.4912675e+00   1.2013916e-01  -3.5484044e-01  -1.3316142e-01   1.1117767e-01   2.7167487e-01   9.4568848e-01
  -4.4805077e+00  -3.0695915e-01   2.3975999e-01   2.6152724e-02   1.2709816e-02   1.7936441e-01   1.4198214e-01
   1.4812967e-01  -8.0719872e-01   3.0948956e+00  -4.1726944e-03  -7.8223355e-04   1.1712041e-01   4.6070650e-01
   6.2452553e-01  -5.6707535e-01  -2.5040274e-01   4.0460534e-02   1.1011537e-02   7.2080943e-01   2.8851281e+00
   9.6880520e-01   6.0619615e-01   7.6350070e-01   5.3240039e-02   4.0898133e-03   1.3573856e-01   1.4456979e+00
  -1.1414679e-01  -5.8073787e-01  -1.0894564e-02  -1.8937011e-02  -3.6301772e-03  -1.1665106e+00   3.1690528e-01
   4.0554218e-02  -4.2358954e-02   3.0129605e-01   3.3255194e-03   3.5411400e-03  -1.7518120e-02  -5.4846647e+00
  -1.6625372e+01  -1.4622501e-01   3.9470457e-01  -7.6961313e-02  -3.5586303e-02  -3.7767436e-02  -2.0841255e-01
  -8.8348956e-01  -1.8747583e-03  -4.8757635e-02  -1.3748479e-02   1.7176851e-02  -3.2409206e-01  -3.0875219e-01
   9.7080876e-02  -1.4526595e-01  -2.4930063e-01   1.2444368e-01   3.8978087e-02   4.7829022e-03   8.8201592e-01
   2.7971436e-02  -3.5846265e-01   1.2408907e+00   1.2454455e-03   2.6591852e-04   3.7646532e-01  -6.6582303e-01
  -3.2118420e+00  -8.2005529e-01  -1.3124842e+00  -3.0345845e-01   2.3540895e-03  -1.6979239e+00  -4.6942343e+00
   2.9449105e-02  -2.5847671e-01   6.4465581e-01  -5.9474755e-03   1.6035489e-02   2.7972876e-01   2.1739975e-01
  -9.0417050e+00   5.4598167e-02   3.0661607e-01  -1.0620120e-01  -3.6872650e-02   8.3601561e-02  -5.1719015e-01
   2.3589011e+00   5.8289473e-02  -5.2028314e-01   1.1907135e-01   4.0710086e-02  -4.5598321e-02   5.5453378e-01
  -4.6334351e+00  -2.7236286e-01  -5.9172634e-01  -1.6055487e-04  -1.8411327e-02  -9.0282156e-02  -1.3587436e-01
   4.9285580e-02  -5.7404223e+00   8.9380075e+00  -3.9857950e-03   1.5032942e-03   1.5135236e-01  -8.0159295e-01
   1.5375947e+00   4.0553592e-02  -2.6658804e-01   9.9799508e-03   1.4015807e-02   2.2733468e-01  -5.9703967e-01
  -6.0503517e-02   6.0654278e-01  -2.0841470e+00   1.6256251e-03   9.2718282e-04  -5.2378909e-02  -8.1125525e-02
   2.1633171e-01  -1.8654715e+00   5.7318989e+00  -2.6615986e-02   8.8722901e-05   1.2835792e+00  -6.5161172e-01
   6.5372845e-01  -1.9141577e-01  -6.6260466e-02   3.4861201e-02   1.0987446e-02   2.1475183e+00   9.3697753e-01
   3.0399957e+00   3.8236416e-01  -3.5797355e-01  -1.0780960e-02  -1.1692279e-02  -2.2847213e-01  -9.1613301e-02
   9.4293444e+00   9.8554563e-02   4.0281830e-01  -4.0576527e-02   4.1247762e-02   2.5342129e-01   5.1499208e-01
   5.5255927e-01  -5.6296342e-01   5.4032307e-01   6.2675703e-02   2.7711623e-02  -1.3426634e+00  -8.7670313e-01
  -2.3564122e+00  -1.1447970e-01  -7.3469855e-01   2.1526674e-01   7.2052594e-02  -8.9075952e-02  -3.2813806e-01
  -3.5362080e-01  -1.2954758e-01  -1.6240606e+00   2.1148788e-02   1.5808139e-03  -2.6599960e-01  -8.0578427e-01
   3.6194864e-02  -5.3064466e-01   8.5358553e-01  -6.7360252e-03   7.4860475e-04  -1.1340430e-02  -9.6670016e-01
  -6.8903011e+00   1.4939259e+00  -1.1983082e-01   4.2287723e-02  -3.4447620e-02  -3.6171301e-01  -1.0422722e+00
  -9.0393029e-02   7.0552982e+00  -1.8221334e-02  -2.0932892e-02  -3.3229899e-03  -1.4958469e+00   1.3030012e+00
   3.3959820e-02  -1.6699984e-01   2.7472273e+00   9.5285755e-03   1.8165741e-04   1.3492964e-01  -2.1871687e-01
   1.0358026e+01  -7.2157356e-03  -4.2758227e-01   8.6869658e-02   2.7935838e-02  -4.9186890e-02   3.3489223e-01
   9.4630835e-02  -6.3484158e-01  -1.8724044e-01  -6.6914748e-03  -1.2532880e-02   2.4215824e-01   7.9953086e+00
  -7.6932166e+00  -2.4350150e-01   7.2747125e-01  -8.1991099e-02  -2.3057340e-02   7.6138640e-02  -4.1503507e-01
  -1.5398207e+00  -7.4316277e-02   1.3960971e-01  -1.3919671e-01  -2.8994547e-02  -6.6944952e-01  -1.3421650e+00
  -1.5127416e+00  -2.2669918e-01   1.0241273e+00   4.6600672e-02  -9.6560461e-02  -2.6314054e+00  -5.2977737e+00
  -1.0295157e-01  -3.3260249e-01  -5.7031521e-01   4.0467973e-04   2.0620752e-03  -6.2705115e-02   3.0474430e-02
  -2.6909807e+00  -8.0701894e-02  -9.8291550e-01   5.1024521e-02  -8.1175810e-03  -8.9602048e-02  -1.6869999e-01
   2.5531226e-02  -8.6561149e-01   2.6824531e+00   1.5118301e-02  -3.7793072e-03   7.9248630e-01  -1.0272194e+00
   3.0132789e-03   4.8660360e-01  -1.7552170e+00  -7.9023638e-03   1.9259531e-04  -5.1361150e-02   2.7240168e-01
   1.9307015e-01   3.2952760e-01   4.6465587e+00  -1.8084664e-02  -2.3087786e-02  -1.2398482e-01   1.1627164e+01
  -9.3130358e-01   1.2402023e-01   2.1843419e-01   1.4046520e-01  -2.2412138e-01  -4.1547850e-01  -8.0128690e-01
  -7.7537580e-01   3.3705935e-03   1.1172207e-01  -3.6451607e-02  -1.5155883e-02  -2.1974677e+00  -2.0195108e+00
   1.7480473e+01   1.5198164e-01  -3.3275369e-01   6.0960007e-02   3.8435609e-02   6.5060378e-02   1.6908065e-01
   3.1404407e-01  -1.7195407e+00   5.9401981e+00   4.8960134e-04  -2.5282782e-03   1.5704596e-01   8.2146673e-01
  -2.6963220e-01   1.0316299e-01  -7.8633342e-03   3.9665297e-02   9.0251337e-04   1.3614607e+00   5.5240904e-02
   2.6122879e-02   6.4941070e-01   3.9142009e-01  -5.4612892e-03   5.0998698e-03   4.7325251e-01   1.9401159e-01
  -2.6262442e+00  -9.9456206e-02  -9.1531572e-01   6.5265764e-02   8.8829989e-03  -6.1107446e-02  -2.8902380e-01
  -8.9422231e-01  -3.2824570e-01   1.8171322e-01   2.3476541e-02  -2.5545885e-02   1.2869660e-02  -6.7547089e-01
  -1.0651669e-01  -3.4639355e-01  -1.9249443e+00   1.3923016e-02   3.7440425e-03   6.2531026e-02  -3.4022750e+00
  -7.5583472e-02   3.5735262e+00  -9.2920972e+00   9.6633285e-04   2.8692954e-03  -5.8752509e-01   8.6848566e-01];


b1 = [  -5.5481738e+00
   3.9259262e+00
   1.2903522e-01
   5.7747632e+00
  -3.4955900e+00
   1.9731688e+00
   2.6186546e+00
   9.6486797e-01
   6.4886240e+00
   1.7422257e+01
  -6.2273779e-01
  -5.6164188e+00
   8.8147782e-02
   1.9090410e+00
  -3.8945939e+00
  -2.5731869e+00
  -3.4574890e+00
   2.9375877e+00
  -2.7923329e+00
  -5.9966955e-02
  -2.5470903e-01
  -3.1523754e+00
   1.6734683e+00
  -1.6850978e+00
   4.9600477e-01
   7.1460770e-01
  -1.6503954e+00
   1.7002254e+00
   6.6826480e+00
  -3.3066180e-01
  -2.3343158e+00
  -2.8836894e+00
  -2.8866139e+00
  -7.5490631e-01
   1.2374446e+00
   1.9656897e+00
  -4.2194173e+00
   4.6965353e-01
  -2.5322403e-01
  -1.3409111e+00
   1.4733955e-01
  -6.4585518e+00
   7.7035274e-01
  -3.4051257e+00
   8.4632954e-01
  -2.4466905e+00
   8.7709190e-02
   8.2858345e-01
   2.6014941e+00
   5.7178763e-01
   9.5499461e-01
  -3.5066699e+00
   5.7141639e+00
  -1.2282062e+00
  -7.0709510e-02
  -1.0320635e+00
   3.6946125e+00
  -5.1425969e+00
   5.8039958e+00
   9.6453709e-01
  -8.7771914e-02
   4.4305016e-01
  -3.9167391e+00
   5.3693636e-01
   3.9559449e-01
  -3.0413259e+00
   1.5686924e+00
   1.4590518e+00
  -2.5007534e-01
  -5.1717349e-01
  -9.4267242e-03
   3.1883431e+00
  -5.1271489e+00
   1.7351447e+00
  -4.3630587e+00
   6.2340466e+00
  -1.2541247e+00
   1.4336153e+00
   3.6713399e-01
   4.5302409e+00
   6.6148257e+00
   1.9939988e+00
   1.9188821e+00
   9.7495215e-02
   3.0831234e+00
  -6.7468324e-01
   3.4251731e-01
  -5.9487855e+00
   1.0563463e+00
   8.9569879e+00
   1.8077609e+00
   6.3080465e-01
  -4.1943947e-01
   3.2394749e+00
  -2.9533144e-01
  -1.0049421e+00
  -5.7914343e+00
  -2.1820610e+00
   7.2663183e-01
  -2.5277851e+00
  -3.6921785e+00
   1.7098985e+00
   8.6796439e-01
   3.4847648e+00
   1.1885506e-01
  -2.4258304e+00
   7.8384750e-01
  -1.5209127e+00
  -2.0123938e+00
   8.7896217e+00
   4.0862815e+00
   1.2694204e+01
   5.6236062e+00
  -1.4390538e-01
   6.1379520e-01
  -8.3365086e+00
   1.7487500e+00
  -6.4238095e-01
  -9.0614889e-01
  -1.0806230e+00
  -1.5049265e-02
  -5.2789586e-02
   2.0301096e+00
   3.7255898e+00
   3.1291891e+00
  -7.6105773e-02
  -6.6583291e-01
   1.0344484e+01
   1.2671066e-01
   1.6368579e+00
   2.1274504e+00
  -4.3137736e+00
   4.1388680e+00
   2.6691611e+00
   2.3639628e+00
   2.4324766e-01
  -4.8662650e+00
  -1.7035996e+01
  -8.1819474e-01
   1.0406554e+00
   1.0212438e+00
  -5.4983094e+00
  -1.2130910e+00
  -9.9552048e+00
   2.5364106e+00
  -5.4081412e+00
   1.4700785e+01
   2.2844158e+00
  -2.2517876e+00
   6.3110574e+00
   1.8515557e+00
   4.3082534e+00
   9.8313725e+00
  -1.1260268e+00
  -3.2519601e+00
  -3.0361690e+00
  -2.7540413e-02
  -9.5756237e+00
  -6.6322780e+00
  -9.8585120e-01
   1.1244647e+01
   8.4663933e+00
  -7.9238526e+00
  -2.3590854e+00
  -3.3792361e+00
  -3.8287921e-01
  -3.9592118e+00
   4.3456090e+00
  -9.1350570e-01
   1.5412949e+01
  -1.9441747e+00
  -2.4888740e+00
   1.7883404e+01
   8.5767523e+00
  -1.5459566e+00
   2.0453268e+00
  -3.4888357e+00
  -1.5259789e+00
  -5.5527006e+00
  -1.2275883e+01];

W2 = [   4.0051335e-01  -1.8668008e+00  -2.4449709e+00   3.5524317e+00   1.7079847e+00  -4.3338332e-01  -2.0120590e+00  -1.6663834e+00   2.9625577e+00  -6.1856568e+00  -5.2609299e+00  -2.3839230e+00   8.2927298e-01   1.4376675e+00  -9.9442478e-01   1.3558373e+00  -9.7257268e-01   9.1198314e-01   5.3063500e+00   1.6709361e+00  -9.4159497e-02   5.6358221e-01  -7.4365759e-02   1.0774980e+00   3.4144870e+00  -4.6889404e-01  -2.6256608e+00  -1.3102189e+00   1.3227473e+00  -1.1739720e-01  -2.5768090e+00   2.8534561e+00  -7.1597615e+00  -6.0735648e-01   2.0123086e+00  -1.3939754e+00  -9.8210578e-01   2.3132509e+00   3.4095134e+00  -5.7646066e-01   1.6412997e+00   3.9648273e+00   1.5763106e+00   8.6477652e-01  -6.7708789e-01   3.0694872e+00  -1.2862605e+00  -2.4711500e+00  -2.3673916e-02   2.3111727e+00  -1.1666559e+00  -3.4948257e+00  -7.8943401e-01   3.8086712e+00   1.8992240e+00  -8.7933218e-01   3.2413581e+00   3.5041918e+00   3.0070329e+00  -2.3920377e+00   2.9350273e+00  -3.4652745e-02   2.8893640e+00  -1.4327895e+00  -1.8702297e+00  -2.9362270e+00  -2.5564287e-01  -2.9213633e+00  -5.5681465e-01   2.1842184e+00  -2.1877956e+00  -2.9563491e-02  -1.3953482e+00   1.4025152e+00   1.2642776e+00   3.2676204e+00  -1.7221933e+00   2.4000951e-01  -1.2475734e+00   2.6407788e+00  -1.3801703e+00  -6.5628459e-01   2.0977593e+00  -1.6668152e-01   5.8359843e-01  -1.3127595e+00   7.4200087e-02   2.1051577e+00  -7.2426263e-01  -1.0959044e+00  -4.7162552e-01   3.4118371e-01  -1.9949195e+00  -1.4555299e-01  -8.2389801e-01   2.1690292e+00   2.7450814e+00  -8.1057946e-01   3.8422445e+00  -2.2732156e+00   5.4573421e-01   3.2747448e+00   3.1128990e+00  -1.4632424e+00  -9.9015895e-01  -3.6746345e+00   7.5257562e-01   3.3992296e+00  -1.5430161e+00  -9.3881065e-02   1.0501996e-01   8.7508654e+00  -1.0119535e+00  -3.7868059e-01   7.2557779e-01  -5.2832410e+00   1.5926276e+00  -1.1770401e+00  -2.5987006e+00  -8.0324571e-01  -6.8335252e-01  -1.4425510e+00  -6.3334696e-01   2.5225613e+00  -7.0825419e-01  -4.0343268e-01   3.0935324e+00  -6.5521722e+00  -1.3148801e+00  -1.1637123e+00   5.8300894e-01  -2.8856674e+00  -9.3452184e-01  -7.5804919e-01   6.5970632e-01  -1.0128075e+00   4.2750246e+00   4.2776113e+00  -8.0805663e-01  -1.2334613e+00   3.0375364e+00  -1.7695951e-02   6.1782596e-01  -5.0020094e+00   1.1966881e+00   7.8535005e-01   1.1862574e+01   1.9620570e+00  -5.5064987e-02  -2.9161469e+00  -3.6154611e+00   2.8174288e+00   7.0385122e-01   2.2790921e-01  -3.0873080e-01   3.2653492e+00  -3.4617421e+00   3.0781994e+00   3.9455831e+00   8.5267604e-02  -4.2221678e+00  -2.8673239e+00  -4.4405557e-02   1.5407482e+00  -9.2303189e-01   2.3324287e+00  -2.0538062e+00  -4.6282526e+00   8.4558750e-01  -3.4403222e+00  -3.9702581e-01  -3.0031833e+00   4.1715506e+00  -1.6511178e-01   1.3453369e+00  -4.2823682e+00   1.6054503e+00  -2.3218583e-01   9.4098870e+00  -8.9086426e+00
   3.4659814e+00  -1.6868253e+00  -1.4857916e+00  -1.2820893e+00   9.7424298e-01   5.9493796e-02  -4.7365029e-01  -1.1064680e+00   1.1331303e+00   7.9448179e-01  -2.0262909e+00  -1.3684994e+00  -3.5802379e+00   7.9491502e-01  -6.2906869e-01   2.9467405e+00   2.7230751e+00  -4.8486171e-01   1.8955042e+00   1.7635870e+00  -1.9145816e-01  -2.8686357e+00   1.5852710e+00   1.1192260e-01   2.5473577e+00   1.4475370e+00  -9.2995167e-02   3.2357207e-01   1.2713028e+00  -1.6474820e-01  -1.7731700e+00  -4.6824645e-01  -2.2250242e+00   9.5490586e-02  -6.0078334e-01  -7.8766890e-01  -4.3287525e+00  -1.1179341e+00   9.4209561e-01  -8.3299613e-01   6.9371591e-01   3.2034936e+00   1.3417669e+00   9.2474657e-01   5.7989613e-01   1.1460032e+00   5.7563722e-02  -8.8332232e-01  -2.3822572e-02   2.2869219e+00  -2.6095592e-01  -1.4361747e+00  -1.7337292e+00   7.7196183e-01  -3.8757090e-01  -3.9464259e-01   7.2202456e-01   3.4954400e-01   2.1592154e+00   3.9559730e+00   1.2861175e+00  -3.1170220e-02  -3.4669003e+00  -4.8484359e-01  -9.8833299e-01   3.4811806e+00  -1.9011638e+00  -7.6287357e-01   3.0237158e-01  -1.9901722e-01  -3.6451027e-01  -4.9558536e-02   2.9152486e-02   1.0492160e+00   1.3630668e+00  -1.4003403e-01  -1.3485711e-01   7.7979248e-02  -8.4341574e-01   7.6531774e-01  -1.3123982e+00  -4.0682382e-01   4.4636549e-01   7.0117189e-02   2.8055727e-01  -2.2104825e-01  -1.1499469e-01  -1.1989322e-01  -2.9172967e-01  -2.1102025e+00  -1.1400687e+00  -1.5901726e+00  -9.7429858e-01   1.4323956e+00  -1.0017952e+00   5.4167462e-01  -7.4976494e-01  -2.9233540e-01   9.7553178e-01  -8.4063147e-01   1.4426585e-01   1.6821396e+00   1.1770997e+00  -8.4413142e-01   1.9697307e+00  -1.6966270e+00   8.4304108e-01   1.4071191e+00   4.2991387e-01  -1.8649945e-03   4.7934584e-02   3.1911443e+00  -7.5809631e-01  -1.4400161e-01   2.0553449e-01  -1.6014492e+00  -3.4223521e-01   3.5391724e-01  -7.1795515e-01   3.2409437e+00  -1.0619946e+00   9.8927685e-01  -9.0046437e-01   1.3619804e+00  -2.0812317e+00  -5.4215728e-01   1.3494155e+00  -1.7802096e-01  -6.0060793e-01  -7.7160925e-01   7.3288144e-01  -3.2877224e+00   2.5722128e+00  -3.7674298e-01   5.0301896e-01  -2.5785280e-01   1.3768267e+00   7.0291312e+00   7.2057180e-01  -1.8257160e+00   3.2756720e+00  -1.7520966e-02   7.8502008e-01  -1.0111239e+00   2.1915733e-01   3.9864930e+00   1.3903068e-02   1.2842439e+00  -2.8077977e+00   1.2737483e-01  -1.5589265e+00  -5.3101889e-01   5.9568399e-01   1.0714928e-01  -3.8397623e-01   1.7852450e+00  -8.8076044e-01   1.2247539e+00   9.7905802e-01   1.2942658e-01   2.4822568e+00  -9.8259966e-01  -1.7379550e+00   1.0923052e+00  -9.7212427e-01   3.3492789e+00   8.5103540e-02   2.0516115e+00  -1.3152380e+00  -1.2395504e+00  -6.3840442e-01  -1.7059156e+00   6.0866950e+00   1.1415179e+00   2.6117314e-01   1.2156630e-02   1.8643292e+00   1.6979459e+00   3.2320472e+00   1.0162774e+00];


b2 = [  -3.9442164e+00
  -2.1429341e+00];


 

return


%------------------------------------------------------------------------
%
%------------------------------------------------------------------------

function [ W1, b1, W2, b2, ti, to ]   = nn_anisotropic;

ti.gain =  [  2.8591851e-03
   2.2471910e-02
   4.0000000e-02
   6.2500000e-02
   5.0632911e-02
   2.3269984e+00
   2.2429120e+00];

ti.xoffset = [ 5.0000000e-01
   0.0000000e+00
   0.0000000e+00
   2.7115000e+02
   0.0000000e+00
   1.4029873e-01
   3.5592941e-03];

ti.ymin = -1.0000000e+00;

to.gain = [  3.0580307e+01
   4.1162395e+01
   5.3868692e+01
   1.1524030e+03
   3.7019172e+01
   4.0606951e+01
   2.2551726e+01
   5.5613649e+01];


to.xoffset = [  -3.3371218e-02
  -5.7938410e-03
  -3.6754564e-02
  -1.7318159e-03
  -2.8659080e-02
  -2.5723239e-02
  -3.6839388e-02
  -1.0555830e-04
];

to.ymin =   -1.0000000e+00;


W1 = [1.6079875e-01   2.2618774e+00   5.9596311e-01  -2.0646864e-02   4.5072675e-03  -1.8807262e-01   1.4647902e+00
   3.5222099e-01  -1.5481654e+00  -5.1262468e-01   1.0379705e-01   1.8388250e-02   1.1623907e-01   3.0270454e+00
  -3.0854046e-01   1.8590806e+00   7.2450309e-01  -2.3979270e-01  -3.4447214e-01  -4.3861193e-01   1.3050657e+00
  -8.6080599e-02  -2.7518169e-02   3.3337040e+00  -2.8270902e-02   5.0550023e-03  -4.5102365e-01  -2.7272360e-01
  -3.6285667e+00   5.4368324e-02   2.9887411e-01  -2.6459723e-02  -2.9365148e-02  -1.3137931e+00  -1.0761417e+00
  -3.9004637e-01  -3.3464486e-01   1.0186497e+01   3.3043395e-02  -9.1932172e-03  -1.3465028e-01   5.7650332e-01
   6.1770358e-01  -1.1067516e+00  -2.3276503e-01  -9.3718410e-02   1.1961423e-03   1.6472286e-01   6.9975468e-01
  -6.8843285e-02  -4.6798938e-01   1.2724358e+00  -8.6706634e-03   4.6379374e-03  -3.5773808e-01  -2.7821031e-01
  -1.0460106e-01  -1.3723440e+00   5.3119622e-01   1.1424593e-02  -6.5447337e-04   8.4309762e-02   3.4771120e-01
   9.4985616e+00  -1.5918728e+00   2.5987021e-02   1.0633786e-02   7.7179442e-03   5.9388072e-02  -2.9462088e-02
   2.6318170e+00   9.3913624e-01  -9.7368051e+00  -9.4562893e-03  -6.5438754e-03   3.7853291e-03  -2.5397318e-01
  -2.3773439e-01   5.9473028e-02   1.3209130e+00  -1.2147994e-02   3.0958857e-03   8.7052986e-03  -4.1787334e-01
  -2.6048031e+00  -9.4515903e-01   9.2270213e+00   6.0887208e-03   7.2098443e-03   9.0114600e-03   1.2021808e-01
  -4.2655960e-02   2.4152649e-01   1.1870559e+00  -2.2485030e-02   3.2349588e-03  -5.1531301e-02  -6.9110936e-01
  -2.1664572e-01   1.6022443e+00  -2.7323285e+00   5.2682159e-03  -6.0266454e-03  -4.6827645e-02  -2.3024925e-01
  -1.1810330e-01   3.5782705e-01   1.0592177e+00   7.7124278e-02  -1.7991485e-02   6.1113809e-01   8.7369182e-01
   2.6238384e-01   3.0054780e-01   4.0426557e-01   4.3603115e-03  -4.2322534e-03   4.2826991e-01   2.3234085e-01
   7.1461494e-01   2.6614395e-01  -1.1116459e+00   4.6483059e-02  -4.7730454e-04   4.4369170e-01   3.9763091e-01
  -1.1316805e-01   2.5936695e+00   7.6698661e-01   2.1843452e-02  -1.3749144e-02  -2.8288597e-02   2.0513294e-01
  -4.3506379e+00  -9.4513188e-02   1.5758635e-01  -3.4087436e-02  -9.5802198e-03   9.9631587e-01   2.7008118e-01
   7.0382651e-01   1.9206089e+00   8.6843047e-01   6.9748745e-01   5.9918492e-01  -2.6835494e-01   6.9502678e-01
   1.2258529e+00  -6.5418384e-01  -2.8924702e-01   2.2194884e-02  -6.2095270e-03   3.9680645e-01   2.5978740e-01
  -3.0887248e-01  -5.5523260e-01   1.2385260e+00   2.0968990e-02  -1.9594990e-02   5.0251653e-01   3.7260052e-01
  -3.6620292e+00   1.4240720e-03   1.3430822e-01  -1.9763243e-02  -3.8110924e-02  -1.3898374e+00  -2.0400448e+00
   3.2512402e-02  -2.0462659e+00  -7.9456998e-01   2.2816905e-01   3.4511744e-01   5.2639752e-01  -1.2781350e+00
   1.1350052e+00  -8.6509366e-01   4.8294059e-01   4.9675365e-02  -1.1463292e-03   9.1015085e-02   3.1555728e-01
   1.9747621e-01   2.6957835e-01  -3.9324431e+00   3.7593467e-02  -7.2598345e-03   3.7833765e-01   4.4871732e-01
   4.1175797e-01   3.4511354e-01  -1.1148213e+01  -3.2012435e-02   1.0195093e-02   1.5230390e-01  -6.7458364e-01
   1.8397877e+00   6.5738279e-01  -9.7342938e-02  -7.7705972e-02   9.9970274e-02   4.3189337e+00   1.0235943e+01
   1.9737500e+00   3.9923186e-01  -5.7481941e-01  -1.4278632e-02   3.7226083e-04   8.4530794e-02   4.0080352e-01
  -6.4644700e-02  -2.8219701e-01   8.1299789e+00   3.0065927e-02  -3.7814952e-03   4.2228771e-01  -2.8787311e-01
  -1.5918794e+00   4.9012946e+00   3.3618307e-01  -3.5508458e-02  -7.6241708e-04  -1.5023384e-01   9.3306876e-01
  -2.5422505e-01  -1.0721748e-01   5.9522114e-01   1.6862277e-03  -5.2645290e-05   3.3705034e-01  -6.5469004e-01
  -1.1198040e-01   2.8192171e-01  -1.8637637e-02  -8.9886080e-04   2.8439296e-03  -2.2671311e-01  -8.8189749e-01
  -2.2908993e-01  -1.7127547e+00  -5.5002089e-01  -6.4037952e-02   8.5509623e-03  -8.9465500e-02  -6.7452282e-01
  -1.0148470e-01  -2.0990826e+00  -5.7252930e-01   1.6984620e-02  -2.4166820e-03   2.0791146e-01  -1.1712421e+00
  -3.4758059e-01  -5.9801831e-01  -7.1151770e-02  -4.7359511e-02  -2.3574675e-03  -9.7670159e-01  -9.0150372e-02
   3.1433590e-01   1.7322085e+00  -3.6608249e-01  -3.7350106e-02   5.0255179e-03  -2.1117589e-01   3.3231080e-01
   4.9410183e-01  -5.2753059e-01   2.3223385e-01   1.0135782e-01  -2.0056687e-03   5.9207518e-01   6.9061340e-01
  -1.6416437e-01  -1.8479596e+00  -7.8317502e-01  -3.9250649e-02   8.5865597e-03  -2.0415203e-01  -8.5867754e-02
  -2.5995129e-01   2.7344901e+00  -1.5906577e-01  -6.3220574e-03   4.1014275e-05  -5.4942303e-01   2.9501305e-01
  -1.7485018e+00   4.1692016e-01  -8.3496719e-02  -6.1155419e-03  -1.9356531e-02  -1.9815596e+00  -1.4933721e-01
  -1.6447586e-01   1.5682772e+00   2.0103355e-01  -8.0202789e-02  -6.6706712e-03  -1.0235897e+00   1.6325259e+00
   2.1112617e-01  -3.0398529e-01   2.5398266e+00  -7.6298013e-02   4.7741383e-03  -9.2523713e-01  -9.5724224e-01
  -6.6964851e-01   7.2224544e-02   1.8133086e+00   3.6443040e-02  -1.1862656e-02  -1.8971007e-01  -5.6549421e-01
   3.1153390e-01  -8.7217069e-01  -4.2250202e-01   1.1567928e-01  -5.2373420e-04   7.2206721e-01  -3.9936796e-01
   1.3314477e-02  -2.4882179e-01   3.5322208e-02   5.8802339e-02   3.2360662e-03  -2.5990725e-01  -1.2131480e+00
   6.6448535e-01  -3.8463100e-01  -6.4454545e-01  -2.2140412e-02   4.8934508e-03   7.4227459e-02   5.5076881e-01
   6.1861893e+00  -1.0429686e+00  -4.4485456e-01  -1.1658355e-03   4.7340678e-03   1.4087819e-01   4.2677752e-02
   2.1288956e-01  -4.2322132e-01   5.2667405e-01   5.0347901e-02  -5.2724152e-04  -1.1831962e-01   5.6145512e-01
   6.9644630e-01   5.9816609e-01  -9.3801611e-01  -1.7500819e-02  -5.0187300e-03  -2.2980368e-01   9.4949613e-01
   2.5243725e-01  -2.4203677e-01  -5.7901497e-01   3.5981000e-02   8.2133473e-05   2.6699234e-01   8.9853015e-01
  -4.3993805e-01  -8.8496479e-01  -7.5933640e-01  -3.5272068e-02  -4.6216608e-03  -1.5270172e+00   3.6057110e-01
   1.4627998e-01  -1.3883531e+00   3.0882053e+00  -4.9536989e-03   5.8281051e-03   5.4962791e-02   2.5453925e-01
   1.6519810e+00  -7.6930722e-03  -6.4754072e-01   3.0131693e-02   6.2731482e-03   9.2142592e-01  -9.5191066e-02
  -1.4135346e+00  -4.0788537e-01   1.1146971e+00  -2.5350719e-02   1.0625266e-03   2.4353280e-01  -7.7238217e-01
  -2.2644608e-01  -1.5925102e+00  -5.8433845e-01  -4.5540836e-02   6.3523611e-03  -2.0329896e-01  -3.0489873e-01
   1.6891196e-01  -1.3254769e+00  -6.8700782e-01  -3.6475152e-02  -1.0330272e-02   4.1739449e-01  -1.8377843e+00
   3.1358392e-01   1.4968871e+00  -2.5981205e-01  -4.2080671e-02   6.1793007e-03  -1.8505454e-01   1.9953860e-02
  -3.4228045e-01  -6.2317902e-02   6.4669564e-01  -2.9002422e-02   8.2668899e-04  -1.7793678e-01  -5.0671032e-01
   1.2537352e-01   7.1101550e-01   1.0565185e+00   3.1241722e-02  -8.9982163e-03  -2.0377190e-01  -4.7737448e-01
  -8.2188166e-01   5.4989705e-02  -1.1994858e+00   8.4837394e-03  -3.4649026e-03   6.7188061e-02   1.4562030e-01
  -5.7047693e-01   8.8102820e-01  -1.8352287e-01   6.0613440e-01  -1.5635642e-01   3.6537211e-01   1.2002526e+00
  -4.3333058e+00  -1.8267243e-01   5.6086744e-01  -1.6144792e-02   5.9882989e-04  -1.2549651e+00  -1.3567194e+00
   1.9088658e+00   3.9699675e-01  -1.2425028e+00  -2.6615628e-02   5.2343657e-03   3.8868879e-02   5.8236447e-01
   9.1894959e-01   4.5589957e-01  -3.1048652e-01   6.9992590e-03   1.4590729e-03   4.4726574e-01  -1.6747763e-01
   3.3350319e+00   4.1423856e-01  -2.2328622e-02  -9.4015436e-02  -9.0557540e-03  -1.6330006e-02   9.7183259e-01
  -9.1525090e-02  -3.9354426e-02   3.2775176e+00   2.9258834e-03  -1.4897316e-03  -4.4617608e-02  -9.3962129e-02
   3.3008067e-01   9.2493183e+00   1.3920845e+01  -4.7024474e-02   1.1874635e-02   8.1044720e-01  -1.1693553e+00
   5.6623050e-03   6.7797206e-01   1.5649287e-01   4.8580500e-03  -1.5225240e-03  -1.4942454e-01  -8.7583316e-01
  -1.6485203e-01   7.8969649e-01   3.4673307e-01   5.5578158e-03  -1.9087628e-03   1.7779989e-01  -7.3686726e-01
  -9.6786382e-02  -1.4952146e+00  -4.4678739e-01  -2.0155191e-03  -2.9307797e-03   3.7243457e-01  -1.3713751e-01
  -3.5495340e-01  -9.0583286e+00  -1.4011772e+01   5.1885379e-02  -1.4455505e-02  -9.6183686e-01   1.3448964e+00
   7.2421214e-01   9.7908485e-01   7.2568280e-01   5.9735179e-02  -1.3473227e-03   7.9794097e-01   1.1123902e-01
  -1.0183905e+00  -2.7611713e-01   2.2446191e-01   2.9505459e-02   9.4470357e-03  -3.9072049e-01  -1.1009382e+00
  -2.7722975e+00  -9.7507820e-01   1.0840342e+01   1.2795248e-02   6.4298535e-03  -1.9742605e-02   4.3847440e-01
   2.9509918e-01   5.3657404e-01   6.0625364e-01   7.2662648e-02  -1.0970171e-02   2.4714993e+00   2.7183288e+00
  -5.8518622e-01   8.3419202e-01   6.6950443e-01  -2.2280866e-03  -4.0782852e-03  -6.4667607e-02  -6.3536648e-03
  -4.8352024e-01  -3.1589389e-01   1.3634696e+00   2.4006473e-02  -4.7006776e-03   2.1590283e-01  -8.3513958e-02
   7.2902612e-01  -1.4852140e+00  -3.3273318e-01  -4.3600398e-02  -3.6919939e-03   2.1093352e-01   1.6151122e-01
   3.4538534e+00  -7.2884067e-01  -3.4670734e-01  -1.2534964e-02  -1.4728539e-02   9.3085703e-01   1.1670247e+00
   1.9659030e+00   4.2244828e-01  -1.4221232e+00  -2.3835994e-02   7.7683077e-03   1.3006537e-02   5.2106202e-01
  -1.0928067e-02  -1.2042445e+00   2.3851958e-01  -1.4245877e-01   1.9539081e-02  -2.3457774e-01   6.5719064e-01
  -6.1452125e-01  -3.5069649e-01   1.2840576e+00   2.1365621e-02  -4.1342107e-03   2.6900196e-01  -2.9550839e-01
   5.6717089e+00  -9.9588079e-01   1.2559624e-02  -2.3061350e-02  -2.1572142e-03  -1.1734377e-02   2.6585966e-01
  -8.7405108e-02  -1.0633954e+00   4.3786716e-01  -1.3569645e-02  -1.4197383e-03  -1.3913697e-01  -4.7613916e-02
  -1.0228413e+00  -1.9269256e-01   1.5669876e-01  -1.4283889e-03   2.7469724e-03  -1.8849224e-01   5.6537562e-02
   2.9438366e+00   1.0262538e+00  -1.1899882e+01  -1.4398693e-02  -6.8767317e-03   3.6666529e-02  -5.7889104e-01
   1.6856760e-01  -1.9872358e+00  -1.2788490e-02  -3.1412264e-02   5.8510152e-03   3.6070305e-01   4.2627940e-01
   3.3255310e+00   1.3533000e-01  -3.2824009e-02   3.0147474e-03   4.5471950e-03   1.0910965e-01   1.4082671e-01
  -9.9663403e-02  -2.3110319e-02   1.2858903e+00   4.6771375e-02  -7.9870529e-03   1.2413487e-01  -2.7720104e-01
  -5.8698135e+00   2.2156751e-01   3.2193453e-01   4.0901304e-03  -1.1845102e-03  -1.4098282e-01  -2.4407969e-01
  -2.2292982e-01   5.6937409e-01  -3.3146530e-01  -3.7603024e-03  -1.0021320e-02  -1.7279446e+00  -8.5022216e-02
   3.3293208e-01  -8.7795495e-01  -7.0426608e-01  -1.2144821e-02   2.2809429e-03  -4.3146874e-02   3.4653207e-01
   1.4846534e+00   4.5982792e-01  -3.4226460e-01  -3.0496983e-02   4.7432857e-03   1.2658311e-01   5.9859299e-02
  -1.3026107e+00   2.3389392e-01  -7.9242370e-01   1.0888280e-02   8.0312714e-03   2.5074423e-02  -5.5077438e-01
  -1.4480577e+00   2.0908665e-01  -1.6386263e-01  -1.4302659e-02  -1.5946683e-02  -1.7324662e+00  -3.0951188e-01
  -2.7384937e+00  -5.5727256e-01   3.5617868e-01   7.0385987e-02  -1.0846817e-01  -4.4028502e+00  -1.0209450e+01
  -1.3516685e-01  -1.0997933e+00  -2.1968288e-01  -9.9249111e-03  -2.9461392e-03   1.0336786e-01   1.8323028e-02
  -2.0499075e+00  -4.0893745e-01  -5.0919507e-01  -2.6995807e-02  -8.5962865e-03   1.5069483e-01  -1.5400682e+00
   6.2394362e+00   1.7463809e-01  -6.6899636e-01  -1.7659282e-02   4.3106153e-03   2.4510703e-01   3.9762159e-02
   1.7649987e+00  -4.2709524e-01   2.4061902e-01  -2.5569806e-02  -2.6366725e-03   2.2893827e-01   4.4948266e-01
   8.1036284e-02   6.6120346e-01   4.8673647e-01   1.1907489e-02  -5.8074712e-03  -3.5179059e-01  -1.1582235e-01
   1.7721467e+00   3.5084223e-01   5.9557320e-01   5.0787117e-02   7.5576582e-03   4.7945360e-01   2.1997527e-01
  -2.5268205e-01  -4.9820321e-01   1.7314307e+00   1.2996305e-03   1.3206292e-03  -3.2424413e-02  -5.0259580e-02
  -1.5865574e+00   3.9485138e+00  -3.7467063e-01  -2.8979900e-02  -2.0884769e-02  -1.9068319e+00  -2.9071358e-01
   9.7229121e-01  -1.1895750e-01   1.3447467e+00  -4.2408020e-03   6.4732986e-03   1.8717926e-02  -2.4255445e-01
   3.8278900e-01   8.8663398e+00   1.4145540e+01  -5.7097344e-02   1.7392441e-02   1.1346957e+00  -1.5447234e+00
   2.5534870e+00  -1.0905359e-01  -7.9613403e-01  -2.4742381e-02   8.4734159e-04  -1.1796278e-02   4.9315402e-01
  -5.2674090e+00  -1.5139754e+00   3.1891312e-01   1.3701594e-02   2.7153047e-03   5.7470576e-01  -7.3468957e-01
   5.8025389e-01  -1.2728848e-01   8.5571359e-01  -1.0271110e-02   1.9709157e-03  -1.1776417e-01   6.6009907e-02
   1.9364057e+00   6.7924785e-01  -1.1314317e-01  -7.7587520e-02   1.0192949e-01   4.2664652e+00   1.0127927e+01
   6.9223721e+00  -3.1539108e-01  -2.8607853e-01  -1.5998290e-02   1.2387676e-03   1.2033192e-01  -2.2097383e-01
  -1.1967652e+00   4.9941025e-01  -2.3160287e-01  -3.8932675e-02   3.4739356e-03  -8.7088632e-03  -5.6666844e-01
  -4.6231634e-01  -4.3034870e-01   1.7605224e+00   3.0322699e-02  -4.1439774e-03   2.2025263e-01   1.6742413e-01
  -1.4447185e+00   1.4540939e+00   7.0417687e-02  -8.1699924e-02   1.2814880e-01  -2.0638852e+00  -2.7173949e+00
   3.4342892e+00   3.0283526e-01  -5.9863862e-01  -6.4934600e-03   2.5323418e-03   6.6371087e-02   3.5204736e-01
   3.0496813e+00   5.1594462e-01   1.1926172e-01  -1.0892746e-01  -1.0828120e-02  -3.1402381e-02   9.5904865e-01
   1.5772022e-02  -8.1938270e-01   2.1777591e-01  -1.5856359e-01   2.1256497e-02  -4.5044686e-01   8.6282283e-01
  -2.2610083e-01  -8.4832642e-01  -4.9911966e-01  -2.3988030e-02  -5.8769662e-03  -9.3947968e-02  -6.3707563e-01
   5.7800940e-01  -7.7956164e-01  -8.4645287e-01   2.9740187e-02   2.5479206e-03   2.2560261e-01  -1.8559523e-01
  -1.1008306e+00  -4.3957866e-02  -3.5222550e-01  -3.7303091e-02  -1.5940095e-02  -1.5577762e+00  -2.7699878e-01
   2.2763547e-01   4.5150927e-01  -6.9825639e-01  -4.7077345e-03   3.6301860e-03   3.3856346e-01   3.0022334e-01
   5.3525187e+00  -3.4477142e-01  -2.7404177e-01   1.3902881e-02   3.8404229e-03   1.9234626e-01   2.6193302e-01
   2.7892314e+00   4.8193162e-01  -4.0441038e-01  -6.5838007e-02   1.0649285e-01   4.4059475e+00   1.0144882e+01
   3.2897177e-02  -6.2988709e-01  -2.7318104e-01   9.5424971e-03   1.4671605e-04   1.5481598e-01   8.1925584e-01
  -1.0569768e+01   1.6043315e-01   2.0451487e-01  -1.2223687e-02  -3.8529953e-03   2.5398907e-02  -5.2476400e-01
  -9.9602624e-02  -8.1220423e-02   2.8939389e+00   1.2940634e-02  -2.9982005e-03   8.9178129e-02  -7.0570392e-02
  -3.0643094e+00  -5.8961105e-01  -2.5589000e-01   1.2458181e-01   1.1419512e-02   3.5341193e-03  -8.5045035e-01
  -1.9339978e+00   1.2645542e-01   2.3992058e-01  -6.6835922e-01  -6.6548848e-01   2.3842168e-01   4.9846459e-01
   1.5434109e+00  -9.7134284e-01   3.3137239e-01  -1.6962902e-02   4.0802164e-02   2.2271881e+00   1.9473687e+00
  -1.0322639e+01  -4.8774697e-02   4.2511552e-01   9.7775667e-03  -3.7667640e-03  -9.4116233e-02   4.3730664e-02
   3.4553538e+00   4.0162394e-02  -1.6609587e-01   6.1417247e-02   1.1008723e-02   3.9299927e-02   9.2102024e-01
   4.3618904e-01  -9.4792316e-01   3.4663745e-01   3.4845290e-01  -1.2887969e-02   4.8065402e-01   1.8419046e+00
   9.2612563e-02   9.1183021e-02   1.2179289e+00   5.3212237e-02  -9.0914049e-03   2.8403662e-02  -4.3177264e-01
  -5.1821853e-01   4.4530956e-01   8.8846662e-01  -7.2052685e-02   3.9385289e-03   2.9150314e-02   5.7697333e-02
  -3.8567120e+00   1.0730456e+00   3.1759397e-01  -4.1220797e-02  -1.6684425e-02  -3.8361612e-01   3.8583490e-01
  -6.2492546e+00   6.1277093e-01   3.5643044e-01  -3.2197213e-02  -1.8424504e-03   1.4937898e-02   3.3211943e-01
   1.4921436e-01   1.1380888e+00  -1.5616567e+00  -4.2749709e-02  -1.1946092e-02  -2.8427960e-01  -7.7427071e-01
   1.3416895e+01  -1.9162935e+00  -3.3570045e-01   2.7045638e-02   1.1929815e-02   4.9987911e-02  -2.9342764e-01
  -1.4019670e+00  -5.4688428e-01   1.3408320e+00  -5.8598368e-02   6.4501436e-03   2.4450647e-01  -9.8151589e-01
  -1.5017836e+00   1.1424357e+00  -1.5579048e+00   2.1382967e-03   2.0894064e-04  -7.2896361e-02  -7.3575748e-02
  -6.2943524e+00   5.0122477e-02  -3.2760482e-01   7.4826551e-04  -3.3954832e-03  -1.8500614e-02   1.9750065e-01
  -1.5644496e+00  -7.4290039e-01  -8.0836335e-03   7.5966204e-02  -8.8661101e-02  -4.2611580e+00  -1.0287330e+01
   5.1258397e+00  -1.0862163e+00   2.7119217e+00   3.5596349e-02   1.0841311e-02   3.9582735e-01   7.8304477e-01
  -2.9233144e-01  -2.4946049e-01  -1.0454179e+00  -5.7545646e-04   1.6380075e-04  -4.8729744e-01   9.7973341e-01
   3.3856059e+00   1.9145871e-01   2.0940848e-01  -1.1100401e-02   1.9470627e-03   9.9879627e-02   6.7819845e-02
   1.6676830e-01   1.2058519e+00  -1.0627885e+00  -1.1312255e-01   2.1944208e-02   3.6820707e-01  -1.7703403e+00
   8.9669544e+00  -4.5562427e-02  -4.0135293e-01   1.2178928e-02   2.7061957e-03   2.2541171e-02   3.7095955e-01
   1.6361637e+00  -4.2427618e-01   4.4177180e-01  -2.4881177e-02  -1.5678178e-03   3.6737746e-01   4.7940038e-01
   9.5194850e+00   4.2774126e-02  -4.2485446e-01  -2.9645451e-03   2.4661474e-03   2.1284154e-01  -1.2709351e-01
   4.2595515e-01   6.3755694e-01  -1.7752514e-01   2.7111160e-02  -3.5012063e-03   7.2754277e-01   4.5222062e-02
   1.1056273e+01  -1.7717747e-01   3.8207519e-01   5.5663036e-03   7.8328456e-03  -3.7461562e-02   5.9239614e-01
  -7.5423283e-03   1.1732833e+00  -2.1491057e+00   4.8931552e-02  -5.0474600e-03   6.2398855e-01   7.0518404e-01
  -3.1177752e+00  -7.6700009e-02   6.0933012e-01   2.1598587e-02   8.4656421e-04   6.8363340e-01  -3.1453347e-01
  -1.0952597e+00   6.0963799e-01  -3.8229664e-01  -7.5195960e-02   1.4945549e-02  -1.4568539e-01  -8.8492540e-01
  -6.7409826e+00   5.4373365e-01   5.1803289e-01  -4.6911247e-02   4.2342947e-03   7.0806030e-03   6.4789039e-01
  -5.1624595e-01   2.0271336e+00  -3.7693348e-01  -1.4697737e-02   2.5650676e-03  -2.0502828e-01  -2.6750354e-01
  -1.6807476e+00   1.5460528e+00   4.0632445e-01   2.8651050e-03  -1.7065101e-02  -2.4688866e+00  -1.0651239e+00
   1.2868672e-01  -8.3886060e-01   6.9614495e-01  -3.7124955e-03   1.9950471e-02   4.8280077e-01   1.1619930e+00
   2.1485584e+01  -5.5808095e-01  -3.0405553e-01  -2.6295814e-03   6.6761982e-03  -1.3230198e-02  -1.2218652e-03
   4.0049501e+00  -4.0893770e-01  -6.1664502e-01  -1.4222051e-03   1.8716488e-03   6.8285690e-02   4.7165916e-01
  -6.1436018e+00   2.1324722e-01   4.8455649e-01   1.0263915e-02  -2.5088297e-04  -1.0086594e-01   1.4415640e-01
  -1.1035755e+00   1.8685458e-02  -1.1987114e+00  -1.0426665e-02   3.1721047e-04  -4.9120001e-02   3.3831895e-01
   1.3184448e+00   1.1103819e-02   1.9585301e+00  -1.4692981e-02  -5.5344003e-04  -2.3123311e-01  -1.3907855e-01
  -5.0468373e+00  -1.4711494e+00   2.8534176e-01   1.7639291e-02   2.3135365e-03   2.6186543e-01  -3.8059194e-01
  -7.7612762e-01   8.1410334e-02  -5.2924027e+00  -8.6117847e-03  -9.6643934e-03  -6.9793380e-01   7.7651371e-02
   1.7924276e+00   3.4609152e-01  -2.2386655e-01  -1.2609197e-01   2.2311891e-02   2.4199463e-01  -3.2626452e-01
   1.8132694e+01  -4.5146183e-01  -3.3475570e-01  -3.9589763e-03   3.3711924e-03  -1.4863433e-03  -4.1266582e-02
   6.7874162e-01  -2.1136114e+00  -8.1855552e-01  -2.8540591e-01   9.8956857e-02  -5.1936941e-01   1.7833453e+00
  -1.3901899e-01   2.0408873e+00   3.6183858e-01   1.7812821e-02  -4.0562379e-03  -4.4874367e-01  -1.0155032e+00
  -6.6442184e+00  -9.0677391e-02   4.4708766e-01  -2.2770800e-02  -3.6177709e-04  -5.7697233e-01   3.0174568e-01
   2.1741024e+00   1.1606564e+00   1.6326010e-01   1.7620961e-02  -3.4166061e-02   2.1748681e+00   3.9327744e+00
   2.3133801e-01  -1.8476848e+00   1.1248976e+01  -6.6580815e-03   5.8145610e-03   2.3411861e-01   1.1481861e-01
  -3.0645908e-01   1.7839102e+00   5.1395075e-01   1.6980614e-01  -7.8907915e-04  -1.2943711e+00   2.6529923e+00
   5.7738807e+00  -1.0525285e+00  -3.5962750e-01   1.8503179e-03   4.6964609e-03   1.3983248e-01   1.0083350e-01
  -3.1907853e+00   1.2994350e+00   5.5953439e-01   4.3119008e-02   9.1952340e-05   5.6738849e-01   3.5811843e-01
  -6.6549833e+00   7.7352017e-02  -4.5651096e-01   1.2274670e-02  -4.7222544e-03  -3.5837316e-02   6.3638227e-01
  -1.2403861e+01   3.1261199e-01  -2.5438453e-01   6.8744326e-03  -7.1437247e-03  -6.1324986e-03  -9.1617615e-02
  -2.6574799e-01   1.6075732e+00  -1.0216848e+01   2.3156481e-03  -4.2162684e-03  -1.3920002e-01  -2.2077556e-01];

b1 =  [-1.9641305e+00
   6.0634385e+00
  -2.6917618e+00
   3.4555850e+00
  -5.9022953e+00
   9.8907194e+00
   2.5889827e+00
   1.9383002e+00
   1.7317644e+00
   1.1334977e+01
  -7.3444923e+00
   9.4736195e-01
   7.0209243e+00
   1.2017111e+00
  -4.3056451e+00
  -2.1026714e+00
  -8.1118486e-01
  -5.3078677e-01
  -2.3788259e+00
  -6.6812263e+00
  -3.2967895e+00
   6.7667306e-01
  -1.6325385e+00
  -6.9872842e+00
   2.3400057e+00
   1.9331680e+00
  -3.8812285e+00
  -1.0720582e+01
   6.6562992e+00
   1.0876691e+00
   7.8879763e+00
  -4.9431161e+00
  -1.1948736e+00
  -1.1166555e-01
   9.9005466e-01
   1.7219734e+00
   1.7658274e-01
  -1.0019979e+00
  -1.0358301e+00
   1.1046754e+00
  -2.3190570e+00
  -1.0182093e+00
  -2.0857566e+00
   2.5291120e+00
   3.4281048e-01
   1.0551091e+00
   6.3554417e-02
   2.7562264e+00
   6.7626099e+00
   1.5630068e-01
  -3.1877531e-04
   4.2802999e-02
   1.5832587e+00
   4.2463906e+00
   3.2002853e+00
  -9.3892945e-01
   8.1996227e-01
  -2.0253547e-01
  -7.0112952e-01
  -1.9847849e-02
   7.9983332e-02
  -1.0260324e+00
  -1.3999837e+00
  -2.8719483e+00
   9.1926764e-01
   6.9361864e-01
   3.0518780e+00
   3.1249964e+00
   3.2543236e+00
   2.1205235e-01
  -2.7042040e-01
  -1.3411568e+00
  -3.2227829e+00
  -2.5768499e-01
  -6.3554839e-01
   8.0911607e+00
   3.0096621e+00
  -6.6176106e-01
   5.1063953e-01
   1.9520113e+00
   3.2829173e+00
   9.9638923e-01
   1.0915880e+00
   5.8188203e-01
   6.2728865e+00
   1.0167847e+00
  -4.9122767e-01
  -8.8297594e+00
   1.3037281e+00
   2.9520425e+00
   7.3135492e-01
  -5.6833504e+00
  -2.9146627e-02
   7.6787372e-01
   1.6872276e+00
  -1.7008671e+00
  -9.7060365e-01
  -7.1107119e+00
   3.8287934e-01
  -2.9498688e+00
   5.5186299e+00
   2.2758984e+00
   2.1669435e-02
   1.1412920e+00
   1.3712648e+00
  -3.8569292e+00
   9.3677190e-01
   3.1943288e+00
   2.1029758e+00
  -5.4429616e+00
   1.1624348e+00
   6.6219675e+00
   6.4081209e+00
  -1.5770647e+00
   7.3249017e-01
  -3.0644386e+00
   2.7139675e+00
   2.9930915e+00
   1.1010092e+00
   9.6701796e-02
   8.6045610e-01
  -9.3959018e-01
   6.5706360e-01
   5.3323127e+00
   7.1435944e+00
  -3.2518370e-01
  -1.1124667e+01
   3.2691174e+00
  -3.1589439e+00
  -1.2241433e+00
   2.3287587e+00
  -1.0353912e+01
   3.3138516e+00
   1.0695740e+00
   8.6432872e-01
  -1.0057410e+00
  -4.7141546e+00
  -6.7874471e+00
  -1.5477423e+00
   1.4823472e+01
  -9.7092495e-01
  -3.5744834e+00
  -6.7181225e+00
  -6.5744306e+00
   8.7513815e+00
  -8.1874333e-01
   3.1413491e+00
   3.7530157e+00
   9.7222610e+00
   2.4274300e+00
   9.6353340e+00
  -4.5699711e-01
   1.2294318e+01
  -3.1519799e+00
  -3.2713250e+00
  -1.2403267e+00
  -7.8047774e+00
  -2.6576320e+00
  -4.6254434e+00
  -9.4023298e-01
   2.2251087e+01
   3.9958426e+00
  -5.8683599e+00
  -2.2094257e+00
   3.0892848e+00
  -5.3763005e+00
  -6.2702709e+00
   2.8773983e+00
   1.9486821e+01
   3.8565125e+00
  -3.9363093e+00
  -7.7498658e+00
   4.9322726e+00
   1.2259897e+01
  -8.5269383e-01
   6.4433446e+00
  -4.0478735e+00
  -7.2798289e+00
  -1.3575911e+01
  -1.1371588e+01];


W2 = [  4.4859959e-01   2.8223066e+00   4.0482956e-02   6.9912169e-02  -7.4536136e-02   4.3291102e-03   3.4419222e-01  -3.6709312e-01   7.6206071e-01  -9.1072440e+00  -9.4066498e-02   1.4529595e-01  -4.5819107e-02  -2.6382870e-01  -4.1208972e-01  -9.7060558e-02   7.1764224e-01  -1.8719611e-01   5.6760623e-01   4.1845140e+00   6.3506404e-03   2.4362293e-01  -4.6957685e-02  -5.0472708e-01  -1.1580961e-03  -5.9554233e-01  -6.3566296e-03   2.4752111e-03   9.6077315e-01  -7.1742477e-01  -2.9049034e-02  -6.3915845e-01  -1.4609030e-01   1.0038827e+00   1.2295745e+00   1.1616559e+00  -7.7916158e-01   1.3421996e+00   1.0179922e-02  -3.4462121e-02   9.2488883e-01   1.0581681e+00   1.4420446e+00  -1.8922641e-02  -2.5541064e-02   4.7695137e-01  -2.1147473e-01  -2.9813380e+00   2.7884462e+00  -2.1413481e-01   1.7312873e-01  -1.5926585e-01   1.4065148e-01  -3.4423136e-01   5.4766021e-02  -4.1853211e-01  -1.2052211e+00   5.5363283e-02  -2.1805417e+00  -1.5249846e+00   9.4318370e-02  -1.2843690e+00   2.1143474e-02  -2.2319172e-02   2.7946473e-01  -1.2521918e+00  -6.5496034e-01  -1.8938612e-01   3.1558388e-02  -1.4506692e+00   1.0754242e-01  -1.4243059e-01   6.6021469e-02   4.2364129e-02   1.9551465e-01  -8.2335775e-02   9.2240256e-03   7.0863092e-02  -5.8460816e-01   5.8049987e-02  -5.0906801e-02  -2.2555567e-01  -3.4380123e-01   1.3733513e-01  -2.0327495e+00   4.3353306e-01  -3.6263492e-01  -3.3749763e-02   5.6896203e-01   3.9164291e-01  -8.4013131e-02   2.6077481e+00  -1.6207791e-01  -5.0525727e-01   1.4339459e+00   2.3199040e-02  -1.8410989e+00  -3.5576581e-01  -2.2679005e+00   2.6151035e-01  -3.2893746e-01  -3.7038837e-01  -1.2648387e-01   9.3144113e-02  -1.5828015e-01   1.4862753e-01  -4.8395144e-01   3.5241337e-02   8.1471115e-01  -1.7989250e+00  -1.2391650e+00  -6.2210837e-01   5.3266183e-01   2.7883541e-01   2.2904268e-01  -2.0547484e-02  -9.7349800e-01   5.7335716e-01   4.1763661e-01   2.5054154e-01   9.1282856e-02   7.5684156e-01  -7.4808255e-01   1.1803589e-02  -3.6428217e-01  -8.5162682e-01   2.6755175e+00   2.1941089e-01   1.3087809e-01   1.3848450e-03   5.7558655e-02   6.2006357e+00  -3.5423275e-01   3.0550403e-02   1.7272999e-01   9.9050395e-02   2.4385990e+00   3.8012501e-01   1.8701211e-02   4.3859632e+00   1.5921024e-01   1.9065461e-01  -6.0199661e+00   3.4738889e-01   7.6250719e-02  -1.0137672e-01   6.7396838e-01   1.9608029e-01   8.5224224e-01   4.5031147e-01   6.6080776e+00  -1.2156399e+00  -1.9742957e+00  -3.3935201e-02  -1.9968384e-04   3.4625465e-02   4.9618914e-01  -9.0441498e-01   9.4334845e-02   3.4512551e-02  -3.7576859e+00   1.6319972e+00   5.5617806e-01   9.9099640e-02   4.9869579e-01   1.9756745e+00  -4.4094402e-02  -1.9744616e-01   5.4056657e+00   4.9967483e-02  -4.9261805e-01  -8.3716778e-01   2.8599871e-02   6.9523882e-02  -1.1571101e-01   3.1551422e+00  -1.8776460e-01   4.3839655e+00  -8.7321812e+00   1.0447575e-01
  -1.6078390e+00   9.2029655e-01  -7.5062398e-02  -5.2246315e-02  -8.8482743e-01   1.9558597e-02   8.3185139e-01  -1.0258824e-02  -2.0181898e-01  -7.7213422e-01   1.4395750e-01  -3.2952800e-01   5.9904947e-02   2.8924845e-01   1.5848230e-01  -2.1202967e-01   9.0003354e-02  -1.0622328e-01  -8.4896976e-02   4.1589686e-01   4.1626947e-03   2.1292052e-01  -5.5308707e-03   6.7441792e-01  -5.1714664e-02  -3.0455141e-01  -1.1916887e-02   1.7538881e-02  -1.8069001e-01  -9.4065899e-02   2.6932302e-02  -9.7404402e-03   3.4477905e-01   4.7276853e-01  -6.3111929e-01  -2.0312542e+00  -1.1402331e-01   2.0693920e-01   1.0893809e-01  -1.2473281e-01  -4.6353449e-01   4.4626404e-02  -1.1096911e-02  -6.3697182e-06  -9.1575422e-03   1.4134330e-01  -5.2092260e-02  -1.5019193e+00  -2.8668857e+00  -1.8804857e-01   1.6304726e-03  -2.5505838e-01   3.4579924e-02   1.1470591e-01  -2.9302956e-01  -4.8391206e-02   5.1675446e-01  -9.8902073e-02  -2.4016267e-01  -2.7097559e-01   8.4562918e-02   3.8782829e-01   5.7508455e-03   2.4083206e-03  -5.4644297e-02   3.0446923e-01  -7.7795732e-01   7.6699168e-02   1.8057545e-02   3.5419875e-01   1.9782665e-01  -1.8172036e-01   3.3386690e-02   5.5192740e-02   5.2445912e-03   1.5024763e-01  -3.8806192e-02  -7.6541173e-02   3.3410020e-01  -1.0509918e+00   1.9329794e-02   1.9425021e-02   1.2845824e-02  -1.0427465e-01  -7.0509644e-02   3.0324748e-01   2.7651224e-01   6.6459015e-02   2.0785180e-01   6.0966185e-01  -1.4139250e-01  -4.3336346e-01   6.5592542e-02   6.6968600e-03   3.7149355e-01  -1.3035553e-01  -1.6981536e-01   9.2409103e-02  -4.5199473e-01   9.7359698e-02  -8.7271248e-02   1.6523602e-01  -4.7467723e-01  -8.3611736e-02   1.5218986e-01   4.9504699e-02   1.9270637e-01   1.4925888e-02  -1.0594934e-01  -3.3338333e-01   2.6081065e-01   1.2603156e-01  -2.3500338e-02  -5.2031609e-01  -1.1194131e-01  -2.0491208e-03  -8.2344120e-02   1.6188966e+00   2.3585660e-02   4.4880339e-01   3.3156548e-02   1.4184302e-01  -2.8455628e-01  -1.0645731e+00   8.5979357e-02   9.2646858e-01   7.0809351e-01  -1.8300333e-01   7.5146608e-01   9.5324221e-04   2.2555211e-02   4.8144803e-01   2.6334702e-01   1.2694432e-02   1.1181519e-01   3.0207331e-02   3.1829210e-01   1.0055247e-01   5.2388029e-03   3.6599894e-01   5.8245799e-02   8.2006936e-02  -4.9196032e-01  -6.2984114e-02  -1.1150492e-02   4.3123716e-02  -5.4220069e-01  -1.5405213e-01   1.6367514e+00  -2.1422605e-01   4.3789259e-01   4.6704546e-02  -3.5367862e-01  -1.1736987e-02   4.7009247e-02   1.7070330e-01  -3.9850663e-01   8.1322668e-02  -1.0418382e-02   1.5904179e-02  -1.2014929e+00  -2.9332860e-01  -2.5986642e-01  -2.9515123e-01  -1.4057544e-01   3.7618521e-01   2.2149217e-02  -4.7019404e-02   2.6075644e+00  -2.9045953e-03   2.8820522e-01   7.8588313e-01   2.3420699e-02  -3.8757459e-02  -1.1642744e-01   3.8973345e+00  -5.0161837e-02   6.4132681e-01  -1.4679000e+00  -5.4329757e-02
   1.3913050e+00   2.5554502e-01   2.3341782e-01   4.3311067e-02  -1.6517203e-02  -1.2661819e-02  -1.5282787e-01   2.4044829e-01  -1.9024182e-01  -1.8834697e+00  -2.7189183e-01   2.0369499e-02  -1.1898793e-01   1.3338586e-01   1.4421315e-01  -3.9761860e-01  -4.4383021e-01   1.1152781e-01   3.4890648e-01   1.4791199e+00  -9.0323058e-03   6.0239746e-02  -2.4240699e-02  -6.8538777e-02   1.0287288e-01  -3.6340488e-01   3.0616496e-03  -1.3827566e-02  -1.8769514e-01   4.3191361e-01  -2.7013418e-02  -6.0607900e-02  -1.4406452e+00   1.0763897e-01  -3.4484802e-02   1.5452227e+00   4.9172440e-01   4.6080759e-01   2.6912879e-01  -1.4220075e-01   1.4773013e-02  -6.5637890e-02   7.6715154e-01  -3.0271624e-02  -5.5590710e-03   7.3253741e-02  -9.5816347e-02  -7.4851202e-01   1.1127621e-01  -3.7191460e-01   5.0516238e-02   3.7163015e-01  -7.4629913e-03   7.4808686e-02  -4.8577531e-01  -3.0636014e-02   1.1542839e-01  -6.2942249e-02  -4.4405618e-01   8.2086325e-01   2.3477348e-01  -4.3982648e-01   5.4955114e-03  -3.7512150e-02  -5.6487621e-02   3.0994030e-02   3.1549048e-01  -9.1590303e-02   2.1867303e-03   3.9348909e-01  -1.6453761e-01   3.3956151e-02   1.9994654e-02  -7.0480707e-02   3.0182528e-01  -2.7011379e-01   7.7588316e-03  -7.1845236e-01  -4.0293016e-01   8.7837102e-02  -8.8366528e-02   9.9550051e-02  -1.4507226e-01   2.5631487e-01  -6.1708097e-01  -3.6633421e-01  -5.7323377e-02  -1.1650805e-01  -1.3683936e-02  -2.0308880e+00   5.3534470e-02  -2.5512544e-01  -4.8603432e-02  -7.8075020e-01  -1.9885653e+00  -2.7004452e-02   2.0018232e-01  -2.7188339e-02   8.4937753e-01  -2.7864117e-01   4.8818731e-02  -3.7423990e-01  -2.0438579e-01   1.6006033e-01  -5.4947514e-02   2.1464813e-03  -1.8538934e-01   1.8368469e-02   7.7578477e-02  -6.4476129e-04  -3.0628669e-01   7.4981472e-02  -8.7632304e-01   1.5489848e-01   1.2945158e-01  -1.2570770e-02   6.0247436e-01  -1.2863085e+00   1.1574314e-01  -2.2087976e-01   1.1985682e-01  -1.7892066e-01   7.6540466e-01  -1.3859674e-01  -1.7272802e-02   3.4194251e-01   3.1478368e-01   1.0382686e-01  -8.5963633e-01   7.1542715e-04   2.4055425e-03  -3.0805576e-01  -7.6507371e-02  -1.4160560e-02  -1.0071741e-01   4.2699031e-01   7.9631408e-02  -1.6660133e+00  -2.3986257e-02   8.8154703e-01  -1.4080886e-02  -4.0887982e-02  -1.5345795e+00  -1.1157097e-01   4.3380203e-02  -2.6444288e-01   1.4961022e+00  -9.2665367e-01  -2.7363473e+00   3.6683249e-01  -5.9653465e-01   4.1956684e-01   8.5587402e-02  -5.6988603e-02  -1.8823488e-01  -1.2711785e-01   2.1984369e+00  -5.0319907e-01   3.5302751e-02  -4.4919871e-02  -7.9019459e-01   2.9024430e-01  -6.4934354e-01   3.7055084e-01   2.2596814e-01  -2.1560387e-02  -1.3422798e-02   4.5108569e-01   2.9271183e+00   1.5879095e-02  -4.0952236e-01  -3.3176095e+00  -1.5063379e-02  -1.0298582e-03   4.8374160e-02   6.5030694e-01   2.8762866e-01   8.5506029e-01  -8.2457466e-01  -2.5480933e-03
   3.6001111e-01   4.2678821e-01  -3.1801909e-01  -1.9536704e-02   7.1093998e-02   1.7482733e-02   6.4936350e-01   2.1975109e-01  -3.6859265e-01   1.6860542e+00  -2.1656778e-02  -2.4961925e-01  -2.3642927e-02   3.4526025e-01   2.1531834e-01  -9.7765379e-01  -3.1410641e-01   3.1416108e-02   2.0770769e-01  -4.1238301e-01  -1.4261739e-02  -5.2956247e-02   2.1276132e-02  -9.7331833e-02  -3.3693691e-02   5.8467162e-01  -9.3877842e-02   5.8109066e-03  -1.1071666e-01   4.9561174e-01  -1.0732947e-02   2.4061525e-01  -6.7574913e-01  -5.4334795e-01  -4.7919151e-01  -1.8704158e-01   3.3547886e-01  -6.1173867e-01  -1.1685929e-02   2.6691023e-01   2.6293102e-02  -1.9567929e-01  -1.7079512e-01   7.5225197e-03   2.8678251e-02   8.0624420e-02   6.8286503e-02  -2.2743302e-01  -1.1222872e+00  -2.7346350e-01   8.2628320e-02  -2.2204647e-01  -4.2724836e-02   9.2749607e-02  -3.0119383e-01   1.4717586e-02   4.7868693e-01  -4.5533873e-02   9.2012920e-01   6.3595422e-01   4.2818925e-01   1.6869475e-01   2.8713109e-02   2.7753527e-03  -1.0746538e-01  -6.4435578e-02   1.4819279e-01  -1.4977711e-01   2.5657357e-02  -1.8995190e-01  -2.7453579e-01  -5.6961173e-02   3.4904743e-02  -1.2786091e-02  -1.3973940e-01   5.2147132e-02   3.1488882e-02  -1.3940611e-01  -6.0386814e-02  -6.2645979e-01   2.0444153e-02  -8.9125524e-03  -4.3695473e-02   7.2241515e-02   2.7067389e-01   3.2703752e-01  -1.3451186e-01   4.9580370e-02  -4.5126907e-03  -1.2871867e-01   1.4303622e-01  -7.5716216e-02  -6.3649546e-02   2.4458090e-01  -3.8347491e-01   7.5229029e-03   2.6077015e-01   5.8455513e-02   6.9630582e-01  -2.4537714e-01   3.4784620e-01   3.4205150e-02  -7.3578396e-01  -1.7924453e-01  -1.0706486e-02  -8.2269173e-02   1.0030927e-01   9.5654444e-03  -2.6609673e-01  -1.3857576e-01   7.9908923e-02  -2.8805362e-01  -8.2390334e-01   5.5930343e-01  -4.9715293e-02   3.6335733e-02   4.2505378e-01  -1.4473619e-01  -2.1175644e-01  -9.9326002e-02  -8.7239391e-03  -3.7648207e-02   4.1876647e-01   6.6780517e-01   2.5456612e-01  -4.6087316e-01   6.7583350e-01   7.9813756e-02   1.0967540e-01  -1.8381670e-02  -6.7145631e-02  -1.0557886e+00  -4.9807523e-01   1.3816038e-02   3.7548723e-02   4.1819679e-01  -6.6629675e-02  -6.6850735e-01   6.0487589e-02  -1.7022166e-01  -1.2182799e-01   9.0347354e-02   5.5730087e-01  -2.2155943e-01  -1.1097164e-02   5.1674004e-02  -3.9152652e-02  -7.7146785e-01   1.0073932e-01   2.8261658e-01  -7.1511806e-01   5.3212574e-01  -3.6578809e-02   4.4895418e-02   1.3918759e-01   3.1110425e-02   5.2848658e-01  -3.1927556e-01   2.3300588e-02   6.8770770e-02  -3.7887972e-01   1.3436440e-01  -8.8375172e-01  -1.3086007e-01   8.4649907e-02   2.7234190e-01  -4.3411940e-03  -6.0978514e-01  -1.2180735e+00  -8.6318821e-02   4.2675441e-01  -4.7338873e-02  -1.0028958e-02   9.4878467e-02  -7.2885637e-02  -6.4719313e-01  -8.4685109e-02  -2.2826137e-01  -5.7344535e-01   1.1208573e-01
  -2.5692224e-01   1.5125199e+00   8.0575213e-02  -1.6551629e+00   3.1161446e-01   1.8264800e+00   6.0387086e-01  -3.4014195e-01  -1.8076474e+00  -2.2045395e+00  -1.9345955e+00   2.0289087e+00  -6.2705823e-01  -2.0150464e+00  -1.9858645e+00  -2.7103905e-01  -7.1785085e-01  -2.6810021e-01  -1.9532234e-01   1.5488287e+00  -2.5903210e-02  -2.0915012e-01   4.1901859e-02  -5.6361463e-01   5.1088006e-02  -1.0672438e+00  -7.1009806e-01   1.3656456e+00  -1.8461511e+00  -1.0559280e+00   2.8153243e-01   2.9571038e-01   1.9127826e+00  -3.3792373e-01  -1.6487095e+00   8.9288347e-02   3.8626664e-01   6.4899199e-01   1.5368874e-02  -4.6983484e-01   7.7518207e-01  -5.0578446e-01   1.3858873e+00   1.5463883e-01   2.2296385e-02   4.5024642e-01  -9.0830633e-02  -7.9061143e-01   3.1862625e+00  -3.9747900e-01   5.1634941e-01   1.3487753e-01   4.2037343e-02  -1.9017827e+00  -1.8862159e+00  -1.3849402e+00   2.2187227e+00  -6.7799153e-02   9.4345518e-03  -7.1834575e-01   9.5469588e-01   4.4799622e-01   1.0496194e-02  -7.1303580e-02   9.6979715e-01  -2.5024181e-02  -1.6753796e+00   2.6918309e+00  -2.5014077e+00   4.2500810e-01  -1.2724043e+00  -2.7736589e-02  -4.7133990e+00  -1.5460214e-01   2.5116631e-02  -2.4178884e+00   4.1371843e-02  -1.5864897e+00  -4.0786715e-02  -7.2829677e-01  -8.0132520e-02  -8.5011579e-01  -2.7434567e-01  -3.5720859e-01  -1.4227762e+00   1.2592639e+00   6.9852818e-01  -1.1169920e+00   6.5205363e-01   3.0644056e-01  -1.5402368e+00  -2.6669423e+00  -1.9486222e-01  -2.6189143e+00   1.5501575e+00  -2.2477646e-01   1.0923901e+00   1.0729543e+00  -1.6686655e-01   1.6925443e-02  -5.9043101e-01  -2.5079180e+00  -1.1581953e+00  -2.1296759e-01  -3.0564101e-01  -4.0643122e-02   4.7601026e-01  -2.2129389e+00   1.2854432e+00  -8.2902737e-01  -1.6735483e+00   1.6707779e+00  -4.5211560e-01  -1.2148039e+00   2.3763770e-01  -9.0192816e-03  -1.2452570e+00   2.6483510e+00   1.1634415e-01  -2.7894301e-01  -8.2403459e-01  -4.5590351e-01  -6.4555135e-01  -2.5587565e+00   9.0318602e-01   9.2590753e-02   4.7648231e-01  -4.1785142e+00   1.1421816e+00   6.5282989e-04   1.3260414e-02   2.2673800e+00  -2.6744852e-02  -1.3613834e-02   1.3842803e+00   1.1971031e-01   1.4714367e+00  -2.2001097e+00  -7.3942916e-02   1.3590389e+00   5.4555789e-01  -2.0707342e-01   1.3061223e+00  -3.6657591e-01   2.1585661e-01  -9.3341347e-01  -1.0087318e-02   1.4606519e-01   2.1702319e+00   2.1362041e+00   3.4040008e-01   4.4949115e-01  -3.2672313e+00  -3.0468373e-01  -9.4215852e-01   1.9318774e-01   2.1511598e+00  -2.3460231e+00   1.0012196e-01   1.8804541e-02  -1.8651721e+00   2.1126784e+00   1.1603765e+00  -3.1935876e+00  -1.2283491e+00   8.5411981e-01   1.3561147e-01   1.5230133e-01   3.0981923e-01  -1.4451134e-02  -1.6729691e+00  -3.6270868e+00  -5.5229419e-02   2.0263410e+00   1.3523907e-02  -1.7415885e-01   9.8232003e-01  -2.6047559e+00  -6.0255880e+00   2.6311645e+00
  -1.0178589e+00   2.0304214e-01   4.1192353e-01  -8.4691517e-01  -1.6497784e+00   1.0637069e+00   6.7378525e-02  -1.5312484e-01  -1.2114742e-01  -2.8457813e-01  -1.5719709e+00   2.0276191e+00  -6.2713036e-01  -1.7351759e+00  -1.5608109e+00   1.6106646e-01   1.5816429e-01   1.9443840e-01   2.5119572e-01  -1.1257213e-01  -1.2222921e-02  -1.8812043e-01  -8.6707870e-02   1.2271510e+00   1.8666324e-01   2.3744931e-01  -6.9319324e-01   6.9669835e-01   7.7222864e-01  -2.4162777e-02  -1.5909673e-01   1.6569251e-01  -9.3506730e-01   4.9677121e-01  -1.8218251e-01  -1.2671604e+00   1.0232278e-01  -1.5976154e-02  -8.1115056e-03  -4.3076828e-01   4.2107766e-01  -1.0989360e-01   8.9368367e-01   1.2775301e-01  -1.5592684e-02   1.8473056e-01  -1.7526740e-01  -4.1306839e-02  -1.5006874e+00  -4.7304432e-02   1.0896354e-01   8.0791089e-01  -6.4487036e-02  -1.4705799e+00  -7.4696242e-01   2.2217487e-01   6.5582069e-01   5.5189395e-02   1.4182230e-01   1.1710699e+00   3.3750001e-01   1.1979621e+00   9.0584763e-03   3.6656275e-03   2.6238016e-01  -2.4194241e-02   1.3264095e+00   3.2595990e-01  -1.1755738e+00   7.5924466e-01  -6.1733595e-01   8.7607882e-03  -2.2433776e+00  -5.2869355e-02  -1.1187196e-01  -1.5418570e+00   4.6688317e-03   2.9147520e-01  -8.9041057e-01  -2.3441351e-01   1.4632118e-03  -1.6391210e-01   8.5746791e-02   3.5223881e-01  -1.5691227e-01  -1.5714658e-01   8.8988684e-02  -6.0521269e-01   4.0557343e-02   1.5480291e-01  -1.3446446e+00  -1.3234933e+00  -1.8791071e-02  -5.8021916e-01   1.2052445e-01  -1.4618855e-01   2.6736133e-01  -8.4468846e-01  -1.6087864e-02  -1.2652253e-01   5.4580389e-02  -2.2615703e-01  -7.2208123e-01  -7.6635391e-02  -4.1417422e-01  -1.3678523e-02   6.0300646e-01  -1.0678104e+00   2.0435047e-01   3.2513737e-01   7.5863968e-01  -1.0834605e+00  -7.7435914e-02   3.8998181e-01   3.2387861e-01  -9.1728715e-03  -5.6322875e-02  -2.2844879e+00  -7.9039545e-02  -1.3625847e-01   5.1830474e-01  -1.5745443e-01   2.8149116e-01  -1.0602608e+00  -6.5258222e-01   5.5333785e-01   5.4168506e-01  -4.1524734e-01  -9.0121801e-01   7.5596803e-04   7.5004318e-03  -6.6017770e-01  -2.7181873e-01   2.1283929e-02   1.2839883e+00   2.3830709e-01   7.4540271e-01  -2.5006911e+00   6.9269247e-02   5.3190631e-01  -2.1453568e-01  -2.3761204e-01   1.6789695e+00  -1.2263047e-01   4.3418025e-02  -7.3602433e-02  -2.2932702e-01  -1.3059041e-02  -1.2778069e+00   2.6003490e-01  -2.1065197e-01   1.1092627e-01  -2.1431179e-02  -1.5074173e-01  -5.3175051e-02  -5.9542766e-02   2.1283368e+00  -8.4029174e-01   3.8357756e-02  -6.2985813e-02  -1.0010094e+00   2.3289235e-01   5.6829933e-01  -1.0114262e+00  -5.4098880e-01  -3.6090461e-01  -9.5892004e-02  -4.1762745e-01   2.4483744e+00   1.6352260e-02  -6.6315873e-01   4.5140575e-01   1.7037494e-02   1.3609704e+00   2.0157600e-02   2.0300973e+00   8.1666539e-02  -1.5673218e+00  -6.3772001e-01   1.8977429e+00
   2.9041951e-01   3.0282848e-02  -2.7639940e-01  -3.0412163e-01  -1.9116702e+00  -1.0946505e-01   8.1379080e-01   1.2954400e+00  -5.1277403e-01   1.6352420e+00  -6.0259795e+00   1.5160836e+00  -2.6998833e+00  -1.0090368e+00  -2.2343756e+00   1.9805387e-01   9.4680407e-01   5.4962146e-01   1.3292576e-01  -2.8020518e-01   3.9540986e-02  -3.9303471e-01  -1.4784787e-01   1.5498600e+00  -2.1238867e-01   1.0918320e+00  -3.2074810e-01  -2.1278551e-01   1.0046815e+00   1.0630906e+00  -2.3935158e-01   1.4384003e-01  -1.9503640e+00   2.9533115e-01  -8.0276594e-01   4.1243078e-01  -4.1277520e-01  -6.4123919e-01  -1.1549054e-01  -2.5558993e-01  -2.1400587e-01  -9.5714142e-02  -1.1089254e+00   9.6803009e-02   9.5059759e-02   6.9177742e-02   3.2047889e-01   1.3053623e+00   4.9401038e-01   1.0453220e-02  -8.0681393e-02   7.3637276e-01  -5.9703740e-02  -1.9748803e+00  -2.0447377e-01   2.2181216e-01   7.5623860e-01  -2.2389677e-02   7.7687976e-01   1.9942678e+00  -3.2443634e-01   2.0874809e+00  -1.1621371e-03   5.9854027e-02  -1.1381179e+00  -4.6230907e-01   1.4560423e+00  -3.2265692e-01  -1.4408204e+00   7.6655645e-02  -3.7226905e-01  -3.6861789e-02  -2.7430744e+00   9.3081457e-02  -4.7120548e-01  -5.3013955e+00  -3.8581400e-02   2.3796107e+00  -2.6125626e+00  -1.0238252e+00   1.0613432e-01   8.3376466e-01  -9.3122504e-01   1.2015333e+00   8.1967618e-01  -6.1322148e-01  -7.9668075e-01  -1.9973995e+00   3.0941705e-01   1.9073165e+00   1.1602423e+00  -3.7518506e-01   8.6725398e-02   4.3669736e-02  -1.5387107e+00  -2.2300659e-01   1.0442598e-01  -8.8247695e-01   8.1013209e-01  -6.0246773e-02   5.1883485e-01  -6.8442591e-01   1.8562667e-01  -9.0610109e-02  -1.3259335e+00  -2.6274360e-02   7.8372559e-01  -1.3027352e+00  -8.1302623e-01   4.7444156e-01   2.9372566e+00  -1.2160146e+00   1.5889218e+00   1.3080084e+00   7.9115830e-01   1.4506684e-02   6.4966780e-01  -1.8300539e+00   8.5563242e-01   6.1180666e-01   2.2580805e+00  -3.9691690e-03   7.3202982e-01   3.9394410e-01  -7.0534222e-01   9.7230607e-02  -3.2992015e+00   1.1409947e+00  -6.2501497e-01  -2.3363790e-04  -2.1324168e-02  -5.6555976e-01  -2.3404408e-01   3.3223901e-03  -6.9187025e-01   7.0187434e-01   3.9740902e-01  -5.3295603e-01   2.1322500e-02  -1.9755864e-01  -1.5363311e-01  -8.6530396e-01   3.1915747e+00  -3.4895361e-02   1.0982985e-01   2.6136235e-01  -1.7666795e+00   6.4129366e-03  -5.6415723e+00   2.0304666e-01   1.4966792e+00  -2.8261403e-01  -8.0268150e-01   2.5069764e-01   1.1971231e-01  -1.4010004e-01   4.4190540e-01  -2.1707099e-01  -3.2290309e-03   1.6545951e-02   1.1839735e+00  -9.0527418e-01   1.8511228e+00   4.4234130e-01  -4.1240619e-01  -4.5363886e-01  -1.8507337e-01  -4.3721233e-01  -1.5303662e+00   5.6736665e-03   1.0266359e+00   1.4870760e+00   3.7965559e-02   1.3018426e+00   6.4059970e-02  -2.3869117e+00  -5.7513609e-02  -1.0374138e+00   1.5212467e+00   1.7959294e+00
   1.8406706e-01   9.6179974e-02   2.9923748e-01   3.2559449e-02   1.4086012e+00  -2.4369824e-01   2.3827844e-02   2.6987321e-01   2.9828749e-01  -1.4581791e+00   1.1445773e+00   4.8001090e-02   4.8728036e-01  -1.0988055e-01  -1.4496207e-02   1.2328160e-02   1.6071725e-01  -6.5511239e-02   2.1521590e-01  -8.7500754e-02   3.9437257e-02   1.6053971e-02   3.1599478e-02  -1.0554204e+00   4.0005799e-02  -2.6161481e-01   7.1961776e-02  -1.9811868e-01  -9.5265536e-02  -7.8659622e-02   6.9575801e-02  -1.5004516e-01   3.3551043e-02  -2.2154992e-02   5.9450836e-02   4.8348751e-01  -6.4878823e-02   4.3494442e-01  -1.1185000e-02  -1.0170250e-01  -1.0438476e-01  -2.3664579e-02  -1.8983858e-01  -9.6928878e-03  -7.2938964e-03   1.3886746e-01   1.4591753e-01  -4.0708087e+00  -4.9058969e-01  -7.6822978e-02  -3.1951276e-02   1.2566221e-01   1.2139909e-02   4.2177220e-02   2.2495978e+00  -5.8604879e-02  -5.5026848e-02  -1.0021463e-01  -3.4555191e-01  -1.7470769e-01  -3.9932435e-02   1.8143048e-01  -3.9578721e-03   3.8603038e-02   8.4419513e-02  -4.9810005e-02   1.0669516e+00   1.9285162e-01  -1.7415507e-02  -2.7957004e-01  -1.3836452e-02  -4.4457030e-02  -3.7942613e-02   4.8273428e-02  -2.6338978e-02   1.0877460e+00  -1.0411855e-02  -6.9195356e-02   3.0296699e-01   2.4880948e-01   7.5701726e-02  -8.3652985e-02  -5.5701692e-01  -1.2951810e-01  -3.1654383e-02   2.9917855e-01  -6.4261028e-02   4.3391894e-01  -1.0920119e-02   5.7497164e-01  -2.7208064e-01   2.2926530e+00   2.3629086e-02  -6.6584539e-02   5.8471379e-02  -6.8394923e-02   1.7772778e-03   1.0587923e-01  -3.5519648e-01  -4.9366386e-01   7.3305892e-03  -2.3849785e-01  -1.2073742e-01   5.5600680e-02   5.4123179e-02  -8.4694023e-03   1.1505568e-01  -2.0119830e-02   5.4296163e-02  -2.3934309e-01  -3.3217530e-02   1.1621584e-01   1.3872963e+00  -1.6740549e-02  -1.0027480e-01  -1.0922583e-02  -1.6454686e-01  -2.2103687e+00   5.0524205e-01   3.4619141e-01  -1.9724205e-01   2.1208096e-02   1.7671213e-02   9.5471442e-01   8.1040790e-02  -3.2413877e-01   1.4339394e+00  -3.9999602e-01  -1.2011187e+00   4.6543561e-03  -1.7319917e-02   5.1753090e-01  -3.8790816e-01   8.9212204e-03   1.5651793e-01  -2.1418356e-02  -3.1773439e-02  -9.6587007e-01   1.2800219e-02   6.8967785e-01   4.9629784e-02   6.5557423e-02   4.8581565e-01  -2.0741530e-03  -4.1746564e-02  -2.2575711e-02  -3.3067519e-01  -2.2599649e-02   3.1588667e-01  -1.6985556e-02   1.6628496e+00   6.0465774e-02   2.4014808e+00   4.0469393e-02   1.3584988e-01  -3.8456627e-02   3.7167867e+00   2.7505514e-01  -4.1011010e-02   2.0888218e-02  -3.2710588e+00   1.4413220e-01   5.5513182e-01  -4.7993179e-02  -5.3328000e-02   2.3697276e-01   6.2064911e-02  -1.4667561e-01   1.4233136e+01   8.3083487e-03   3.3526710e-01   1.6534487e+00  -1.0197704e-02  -8.5443607e-02   3.7784463e-02   1.0385661e+00  -1.3692545e-01  -6.6563395e-02   2.4313127e+00  -1.1467764e-01];

b2  =  [-7.9389554e-01
  -1.6530464e+00
   1.2170673e+00
   7.7658185e-01
  -1.2672478e+00
   2.7329894e-01
   2.9603450e+00
  -4.1883279e+00];



return
