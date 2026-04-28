%=============================================================================
%
% FORMAT 
%
%    A function to simulate atmosphere brightness temperatures,
%    based on the Rosenskranz 2017 model. Notice that tbdown
%    includes the cosmic background TBC, but the TBC is also
%    output so it can be removed if desired to do a 
%	    tbdown[noTBC] = tbdown - tra * TBc   
%
% FORMAT   [ tra, tbdown, tbup, origTCWV, origTCLRW, RHO, CLRW, H, TBC   ] = 
%	     fom_atmos_tbs_rosenkrantz( ANGLE, FREQ, SH, P, SP, T, Q, CLRWC, LOWFREQ, TCWV, TCLRW )
%
%
%        
% OUT	                   tra [-]              Transmitance
%			   tbdown [K]           Downward atmospheric brightness temperature 
%			   tbup [K]	        Upward atmospheric brightness temperature 
% 		       	   origTCWV [Kg/m2]     Total Column Water Vapour of the original Q
%			   origTCLRW [Kg/m2]	Total Column Liquid and Rain Water of the
%						original CLRWC 
%			   RHO [Kg/m3]		Surface pressure adjusted Q  converted to density
%			   CLRW [Kg/m3]		Surface pressure adjusted CLRWC converted to density
%			   H [km]		Altitude corresponding to P
%			   TBC			Cosmic background brightness temperature 	
% 		
%
% IN    
%			   ANGLE [degrees]	Viewing zenith angle
%			   FREQ [GHz]		Frequency		
%			   SH [km]		Sensor height from ground	
%
%			   SP [Pa]		Surface pressure	
%			   P [Pa]		Pressure profile 
%			   T [K]		Temperature profile
%			   Q [Kg/Kg]		Specific Humidity profile
%			   CLRWC [Kg/Kg]	Specific Cloud Liquid and Rain Water Content profile 
%
%			   LOWFREQ              Flag, if 1 it will do a simplified calculation with only
%						lower freq lines of the O2 absorption and only the 22 
%						GHz H2O line. If 0, full calculations. With 1
%						is faster but not that accurate, not recommended
%						apart from some specific exercises. 
%
%						and the OPTIONAL altitude-integrated contents. If given
%						the corresponding profiles are adjusted to the
%						colum content, if empty not done.
%
%			   TCWV	[kg/m2]		Total Column Water Vapour 
%			   TCLRW [kg/m2]	Total Column Cloud Liquid and Rain Water		
%
%
%
% REFERENCE
%
%              Adapted from Rosenkranz, P. W. (2017). Line-by-line microwave 
%	       radiative transfer (non-scattering). Code Library. 
%	       https://doi.org/10.21982/M81013
%
% HISTORY
%
%		2019-01-08   Developed by catherine.prigent@obspm.fr and lise.kilic@obspm.fr
%
%=============================================================================


function [ tra, tbdown, tbup, origTCWV, origTCLRW, RHO, CLRW, H, TBC   ] = fom_atmos_tbs_rosenkrantz( ANGLE, FREQ, SH, P, SP, T, Q, CLRWC, LOWFREQ, TCWV, TCLRW )



%===== Unit conversions (if needed) and from structure
%      to original code variable names

% Incidence Angle in degrees
nza         = length(ANGLE);


% Frequency in GHz
nfr        = length(FREQ);


% Pressure profile in Pa

npr       = length(P);


if P(1) < P(end)
  error('Atmospheric profiles should be given from ground up, this one looks like from top pof the atmosphere to ground');
end

% Surface pressure in Pa


% adjusting profiles with surface pressure
% and derivation of density profiles from
% mixing ratios for absoprtion calculations

CLWC       = CLRWC/2;
CRWC       = CLRWC/2;


[P, H, T, Q, CLWC, CRWC, origTCWV, origTCLW, origTCRW, RHO, CLW, CRW ] = fom_aux_atmos_adjust_slp( P, SP, T, Q, CLWC, CRWC ); 


origTCLRW  = origTCLW+origTCRW;
CLRW       = CLW + CRW;


% Number of atmospheric levels
INLEV	   = length(P);



%= P in mBar, H in km, RHO, CLRWC in g/m3 from here on

P     = 1e-2 * P;
H     = 1e-3 * H;
RHO   = 1e3  * RHO;
CLRW  = 1e3  * CLRW;



%=== adapting to sensor height if inside the atmosphere

if SH <= H(end)
  ind = find( H <= SH );
  P     = P(ind); 
  H     = H(ind);
  RHO   = RHO(ind);
  CLRW  = CLRW(ind);
end



%=== adjustment of RHO profile to CWVC

%= column integrating RHO if needed

if exist( 'TCWV', 'var')

  if ~isempty( TCWV ) 

    if isnan( TCWV ) == 0

      %= scaling RHO
      RHO  = (TCWV/origTCWV) * RHO;

    end

  end

end

%=== adjustment of CLRWC profile to TCLRW

%= column integrating CLRWC

if exist( 'TCLRW', 'var')

  if ~isempty( TCLRW ) & TCLRW > 0 & origTCLRW > 0

    %= scaling LW+RW
    CLRW   = (TCLRW/origTCLRW) * CLRW;

  end

end


%=== Cosmic background

% including zero-point fluctuation in cosmic background compensates
% for non-linearity of Planck function

rcons		    = .0479923;
efac		    = exp(rcons*FREQ/2.73);
TBC		    = .5*rcons*FREQ.*(efac+1.)./(efac-1.);
TBC                 = repmat( TBC', 1, nza);




%=== Opacity calculation

OPACITY = nan(npr, nfr);



if LOWFREQ

  % simplifid with only lower freq lines on the O2 absorption and only the 22 GHz H2O line
  % faster but not that accurate 

  for F = 1:nfr

    for I=2:INLEV

      % use the 'absorption-of-averages' method to compute optical
      % depth of each slab; see M.J. Schwartz, Ph.D. thesis pp. 84-87.

      TAV = (T(I) + T(I-1))/2.;
      PAV = sqrt(P(I)*P(I-1));
      WVAV = (RHO(I) + RHO(I-1))/2.;
      WLAV = (CLRW(I) + CLRW(I-1))/2.;
      ABSCOEF =   o2abssim(TAV,PAV,WVAV,FREQ(F)) + abh2osim(TAV,PAV,WVAV,FREQ(F)) + absn2(TAV,PAV,FREQ(F)) + abliq(WLAV,FREQ(F),TAV);
      OPACITY(I,F) = ABSCOEF.*abs(H(I)-H(I-1));

    end

  end

else


  % full calculations all original H2O and O2 lines in the absorprion 

  for F = 1:nfr

    for I=2:INLEV

      % use the 'absorption-of-averages' method to compute optical
      % depth of each slab; see M.J. Schwartz, Ph.D. thesis pp. 84-87.

      TAV = (T(I) + T(I-1))/2.;
      PAV = sqrt(P(I)*P(I-1));
      WVAV = (RHO(I) + RHO(I-1))/2.;
      WLAV = (CLRW(I) + CLRW(I-1))/2.;
      ABSCOEF =  o2abs(TAV,PAV,WVAV,FREQ(F)) + abh2o(TAV,PAV,WVAV,FREQ(F)) + absn2(TAV,PAV,FREQ(F)) + abliq(WLAV,FREQ(F),TAV);
      OPACITY(I,F) = ABSCOEF.*abs(H(I)-H(I-1));

    end

  end

end





%= all layers

tbup = nan( nfr, nza );
tb1  = nan( nfr, nza );
tra  = nan( nfr, nza );



for A = 1:nza

  for F = 1:nfr

    [ tbup(F,A), tb1(F,A), tra(F,A) ] = calcultra( 2, INLEV, squeeze(OPACITY(:,F)), T, ANGLE(A) );

  end

end


% TBdown
tbdown	= tb1 + tra .* TBC;


return





%==========================================================================

function [ tbup, tb1, tra ] = calcultra( IA, INLEV, OPACITY, TEMP, ANGLE )


%=== Microwave radiative transfer TB down and TB up

%  COMPUTES MICROWAVE RADIATIVE TRANSFER
%  PATH 1: PROPAGATION FROM LEVEL NLEV TO LEVEL 1
%  PATH 2: PROPAGATION FROM LEVEL 1 TO LEVEL NLEV
%  REFLECTION AT BOUNDARIES AND THE COSMIC BACKGROUND ARE TO BE
%  CONSIDERED IN THE CALLING PROGRAM, E.G. FOR AN OBSERVER AT LEVEL 1
%  WITH SURFACE AT LEVEL NLEV,
%  TB = TB1 + E1*( (1.-Refl)*Tsurf + Refl*(TB2 + E2*TBcosmic) );
%  OR, FOR AN OBSERVER AT LEVEL NLEV WITH SURFACE AT LEVEL 1,
%  TB = TB2 + E2*( (1.-Refl)*Tsurf + Refl*(TB1 + E1*TBcosmic) );
%  Equivalent to
%  [tb1,tb2,e1,e2] = tbarray2(NLV,H,TEMP,PRES,RHO,CLD,OZONE,1,SECANT,SECANT,FREQ);
%  in the original code



tb1	= 0.;
tbup	= 0.;
tra     = 1.0;

CTHN	= cos(ANGLE/57.296);
SECANT	= 1./CTHN;

% Loops on the atmosphericlevels

for I=IA:INLEV
    % trace path 1 using integral form of RTE
    EM = tra;
    tra = tra  .*exp(-SECANT.*OPACITY(I,:));
    TAV = (TEMP(I,:) + TEMP(I-1,:))/2.;
    tba = TAV.*(EM-tra); 
    tb1 = tb1 + tba;
    % trace path 2 using differential form of RTE
    TRAN_SLAB = exp(-SECANT.*OPACITY(I,:));
    tba  = TRAN_SLAB.*(tbup-TAV);
    tbup = TAV + tba;
end 







return




%==========================================================================

function absor=abh2osim(T,P,RHO,F,do_lowfreq)

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

FL=[22.2351 183.3101 321.2256 325.1529 380.1974 439.1508 443.0183 448.0011 470.8890 474.6891 488.4901 556.9360 620.7008 658.0060 752.0331 916.1716]';
S1=[0.1335E-13 0.2319E-11 0.7657E-13 0.2721E-11 0.2477E-10 0.2137E-11 0.4440E-12 0.2588E-10 0.8196E-12 0.3268E-11 0.6628E-12 0.1570E-08 0.1700E-10 0.9033E-12 0.1035E-08 0.4275E-10]';
B2=[2.172 0.677 6.262 1.561 1.062 3.643 5.116 1.424 3.645 2.411 2.890 0.161 2.423 7.921 0.402 1.461]';
W3=[0.2699E-02 0.2945E-02 0.2426E-02 0.2847E-02 0.2868E-02 0.2055E-02 0.1819E-02 0.2612E-02 0.2169E-02 0.2366E-02 0.2616E-02 0.3115E-02 0.2468E-02 0.3154E-02 0.3114E-02 0.2695E-02]';
X=[0.76 0.77 0.73 0.64 0.54 0.69 0.70 0.70 0.73 0.71 0.75 0.75 0.79 0.73 0.77 0.79]';
WS=[0.1329E-01 0.1478E-01 0.1065E-01 0.1395E-01 0.1440E-01 0.9060E-02 0.7960E-02 0.1301E-01 0.9700E-02 0.1124E-01 0.1358E-01 0.1424E-01 0.1194E-01 0.1384E-01 0.1358E-01 0.1355E-01]';
XS=[1.20 0.78 0.54 0.74 0.89 0.52 0.50 0.67 0.65 0.64 0.72 1.00 0.75 1.00 0.84 0.48]';
XH=[2.60 0.77 0.73 0.64 0.54 0.69 0.70 0.70 0.73 0.71 0.75 0.75 0.79 0.73 0.77 0.79]';
XHS=[1.20 0.78 0.54 0.74 0.89 0.52 0.50 0.67 0.65 0.64 0.72 1.00 0.92 1.00 0.84 0.47]';
SH=[-.3300E-04 -.7200E-04 -.1430E-03 -.1300E-04 -.7400E-04 0.5100E-04 0.1400E-03 -.1160E-03 0.6100E-04 -.2700E-04 -.6500E-04 0.1870E-03 0.0000E+00 0.1760E-03 0.1620E-03 0.0000E+00]';
SHS=[0.8140E-03 0.1730E-03 0.2780E-03 0.1325E-02 0.2400E-03 0.1650E-03 -.2290E-03 -.6150E-03 -.4650E-03 -.7200E-03 -.3600E-03 -.1693E-02 0.6870E-03 -.1496E-02 -.8780E-03 0.5210E-03]';
REFTLINE = 296.; %reference T for lines
%    read continuum parameters; units: Kelvin, 1/(km*mb^2*GHz^2)
REFTCON=300.;
CF=0.59E-09;
XCF=3.;
CS=0.142E-07;
XCS=7.5;

nc = length( RHO );  
FL = repmat( FL, 1, nc);
S1 = repmat( S1, 1, nc);
B2 = repmat( B2, 1, nc); 
W3 = repmat( W3, 1, nc);
X  = repmat( X, 1, nc);
WS = repmat( WS, 1, nc);
XS = repmat( XS, 1, nc);
XH = repmat( XH, 1, nc);
XHS = repmat( XHS, 1, nc);
SH = repmat( SH, 1, nc);
SHS = repmat( SHS, 1, nc);


PVAP = RHO.*T/216.68;
PDA = P -PVAP;
DEN = 3.344E16*RHO;  
% CONTINUUM TERMS
TI = REFTCON./T;
% Xcf and Xcs include 3 for density & stimulated emission
CON = (CF.*PDA.*TI.^XCF + CS.*PVAP.*(TI.^XCS)).*PVAP.*F.*F;
% ADD RESONANCES
TI = REFTLINE./T;
TI2 = TI.^2.5;

I = 1; % only selecting 22.23GHz line

       WIDTHF = W3(I,:).*PDA.*(TI.^X(I,:));
       WIDTHS = WS(I,:).*PVAP.*(TI.^XS(I,:));
       WIDTH = WIDTHF + WIDTHS;
       SHIFTF = SH(I,:).*PDA.*(TI.^XH(I,:));
       SHIFTS = SHS(I,:).*PVAP.*(TI.^XHS(I,:));
       SHIFT = SHIFTF + SHIFTS;
       WSQ = WIDTH.^2;
       % line intensities include isotopic abundance
       S = S1(I,:).*TI2.*exp(B2(I,:).*(1.-TI)); 
       DF      = F - FL(I,:) - SHIFT;
       DF(2,:) = F + FL(I,:) + SHIFT;
       % USE CLOUGH'S DEFINITION OF LOCAL LINE CONTRIBUTION
       BASE = WIDTH./(562500. + WSQ);
       % DO FOR POSITIVE AND NEGATIVE RESONANCES
       RES = 0.;
       for J=1:2
          if (abs(DF(J))<750.) 
             RES = RES + WIDTH./((DF(J,:).^2)+WSQ) - BASE;
          end  
       end
       SUM = S.*RES.*((F./FL(I,:)).^2);
       absor = .3183E-4*DEN.*SUM + CON;   

absor( RHO <= 0 ) = 0; 

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

FL=[22.2351 183.3101 321.2256 325.1529 380.1974 439.1508 443.0183 448.0011 470.8890 474.6891 488.4901 556.9360 620.7008 658.0060 752.0331 916.1716]';
S1=[0.1335E-13 0.2319E-11 0.7657E-13 0.2721E-11 0.2477E-10 0.2137E-11 0.4440E-12 0.2588E-10 0.8196E-12 0.3268E-11 0.6628E-12 0.1570E-08 0.1700E-10 0.9033E-12 0.1035E-08 0.4275E-10]';
B2=[2.172 0.677 6.262 1.561 1.062 3.643 5.116 1.424 3.645 2.411 2.890 0.161 2.423 7.921 0.402 1.461]';
W3=[0.2699E-02 0.2945E-02 0.2426E-02 0.2847E-02 0.2868E-02 0.2055E-02 0.1819E-02 0.2612E-02 0.2169E-02 0.2366E-02 0.2616E-02 0.3115E-02 0.2468E-02 0.3154E-02 0.3114E-02 0.2695E-02]';
X=[0.76 0.77 0.73 0.64 0.54 0.69 0.70 0.70 0.73 0.71 0.75 0.75 0.79 0.73 0.77 0.79]';
WS=[0.1329E-01 0.1478E-01 0.1065E-01 0.1395E-01 0.1440E-01 0.9060E-02 0.7960E-02 0.1301E-01 0.9700E-02 0.1124E-01 0.1358E-01 0.1424E-01 0.1194E-01 0.1384E-01 0.1358E-01 0.1355E-01]';
XS=[1.20 0.78 0.54 0.74 0.89 0.52 0.50 0.67 0.65 0.64 0.72 1.00 0.75 1.00 0.84 0.48]';
XH=[2.60 0.77 0.73 0.64 0.54 0.69 0.70 0.70 0.73 0.71 0.75 0.75 0.79 0.73 0.77 0.79]';
XHS=[1.20 0.78 0.54 0.74 0.89 0.52 0.50 0.67 0.65 0.64 0.72 1.00 0.92 1.00 0.84 0.47]';
SH=[-.3300E-04 -.7200E-04 -.1430E-03 -.1300E-04 -.7400E-04 0.5100E-04 0.1400E-03 -.1160E-03 0.6100E-04 -.2700E-04 -.6500E-04 0.1870E-03 0.0000E+00 0.1760E-03 0.1620E-03 0.0000E+00]';
SHS=[0.8140E-03 0.1730E-03 0.2780E-03 0.1325E-02 0.2400E-03 0.1650E-03 -.2290E-03 -.6150E-03 -.4650E-03 -.7200E-03 -.3600E-03 -.1693E-02 0.6870E-03 -.1496E-02 -.8780E-03 0.5210E-03]';
REFTLINE = 296.; %reference T for lines
%    read continuum parameters; units: Kelvin, 1/(km*mb^2*GHz^2)
REFTCON=300.;
CF=0.59E-09;
XCF=3.;
CS=0.142E-07;
XCS=7.5;

nc = length( RHO );  
FL = repmat( FL, 1, nc);
S1 = repmat( S1, 1, nc);
B2 = repmat( B2, 1, nc); 
W3 = repmat( W3, 1, nc);
X  = repmat( X, 1, nc);
WS = repmat( WS, 1, nc);
XS = repmat( XS, 1, nc);
XH = repmat( XH, 1, nc);
XHS = repmat( XHS, 1, nc);
SH = repmat( SH, 1, nc);
SHS = repmat( SHS, 1, nc);


   PVAP = RHO.*T/216.68;
   PDA = P -PVAP;
   DEN = 3.344E16*RHO;  
   % CONTINUUM TERMS
   TI = REFTCON./T;
   % Xcf and Xcs include 3 for density & stimulated emission
   CON = (CF.*PDA.*TI.^XCF + CS.*PVAP.*(TI.^XCS)).*PVAP.*F.*F;
   % ADD RESONANCES
   TI = REFTLINE./T;
   TI2 = TI.^2.5;
   SUM = 0.;

NLINE = 16; 

   for I=1:NLINE
       WIDTHF = W3(I,:).*PDA.*(TI.^X(I,:));
       WIDTHS = WS(I,:).*PVAP.*(TI.^XS(I,:));
       WIDTH = WIDTHF + WIDTHS;
       SHIFTF = SH(I,:).*PDA.*(TI.^XH(I,:));
       SHIFTS = SHS(I,:).*PVAP.*(TI.^XHS(I,:));
       SHIFT = SHIFTF + SHIFTS;
       WSQ = WIDTH.^2;
       % line intensities include isotopic abundance
       S = S1(I,:).*TI2.*exp(B2(I,:).*(1.-TI)); 
       DF      = F - FL(I,:) - SHIFT;
       DF(2,:) = F + FL(I,:) + SHIFT;
       % USE CLOUGH'S DEFINITION OF LOCAL LINE CONTRIBUTION
       BASE = WIDTH./(562500. + WSQ);
       % DO FOR POSITIVE AND NEGATIVE RESONANCES
       RES = 0.;
       for J=1:2
          if (abs(DF(J))<750.) 
             RES = RES + WIDTH./((DF(J,:).^2)+WSQ) - BASE;
          end  
       end
       SUM = SUM + S.*RES.*((F./FL(I,:)).^2);
       absor = .3183E-4*DEN.*SUM + CON;   
   end


absor( RHO <= 0 ) = 0; 

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


f = f * ones(size(t));
eps=dilec(f,t);
re = (eps-1.)./(eps+2.);
absor = -.06286*imag(re).*f.*clo;
ind  = find( clo <=0 & t < 233 );
absor( ind ) = 0;


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
absor = 9.95e-14*fdepen.*p.*p.*f.*f.*(th.^3.22);
      
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
kappa = -43.7527.*theta.^0.05 +299.504.*theta.^1.47-399.364.*theta.^2.11 +221.327.*theta.^2.31;
%  Debye term from 
%  W. Ellison, J. Phys. Chem. Ref. Data, 36, 1-18 (2007).
delta = 80.69715*exp(-tc/226.45);
sd = 1164.023*exp(-651.4728./(tc+133.07));
kappa = kappa -delta.*z./(sd+z);
%  B band from P.W. Rosenkranz, IEEE Trans. Geosci. & Remote Sens. v.53(3) pp.1387-93 (2015).
delta = 4.008724*exp(-tc/103.05);
hdelta = delta/2.;
f1 = 10.46012+0.1454962*tc+6.3267156E-02*tc.^2+9.3786645E-04*tc.^3;
z1 = complex(-.75,1.)*f1;
z2 = complex(-4500.,2000.);
cnorm = log(z2./z1);
chip = hdelta.*log((z-z2)./(z-z1))./cnorm;
chij = hdelta.*log((z-conj(z2))./(z-conj(z1)))./conj(cnorm);
dchi = chip+chij-delta;
kappa = kappa + dchi;
diel=kappa;
%fprintf('%.2f %.2f %f%+fj \n',f,t,real(diel),imag(diel));

return



%==========================================================================

function absor=o2abssim(TEMP,PRES,VAPDEN,freq)

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
     
     %TEMP=t;
     %PRES=p;
     %VAPDEN=rho;

      FREQ=freq*ones(size(TEMP)) ;;
      TH = 300./TEMP;
      TH1 = TH-1.;
      B = TH.^X;
      PRESWV = VAPDEN.*TEMP/216.68;
      PRESDA = PRES -PRESWV;
      DEN = .001*(PRESDA.*B + 1.2*PRESWV.*TH);
      DFNR = WB300.*DEN;
      SUM = 1.584E-17*FREQ.*FREQ.*DFNR./(TH.*(FREQ.*FREQ + DFNR.*DFNR));

      %ind = find( S300 > 2e-15 & F < 200 );

       NL = 38;  
       for K=1:NL
       %for K=ind
          DF = W300(K)*DEN;
          FCEN = F(K);
          Y = DEN.*(Y300(K)+V(K)*TH1);
          STR = S300(K).*exp(-BE(K)*TH1);
          SF1 = (DF + (FREQ-FCEN).*Y)./((FREQ-FCEN).^2 + DF.*DF);
          SF2 = (DF - (FREQ+FCEN).*Y)./((FREQ+FCEN).^2 + DF.*DF);
          SUM = SUM + STR.*(SF1+SF2).*(FREQ./F(K)).^2;
      end

      absor = 1.6097E11.*SUM.*PRESDA.*TH.^3;
      absor(absor<0) = 0.;


return



%==========================================================================

function absor=o2abs(TEMP,PRES,VAPDEN,freq)

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
     

      FREQ=freq*ones(size(TEMP)) ;;
      TH = 300./TEMP;
      TH1 = TH-1.;
      B = TH.^X;
      PRESWV = VAPDEN.*TEMP/216.68;
      PRESDA = PRES -PRESWV;
      DEN = .001*(PRESDA.*B + 1.2*PRESWV.*TH);
      DFNR = WB300.*DEN;
      SUM = 1.584E-17*FREQ.*FREQ.*DFNR./(TH.*(FREQ.*FREQ + DFNR.*DFNR));

      for K=1:NL
          DF = W300(K)*DEN;
          FCEN = F(K);
          Y = DEN.*(Y300(K)+V(K)*TH1);
          STR = S300(K).*exp(-BE(K)*TH1);
          SF1 = (DF + (FREQ-FCEN).*Y)./((FREQ-FCEN).^2 + DF.*DF);
          SF2 = (DF - (FREQ+FCEN).*Y)./((FREQ+FCEN).^2 + DF.*DF);
          SUM = SUM + STR.*(SF1+SF2).*(FREQ./F(K)).^2;
      end
      absor = 1.6097E11.*SUM.*PRESDA.*TH.^3;
      absor(absor<0) = 0.;


return


%==============================================================================
% surface pressure adjustment
%
% IN    P         Column vector of pressures [Pa].
%	SP        Surface pressure [Pa]
%       T         Column vector of temperatures [K].
%       Q         Column vector of water vapour mixing ratio  [Kg/Kg]. 
%       CLWC      Column vector of cloud liquid water content mixing ratio  [Kg/Kg]. 
%       CRWC      Column vector of rain content mixing ratio  [Kg/Kg]. 
%
% OUT   P         Adjusted column vector of pressures [Pa]. 
%	Z	  Adjusted column vector of corresponding altituded [Km]
%       T         Adjusted column vector of temperatures [K].
%       Q         Adjsuted column vector of water vapour mixing ratio  [Kg/Kg]. 
%       CLWC      Adjusted column vector of cloud liquid water contentmixing ratio  [Kg/Kg]. 
%       CRWC      Adjusted column vector of rain content mixing ratio  [Kg/Kg].
%       TCWV      Adjusted total column water vapor [Kg/m2].
%       TCLW      Adjusted total column cloud liquid water [Kg/m2].
%       TCRW      Adjusted total column rain water [Kg/m2].
%       RHO       Adjusted column vector of water vapour density  [Kg/m3]. 
%       CLW       Adjusted column vector of cloud liquid water density  [Kg/m3]. 
%       CRW       Adjusted column vector of rain content density  [Kg/m3].

function [P, Z, T, Q, CLWC, CRWC, TCWV, TCLW, TCRW, RHO, CLW, CRW ] = fom_aux_atmos_adjust_slp( P, SP, T, Q, CLWC, CRWC ) 

%= first cutting profiles if needed

x       = log10(  P );
xq      = log10( SP );


if xq == x(1) 

  ima   = -1;

elseif xq > x(1)

  ima   = 0;
  P     = [SP;P];

else

  ind     = find( (xq-x) < 0 );
  ima     = ind(end);
  P       = [ SP;P(ima+1:end)]; 

end



if ima ~= -1

  v       = T;
  vq      = interp1( x, v, xq, 'linear', 'extrap' );
        if ima > 0
          T(ima) = vq;
          T      = T(ima:end); 
        else
          T      = [ vq; T];
        end

  v       = Q;
  vq      = interp1( x, v, xq, 'linear', 'extrap' );
        if ima > 0
          Q(ima) = vq;
          Q      = Q(ima:end); 
        else
          Q      = [ vq; Q];
        end


  v       = CLWC;
  vq      = interp1( x, v, xq, 'linear', 'extrap');
        if ima > 0
          CLWC(ima) = vq;
          CLWC      = CLWC(ima:end); 
        else
          CLWC      = [ vq; CLWC];
        end


  v       = CRWC;
  vq      = interp1( x, v, xq, 'linear', 'extrap');
        if ima > 0
          CRWC(ima) = vq;
          CRWC      = CRWC(ima:end); 
        else
          CRWC      = [ vq; CRWC];
        end

end


% last atmospheric layer to surface pressure  



% to volume mixing ratio in Kg/m3
RHO	   =  mixing_ratio_to_density(Q, T, P); 
CLW	   =  mixing_ratio_to_density(CLWC, T, P); 
CRW	   =  mixing_ratio_to_density(CRWC, T, P); 

cvar = [ RHO  CLW  CRW ];

% column integrations Kg/m3 * dz(m) -> Kg/m2
[cvar,  Z] = fom_aux_atmos_column_integration(cvar, P, T, Q);


TCWV   = cvar(:,1); 
TCLW   = cvar(:,2); 
TCRW   = cvar(:,3); 


return


%==========================================================================
% column integration
% IN    p         Column vector of pressures [Pa].
%       t         Column vector of temperatures [K].
%       q         Water vapour mixing ratio  [Kg/Kg]. Vector or a scalar, e.g. 0. 
%     var         npre x nproducts
%     cvar        integrated as var-unit * m e.g. Kg/m3 o Kg/m2 colum

function [cvar,z] = fom_aux_atmos_column_integration(var,p,t,q)

npro = size(var,2);


% approximated drop in pressure
% using a more elaborated one

% lp = log10(p)*16;
% z  = -lp;

% using hydrostatic equilibrium
% in very rare occasions it fails
% and we revert to a more robust
% but less accurate version

try
 
 z= pt2z(p,t,q,1e5,0);

catch

 Rd = 287; 
 Rv = 461;
 g0 = 9.78; 
 Tv = t .* (1 + (Rv*q/Rd)) ./ (1 + q);
 z = (Rd * Tv / g0) .* log(1e5 ./ p);

end  


if ~isempty( var )

  np   = length(p);
  cvar = 0;

  for p = 1:np-1

    if p == 1
      dz = (z(2)-z(1))/2;
    elseif p == np-1
      dz = (z(p)-z(p-1))/2;
    else
      dz = (z(p) - z(p-1))/2 + (z(p+1) - z(p))/2;  
    end 
   
    cvar = cvar + var(p,:) * dz;  

  end


else

  cvar = [];

end

return




%==========================================================================
% PT2Z   Hydrostatic altitudes
%
%    Calculates altitudes fulfilling hydrostatic equilibrium, based on
%    vertical profiles of pressure, temperature and water vapour. Pressure
%    and altitude of a reference point must be specified.
%
%    Molecular weights and gravitational constants are hard coded and 
%    function is only valid for the Earth.
%
%    As the gravitation changes with altitude, an iterative process is
%    needed. The accuracy can be controlled by *z_acc*. The calculations
%    are repeated until the max change of the altitudes is below *z_acc*. If
%    z_acc<0, the calculations are run twice, which should give an accuracy
%    better than 1 m.
%
% FORMAT   z = pt2z( p, t, h2o, p0, z0 [,lat,z_acc,refell] )
%        
% OUT   z         Altitudes [m].
% IN    p         Column vector of pressures [Pa].
%       t         Column vector of temperatures [K].
%       h2o       Water vapour [VMR]. Vector or a scalar, e.g. 0. 
%       p0        Pressure of reference point [Pa].
%       z0        Altitude of reference point [m].
%       lat       Latitude. Default is 45.
%       z_acc     Accuracy for z. Default is -1.
%       ellipsoid Reference ellipsoid data, see *ellipsoidmodels*. 
%                 Default is data matching WGS84.

% 2005-05-11   Created by Patrick Eriksson.


function z = pt2z(p,t,h2o,p0,z0)


%= Make rough estimate of *z*
%
z = p2z_simple( p );
z = shift2refpoint( p, z, p0, z0 );


%= Set Earth radius and g at z=0
%
re =  6378137;
g0 =  9.780327;


%= Gas constant and molecular weight of dry air and water vapour 
%
r  = 8.3145; %constants( 'GAS_CONST' );
md = 28.966;
mw = 18.016;
%
k  = 1-mw/md;        % 1 - eps         
rd = 1e3 * r / md;   % Gas constant for 1 kg dry air

np = length(p);

for iter = 1:2

  zold = z;
  
  g = z2g( re, g0, z );

  for i = 1 : (np-1)
      
	gp  = ( g(i) + g(i+1) ) / 2;
  
	%-- Calculate average water VMR (= average e/p)
	hm  = (h2o(i)+h2o(i+1)) / 2;
  
	%--  The virtual temperature (no liquid water)
	tv = (t(i)+t(i+1)) / ( 2 * (1-hm*k) );   % E.g. 3.16 in Wallace&Hobbs
  
	%-- The change in vertical altitude from i to i+1 
	dz = rd * (tv/gp) * log( p(i)/p(i+1) );
	z(i+1) = z(i) + dz;
      
  end
  
end

return




%==========================================================================

function z = shift2refpoint( p, z, p0, z0 )
  %
  z = z - ( interpp( p, z, p0 ) - z0 );
  %
return

%==========================================================================

function g = z2g(r_geoid,g0,z)
  %
  g = g0 * (r_geoid./(r_geoid+z)).^2;
  %
return

%==========================================================================

% NAME:    interpp
%
%          Makes log-linear interpolation. That is, the x-dimension is
%          converted to log before doing the interpolation.
%
%          The profiles are assumed to be constant outside the end points.
%
%          A typical application of the function is interpolation of
%          atmospheric vertical profiles.
%
% FORMAT:  X = interpp(pp,Xp,p)
%
% RETURN:  X          interpolated profiles
% IN:      pp         original pressure levels
%          Xp         original profiles
%          p          new pressure levels 
%------------------------------------------------------------------------

% HISTORY: 2005-05-11  Moved to Atmlab by PE.
%          2000-01-04  Moved from Norns to AMI
%          1999-11-02  Created by Patrick Eriksson.


function X = interpp(pp,Xp,p)

pp	= vec2col(pp);
np	= length(pp);

X	= interp1([1e3;log(pp);-1e3],[Xp(1,:);Xp;Xp(np,:)],log(p));

return

% VEC2COL   Ensures that a variable not has less rows than columnss.
%
%    The most common application of this function is to ensure that a 
%    vector is a column vector.
%
% FORMAT   v = vec2col(v)
%        
% OUT   v   A variable of any type.
% IN    v   The variable possible transposed.

% 1993        Created by Patrick Eriksson. 
% 2002-12-10  Adapted to Atmlab from arts/ami.
% 2013-03-04  Bugfix by Gerrit Holl for empty vector + do not get conj

function v = vec2col(v)

wid = ['atmlab:' mfilename ':zero'];
[rows,cols] = size(v);

if all(size(v)==0)
    warning(wid, 'Cannot convert zero-sized vector to column');
end

if (isempty(v) && ~iscolumn(v)) || (~isempty(v) && cols > rows)
    v = v.';
end

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
% OUT   z   Altitudes [m]
% IN   p   Pressure [Pa]

% 2005-05-11   Created by Patrick Eriksson.
%              from atmlab package http://www.radiativetransfer.org/tools/

function z = p2z_simple( p )

z = 16e3 * ( 5 - log10(p) );

return





%==========================================================================
% SUBMODULE mixing_ratio_to_density
%
%       Calculates water vapor density vd from mixing ratio
%
%          r  = density /density air
%
% FORMAT   vd = mixing_ratio_to_density(rh, T, p)
%        
% OUT   vd  density [kg/m3]
% IN    r   mixing ratio [kg/kg], it can be a scalar or a tensor
%       T   air temperature [K], it can be a scalar or a tensor
%       p   air pressure [Pa], it can be a scalar or a tensor
%
%-------------------------------------------------------------------------------
% Project:	  CIMR Algorithm Performance Evaluation
% Package:	  CIMR Scientific Work Bench
% Developer:	  Estellus 
% Contact:	  carlos.jimenez@estellus.fr 
% Initiated:	  2019-03-08
%-------------------------------------------------------------------------------

function vd = mixing_ratio_to_density(r, T, p)

% vd = r * air density
% set gass constant for dry air in J.kg-1.K-1
Rd = 287.0580; 


vd = r .* p ./ ( Rd * T);

return




