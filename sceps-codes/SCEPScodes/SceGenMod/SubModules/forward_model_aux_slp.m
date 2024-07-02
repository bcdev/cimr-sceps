%-------------------------------------------------------------------------------
%
% SUBMODULE  forward_model_aux_slp 
%
%  Adjusting atmopsheric profiles to make them consistent with 
%  surface pressure. If surpre > p(1) we just extrapolate to
%  the P(1) atmospheric quantity value. Otherwise, we linearly
%  interpolate with log10(P)
%
% FORMAT  atmos_input = forward_model_aux_slp( atmos_input, surpre ) 
%        
% OUT   atmos_input   structure with changed atmospheric data
%
% IN    atmos_input		structure with original atmospheric data;
%				see forward_model_core.m
%	sur_pre	      [mbar]	surface pressure 
%
%-------------------------------------------------------------------------------
% Project:	  CIMR Algorithm Performance Evaluation
% Package:	  CIMR Scientific Work Bench
% Developer:	  Estellus 
% Contact:	  carlos.jimenez@estellus.fr 
% Initiated:	  2020-05-30
%-------------------------------------------------------------------------------

function atmos_input = forward_model_aux_slp( atmos_input, surpre ) 


%= first cutting profiles if needed

x       = log10 (atmos_input.P);
xq      = log10( surpre );

xa      = x - xq;
ind     = find( (xq-x) < 0 );
if isempty(ind)
  ima     = 1;
else
  ima     = ind(end);
end

x       = log10 (atmos_input.P);
xq      = log10( surpre );

v       = atmos_input.T;
vq      = interp1( x, v, xq, 'linear', nan );

        if ~isnan(vq)
          atmos_input.T(ima) = vq;
        end

v       = atmos_input.RHO;
vq      = interp1( x, v, xq, 'linear', nan );
        if ~isnan(vq)
          atmos_input.RHO(ima) = vq;
        end

v       = atmos_input.LW;
vq      = interp1( x, v, xq, 'linear', nan );
        if ~isnan(vq)
          atmos_input.LW(ima) = vq;
        end

% last atmospheric layer to surface pressure  


atmos_input.P(ima) = surpre; 
atmos_input.P      = atmos_input.P(ima:end); 
atmos_input.T      = atmos_input.T(ima:end); 
atmos_input.RHO    = atmos_input.RHO(ima:end); 
atmos_input.LW     = atmos_input.LW(ima:end); 


if isfield( atmos_input, 'CWVC' ) | isfield( atmos_input, 'CLWC' )

    H    = p2z_barometric( atmos_input.P, atmos_input.T, 0 );
    NLEV = length(atmos_input.P);

end


if isfield( atmos_input, 'CWVC' )

    RHO = atmos_input.RHO;

    %= column integrating RHO

    CRHO = RHO(1) * 0.5 * (H(2) - H(1)) ;

    for p = 2:NLEV-1

      hu = H(p+1) - (H(p+1)-H(p))/2;
      hl = H(p-1) + (H(p)-H(p-1))/2;

      CRHO = CRHO + RHO(p) * ( hu - hl );  

    end

    CRHO = CRHO + RHO(NLEV) * 0.5 * (H(NLEV) - H(NLEV-1));   %g/m2
    atmos_input.CWVC = CRHO;

end


  
%=== adjustment of LW profile to CLWC

if isfield( atmos_input, 'CLWC' )

    CLD = atmos_input.LW;

    %= column integrating LW

    CLW = CLD(1) * 0.5 * (H(2) - H(1)) ;

    for p = 2:NLEV-1

      hu = H(p+1) - (H(p+1)-H(p))/2;
      hl = H(p-1) + (H(p)-H(p-1))/2;

      CLW = CLW + CLD(p) * ( hu - hl );  

    end

    CLW = CLW + CLD(NLEV) * 0.5 * (H(NLEV) - H(NLEV-1)); 
    atmos_input.CLWC = CLW;

end



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
R    = 8.3144598;   % universal gas constant: J/(mol·K)
g    = 9.80665;     % gravitational acceleration: m/s2

% assuming constant gravitational acceleration and molar mass
z    = R.*(t-t(1))./Mdry./g.*log(p./p(1))./log(t(1)./t); 

% reference altitude
z(1) = zo;

return

