%-------------------------------------------------------------------------------
%
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

% set gass constant for dry air in J.kg-1.K-1
Rd = 287.0580; 

vd = r .* p ./ ( Rd * T);  

return
