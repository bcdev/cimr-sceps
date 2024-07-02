%-------------------------------------------------------------------------------
%
% SUBMODULE   Snow depth from AMSR observations, parameterized
%	      from Killic et al., 2019
%	      https://doi.org/10.5194/tc-13-1283-2019
%
% FORMAT   snd = snow_depth_from_AMSR( tb6v, tb18v, tb36v)
%        
% OUT   snd	Snow depth in meters
% IN    tb6v	6.9 GHz V-pol brightness temperatures
%	tb18v	18.7 GHz V-pol brightness temperatures
%	tb36v	36.9 GHz V-pol brightness temperatures
%
%-------------------------------------------------------------------------------
% Project:	  CIMR Algorithm Performance Evaluation
% Package:	  CIMR Scientific Work Bench
% Developer:	  Estellus 
% Contact:	  carlos.jimenez@estellus.fr 
% Initiated:	  2019-11-04
%-------------------------------------------------------------------------------

function snd = snow_depth_from_AMSR( tb6v, tb18v, tb36v)

snd    = 1.7701 + 0.0175 * tb6v - 0.0280 * tb18v + 0.0041 * tb36v;

snd( snd < 0 ) = 0; 

return
