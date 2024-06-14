%-------------------------------------------------------------------------------
%
% SUBMODULE   npi2pi
%
%    Cobeting lon 0 360 to lon -180 180
%
% FORMAT  lon = npi2pi(lon) 
%        
% OUT   lon   converted longitude
% IN    lon   input longitude
%
%-------------------------------------------------------------------------------
% Project:	  CIMR Algorithm Performance Evaluation
% Package:	  CIMR Scientific Work Bench
% Developer:	  Estellus 
% Contact:	  carlos.jimenez@estellus.fr 
% Initiated:	  2019-11-30
%-------------------------------------------------------------------------------

function lon = npi2pi(lon)

ind = find( lon > 180 );
lon(ind) = lon(ind)-360;

return
