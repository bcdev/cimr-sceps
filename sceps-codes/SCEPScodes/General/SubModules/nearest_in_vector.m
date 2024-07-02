% SUB MODULE nearest_in_vector
%
%          Find the nearest value in a vector z to a given value zi
%
% FORMAT   [ind,zn,zd] = nearest_in_vec(z,zi)
%
% OUT      ind  the index in vector z of closest value to zi
%          zn   the closest value in z to zi
%          zd   the distance between zn and zi
% IN       z    the vector
%          zi   the value

% 2019-01-08   Created by carlos.jimenez@obspm.fr

function [ind,zn,zd] = nearest_in_vector(z,zi);

[ aux, ind ] = min( abs( z -zi ) );
zn           = z( ind );
zd           = zi - zn; 

return





