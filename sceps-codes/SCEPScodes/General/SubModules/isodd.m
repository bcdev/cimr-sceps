% ISODD   True for odd integers.
%
%    Determines if all the elements of an array are odd integers. 
%    
%    Note that this function differs from *isinteger*, which is a type 
%    test. Here also double numbers are accapted, as long as they are
%    odd integers.
%
% FORMAT   bool = iseven(n)
%        
% OUT   bool   true if odd integer, else false.
% IN    n      A variable of any valid numeric data type.

% 2002-12-10  Created by Patrick Eriksson. 


function bool = isodd(n)

bool = false;

if isempty(n) | ~isnumeric(n), return, end

if any( rem(n(:),2) == 0 ), return, end

bool = true;
