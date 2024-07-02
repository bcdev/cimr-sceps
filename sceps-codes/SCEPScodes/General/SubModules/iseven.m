% ISEVEN   True for even integers.
%
%    Determines if all the elements of an array are even integers. 
%    
%    Note that this function differs from *isinteger*, which is a type 
%    test. Here also double numbers are accapted, as long as they are
%    even integers.
%
% FORMAT   bool = iseven(n)
%        
% OUT   bool   true if even integer, else false.
% IN    n      A variable of any valid numeric data type.

% 2002-12-10  Created by Patrick Eriksson. 


function bool = iseven(n)

bool = false;

if isempty(n) | ~isnumeric(n), return, end

if any( rem(n(:),2) == 1 ), return, end

bool = true;
