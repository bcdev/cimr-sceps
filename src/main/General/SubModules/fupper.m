%-------------------------------------------------------------------------------
% SUBMODULE   fupper
%
%    Simple function that works like upper but only
%    to the first letter.
%    
% FORMAT   wordout = fupper( wordin )
%        
% OUT   wordout    string	   String converted
% IN    wordin     string	   String to be converted
%
%-------------------------------------------------------------------------------
% Project:	  CIMR Algorithm Performance Evaluation
% Package:	  CIMR Scientific Work Bench
% Developer:	  Estellus 
% Contact:	  carlos.jimenez@estellus.fr 
% Initiated:	  2019-01-02
%-------------------------------------------------------------------------------

function wordout = fupper( wordin )

wordout = [ upper(wordin(1)), wordin(2:end) ];

return
