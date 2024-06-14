%-------------------------------------------------------------------------------
%
% SUBMODULE   osfi_error
%
%    Simple function to report an error by using the OSFI
%    class LOG and halting or not the matlab execution. It
%    requires that a global og class named LOG has been 
%    previoulsy created for the current matlab session, 
%    
% FORMAT   osfi_error( msgd, do_halt )
%        
% OUT   -
% IN    msgd     string	   Message to be displayed
%	do_halt  boolean   If 1(0) it (does/does not) halt the current execution
%
%-------------------------------------------------------------------------------
% Project:	  CIMR Algorithm Performance Evaluation
% Package:	  CIMR Scientific Work Bench
% Developer:	  Estellus 
% Contact:	  carlos.jimenez@estellus.fr 
% Initiated:	  2019-01-02
%-------------------------------------------------------------------------------


function osfi_error( msgd, do_halt )

if nargin == 1
  do_halt = 1;
end

global LOG

LOG.error ( msgd );
if do_halt
  error( msgd );
end

return
