% SUBMODULE   save_data_single
%
%    Simple function to save the data passed in a matlab
%    file with data stored in a variable named 'data'.
%    
% FORMAT   save_data_single( outfile, data )
%        
% OUT   -
% IN    outfile	string		File name, full path
%	data	[as stored]	Variable to save
%
%-------------------------------------------------------------------------------
% Project:	  CIMR Algorithm Performance Evaluation
% Package:	  CIMR Scientific Work Bench
% Developer:	  Estellus 
% Contact:	  carlos.jimenez@estellus.fr 
% Initiated:	  2019-01-02
%-------------------------------------------------------------------------------

function save_data_single( outfile, data )

save( outfile, 'data', '-v7.3' );

return
