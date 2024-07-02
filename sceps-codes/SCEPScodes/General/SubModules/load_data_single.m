%-------------------------------------------------------------------------------
% SUBMODULE   load_data_single
%
%    Simple function to load the data stored in a matlab
%    file with data stored in a variable named 'data'.
%    
% FORMAT   load_data_single( infile )
%        
% OUT   data	[as stored]	Stored variable
% IN    infile	string		File name, full path
%
%-------------------------------------------------------------------------------
% Project:	  CIMR Algorithm Performance Evaluation
% Package:	  CIMR Scientific Work Bench
% Developer:	  Estellus 
% Contact:	  carlos.jimenez@estellus.fr 
% Initiated:	  2019-01-04
%-------------------------------------------------------------------------------

function data  = load_data_single( infile )

data = load( infile, 'data' );
data = data.data;

return
