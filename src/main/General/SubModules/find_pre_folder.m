%-------------------------------------------------------------------------------
%
% SUBMODULE   find_pre_folder
%    For a given path, it finds the folder one step above,
%    or if the path is a file, the folder where the file is.
%    It also outputs the remaining part of the path.
%    
% FORMAT   dirout = find_pre_folder( dirin )
%        
% OUT   dirout	Directory path, or directory and file name
%	remout	Reamining part in dirin 
% IN    dirin	Directory path as explained above.
%
%-------------------------------------------------------------------------------
% Project:	  CIMR Algorithm Performance Evaluation
% Package:	  CIMR Scientific Work Bench
% Developer:	  Estellus 
% Contact:	  carlos.jimenez@estellus.fr 
% Initiated:	  2019-01-02
%-------------------------------------------------------------------------------

function [ dirout, remout ] = find_pre_folder( dirin )

aux    = strfind( dirin, '/' );
dirout = dirin(1:aux(end)-1);
remout = dirin(aux(end)+1:end);

return
