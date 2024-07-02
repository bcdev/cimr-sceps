%-------------------------------------------------------------------------------
%
% SESSION   An (almost) empty session header.
%
%    XXX
%
% FORMAT   xxx
%        
% OUT   x   xxx
% IN    x   xxx
% OPT   x   xxx
%
%-------------------------------------------------------------------------------
% Project:	  CIMR Algorithm Performance Evaluation
% Package:	  CIMR Scientific Work Bench
% Developer:	  Estellus 
% Contact:	  carlos.jimenez@estellus.fr 
% Initiated:	  YYYY-MM-DD
%-------------------------------------------------------------------------------


global E2E_HOME
global SOF_CIMR_PATH
global LOG

session_id  = '';

LOG.info( ['Starting session ', session_id ] );


%= changing to session folder 

sfo = [ E2E_HOME, '/sessions/', session_id ];

if ~exist( sfo, 'dir' )

  LOG.info( [ 'Creating session folder ', sfo ]);
  mkdir( sfo);

end  

cd( [ E2E_HOME, '/sessions/', session_id ] );



%= copying global configuration file to session folder

[success] = copyfile( [ SOF_CIMR_PATH, '/ConfigFiles/Global_Configuration.xml'] );

if success == 0
  LOG.error('Problem with session folder, does it exist?')
  error('');
end



%= copying local configuration files to session folder
%  as there is only one module in the simulation we copy
%  just one

[success] = copyfile( [ SOF_CIMR_PATH, '/ConfigFiles/', session_id, '/', lower( session_id ), '_Local_Configuration.xml'] );

if success == 0
  LOG.error('Problem with local configuration file, does it exist?')
  error('');
end



%= starting a simulation




LOG.info( ['Finishing session ', session_id ] );

