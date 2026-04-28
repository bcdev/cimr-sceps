%=================================================================
%
% Mscript to initilaize the use of the SCEPSscd matlab simulations
%
% FORMAT   startup_matlab_SCEPScodes( devSCEPSpath, dataSCEPSpath )
%        
% OUT   
%
% IN    codeSCEPSpath	Absolute path to SCEPScodes folder	
%	dataSCEPSpath	Absolute path to SCEPSdata folder		
%	wspaceSCEPSpath	Absolute path to SCEP workspace folder		


function startup_matlab_SCEPScodes( codeSCEPSpath, dataSCEPSpath, wspaceSCEPSpath  ) 


%= in case any parpool remains opened

delete(gcp('nocreate'))

%= who am i


global UID
[s,UID] = system('echo $HOME');
UID = UID( 1 : (length(UID)-1) );


if nargin == 0
  codeSCEPSpath  = [ UID, '/Work/Studies/SCEPS2OpenSF/SCEPSCodes' ];
  dataSCEPSpath    = [ UID, '/Work/DataE/SCEPS2OpenSF/SCEPSData' ];
  wspaceSCEPSpath  = [ UID, '/Work/DataE/SCEPS2OpenSF/MatlabWorkSpace' ];
end


%= saving lobal variables for both matlab and opensf paths

global SCEPS_CODES_PATH
SCEPS_CODES_PATH = codeSCEPSpath;





%= adding OSFI library for Matlab

addpath([ codeSCEPSpath, '/OSFI' ]);



%= setting log class and making it global 

LOG = Logger();

disp(' ');
LOG.info(['Starting Software Workbench for an End-to-End simulation chain for CIMR']);




%= setting OSF home variable to the data package location
%  as in OpenSF, and making system environmental variable
%  to mimic OpenSF

E2E_HOME = [ UID, '/Work/DataE/SCEPS2OpenSF/SCEPSData' ];
setenv( 'E2E_HOME', E2E_HOME )



%= adding path to general folder with modules, simulations, and sessions
%  user can define new paths to use other folders

% General folder from OpenSF-github-based codes

addpath([ codeSCEPSpath, '/General/ConfigFiles']);
addpath([ codeSCEPSpath, '/General/SubModules']);



% Scene Generation Module

addpath([ codeSCEPSpath, '/SceGenMod/ConfigFiles']);
addpath([ codeSCEPSpath, '/SceGenMod/Modules']);
addpath([ codeSCEPSpath, '/SceGenMod/SubModules']);

% this is a link from MwInversionSuite where we are storing forward models
% when distributed we need to copy this

addpath([ codeSCEPSpath, '/SceGenMod/SubModules/ForwardModels']);

addpath([ codeSCEPSpath, '/SceGenMod/Simulations']);


% Observing System  Module

addpath([ codeSCEPSpath, '/ObsSimMod/ConfigFiles']);
addpath([ codeSCEPSpath, '/ObsSimMod/Modules']);

addpath([ codeSCEPSpath, '/ObsSimMod/Simulations']);



%= adding path to general folder with sessions
%  user can define a new path to use other folders
%  we just recopy the general one here
%  But better to remove so the user needs to 
%  think about this and go to his session

addpath([ codeSCEPSpath, '/Sessions']);


%= adding path to compilations folder

addpath([ codeSCEPSpath, '/RunTimeEnv']);


%= simulation folders for each session will be 
%  saved in a WorkSpace as OpenSF, NOT in sessions
%  as before, /sessions now is just for the 
%  mscripts simulating an OpenSF run
%  We make a global variable to pass this WorkSPace

global SCEPS_WORK_SPACE
SCEPS_WORK_SPACE = wspaceSCEPSpath;


return


