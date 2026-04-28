classdef Logger < handle
% Logger module provides routines for presenting messages in the standard output.
% Any E2E-ICD-complieant framework can later intercept these messages.
% 
% Logger formats the logging messages as follows:
% 
%     <message>   ::= (<progress> | <log>) <EOL>
%     <progress>  ::= "Progress" <whitespaces> <delimiter> <whitespace> <progress_body>
%     <delimiter> ::= "|"
%     <progress_body> ::= <integer> " of " <integer>
%     <log> ::= <type> <whitespaces> <delimiter> <whitespaces> <text> [<whitespaces>
%               <delimiter> <whitespaces> <version>]
%     <type>    ::= "Error" | "Warning" | "Info" | "Debug"
%     <version> ::= <digit>("." <digit>)*
%     <digit>   ::= "0" | "1" | "2" | "3" | "4" | "5" | "6" | "7" | "8" | "9"
%     <whitespace>  ::= (" ")
%     <whitespaces> ::= <whitespace>+
% 
% This format defines five different types of messages:
% - Information. This is a message raised by the model describing a harmless event,
%     giving knowledge to the user. Model execution should continue with no interruptions.
% - Warning. The model has detected a non-fatal error or situation that may cause a
%     fatal error. This is a harmless event, thus, the execution should continue with no
%     interruption.
% - Error. A fatal error has happened in the model execution and the model itself
%     informs the user about it, so the model has time to gracefully close the execution.
% - Debug. Detailed information of the model execution given to the user. Information
%     is intended to lead the user (or model developer) upon fixing a problem. This is a
%     harmless event so model execution should continue with no interruptions. Debug messages
%     are only shown if an environment variable named DEBUG_MODE is defined and set as On in
%     the model execution context.
% - Progress. Numerical information on the amount of model execution performed.
% 
% Logger also provides a way to finish the model execution (after a fatal error or at
% the expected end of the execution) and return a non-zero code to the operating system.

%
% Copyright notice - must go below the doc and divided by a non-comment line to avoid replacing the "help" functionality
% openSF Integration Libraries (OSFI)
% Deimos Space, S.L.U.
% 
% This file is part of OSFI. OSFI is free software; you can redistribute it
% and/or modify it under the terms of the 'ESA Software Community Licence Permissive' as
% published by the European Space Agency; either version 2.4 of the License,
% or (at your option) any later version. You should have received a
% copy of the 'ESA Software Community Licence Permissive - v2.4' along with this program
% or one can be found at <http://eop-cfi.esa.int/index.php/docs-and-mission-data/licensing-documents>.
%

    properties (Access = private)
        debugging      = false;
        delimiter      = '|';
        logFileName    = '.tmpLogFile';
        standAloneMode = true;
    end
    
    methods
        function Logger = Logger (standAlone)
            if (nargin>0)
                warning('OSFI:EHLogger:DeprecatedAPI', 'Calling Logger with nargin > 0 is deprecated and will be removed in the future.');
            end
            
            Debug_Mode = getenv('DEBUG_MODE');
	        if (strcmpi(Debug_Mode,'On')) 
            	Logger.debugging = true; 
	        else
	    	    Logger.debugging = false;
	        end
        end
        function error(LOG, message)
            LOG.write('Error', message);
        end
        function warning(LOG, message)
            LOG.write('Warning', message);
        end
        function info(LOG, message)
            LOG.write('Info', message);
        end
        function debug(LOG, message)
            if (LOG.debugging)
                LOG.write('Debug', message);
            end
        end
        function progress(LOG, step, nSteps)
            LOG.write('Progress', sprintf('%d of %d', step, nSteps));
        end
        function finishExecution(LOG, exitCode)
            LOG.info('Finishing module execution');
            if nargin < 2 % previous behaviour - just raise an error
                error('Finishing module execution');
            else % Closer alignment with other OSFIs: terminate
                exit(exitCode)
            end
        end
        function qualityReport(LOG, name, value)
           if ischar(value)
                LOG.write('Quality', sprintf('%s: %s', name, value));
           else
                LOG.write('Quality', sprintf('%s: %f', name, value));
           end
        end
        
        function setDebugMode (LOG, debugging)
            LOG.debugging = debugging;
        end
        function setStandAlonMode (LOG, standAlone)
            LOG.standAloneMode = standAlone;
        end
        
        
    end
    
    methods (Access = private)   
        function putLog(LOG, message)
            if (LOG.standAloneMode)
                disp (message);
            else
                log_file = fopen(LOG.logFileName,'a+');
                fprintf(log_file, message);
                fclose(log_file);
            end
        end 
        
        function write(LOG, type, message)
            log = sprintf('%-8s %s %s', type, LOG.delimiter, message);
            LOG.putLog (log);
        end
        
        function version = getVersion(LOG) 
		   persistent osfiVersion;
           if isempty(osfiVersion)
                try
                    osfiVersion = osfi_version();
                catch
                    osfiVersion = '(unknown)';
                    LOG.warning('No osfi_version function found in the path! Running OSFI from the development tree?')
                end
           end 
		   version =  osfiVersion;      
       end
                              
    end
    
end

