classdef CLP
% CLP provides routines to parse the command line arguments, generating
% a list of tokens grouped as:
% - Configuration file
% - Input files
% - Output files
% 
% Command line must follow the following format (expressed in Extended Backus-Naur form):
% 
%     <command_line> ::= <executable_name> <whitespaces> <configuration_files> <whitespaces>
%                        <input_files> <whitespaces> <output_files> <EOL>
%     <executable_name> ::= <file_name>
%     <whitespace>  ::= (“ “)
%     <whitespaces> ::= <whitespace>+
%     <file_name>   ::= (<alphanumeric>)+
%     <configuration_files> ::= <list_of_files>
%     <input_files>   ::= <list_of_files>
%     <output_files>  ::= <list_of_files>
%     <list_of_files> ::= <file_name> (“,” <file_name>)*
% 
% File names must be valid OS-dependant file locations. Per the E2E-ICD, the configuration, input
% and output file names can be either a full-path name or a relative path name. Relative path names
% are interpreted by the module as relative to the current working directory. In particular, no
% reference to the possible E2E_HOME variable is made.
%
% In order to increase the flexibility of the OpenSF/Modules CLI, the following new CLI is also available:
% 
% <command_line> ::= <executable> <options>
% <options> ::= (<global_configuration> | <local_configuration> | <input_file> | <output_file>)*
% <executable>> ::= <file_name>
% <global_configuration> ::= (--global|-g) <file_name>
% <local_configuration> ::= (--local|-l) <file_name>
% <input_file> ::= (--input|-i) <file_name>
% <output_file> ::= (--output|-o) <file_name>
% <file_name> ::= (<alphanumeric>)+

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

    properties (GetAccess = public, SetAccess = immutable)
        % These properties will be strings if the inputs are strings, and char cells otherwise
        confFiles; % Absence of GCF/LCF is represented by missing (for strings) or '' (for chars)
        inputFiles;
        outputFiles;
    end
    properties (Constant)
        delimiter = ',';
    end
    
    methods
        function obj = CLP(varargin)
            % CLP Create an instance of the CLP class.
            %   obj = CLP('g.xml,l.xml', '', 'out.dat') uses CLIv1 inputs
            %   obj = CLP('-g', "a.xml", '-o', "out.dat") uses CLIv2 inputs
            %
            % If succesful, the object properties (confFiles, inputFiles, outputFiles) will be
            % filled according to the inputs.
            %
            % Note that all inputs (or values, for v2) must be of the same type, that is, either
            % all strings or all character vectors. The properties of the resulting object will be
            % of the same type as the inputs/values passed.
            parsefun = @CLP.parseCLIv1Args;
            if CLP.looksLikeCLIv2(varargin{:})
                parsefun = @CLP.parseCLIv2Args;
            end
            [obj.confFiles, obj.inputFiles, obj.outputFiles] = parsefun(varargin{:});
        end
        function [conf_files]  = getConfFiles (obj)
            CLP.deprecated();
            conf_files = strjoin(obj.confFiles, ',');
        end
        function [conf_file]   = getConfFile (obj, index)
            CLP.deprecated();
            conf_file = CLP.getFile(obj.confFiles, index);
        end
        function [input_files] = getInputFiles (obj)
            CLP.deprecated();
            input_files = strjoin(obj.inputFiles, ',');
        end
        function [input_file]  = getInputFile (obj, index)
            CLP.deprecated();
            input_file = CLP.getFile(obj.inputFiles, index);
        end
        function [out_files]   = getOutputFiles (obj)
            CLP.deprecated();
            out_files = strjoin(obj.outputFiles, ',');
        end
        function [out_file]    = getOutputFile (obj, index)
            CLP.deprecated();
            out_file = CLP.getFile(obj.outputFiles, index);
        end
        
        function [nConf]       = nConfFiles (obj)
            CLP.deprecated();
            nConf = length(obj.confFiles);
        end
        function [nIn]         = nInputFiles (obj)
            CLP.deprecated();
            nIn = length(obj.inputFiles);
        end
        function [nOut]        = nOutputFiles (obj)
            CLP.deprecated();
            nOut = length(obj.outputFiles);
        end
    end

    methods(Static,Access=private)
        function deprecated() % Helper for the old API functions to throw a deprecation warning
            warning('OSFI:CLP:DeprecatedAPI', ...
                    'This function is deprecated and will be removed in the future - use the class properties instead');
        end
        function f = getFile(src, n)
            % Implementation for getXFile usable both with cellstr and string array properties
            if isstring(src)
                f = src(n);
            else
                f = src{n};
            end
        end

        function f = parseSingleCLIv1Arg(str)
            if strlength(str) == 0 % Fully empty argument -> empty string array or cellstr
                if isstring(str)
                    f = strings(0,1); 
                else
                    f = {};
                end
            else % strsplit keeps the type of the argument
                f = strsplit(str, CLP.delimiter, 'CollapseDelimiters', false);
                if isstring(str) % Replace "" with missing string
                    f = standardizeMissing(f, "");
                end
            end
        end
        function [conf, inp, out] = parseCLIv1Args(varargin)
            if ~all(cellfun(@isstring, varargin)) && ~all(cellfun(@ischar, varargin))
                error('OSFI:CLP:InvalidArgTypes', 'Inputs must be either all char vectors or all string scalars');
            end
            if nargin == 2
                args = [{''}, varargin]; % Legacy behaviour, for modules without config files
            elseif nargin == 3
                args = varargin;
            else
                error('OSFI:CLP:InvalidNargin', 'CLIv1: either 2 or 3 arguments expected, got %d', nargin);
            end
            conf = CLP.parseSingleCLIv1Arg(args{1});
            inp = CLP.parseSingleCLIv1Arg(args{2});
            out = CLP.parseSingleCLIv1Arg(args{3});
        end

        function [conf, inp, out] = parseCLIv2Args(varargin)
            if mod(nargin,2) ~= 0;
                error('OSFI:CLP:InvalidNargin', 'Unpaired key/value inputs');
            end
            % If all inputs are strings (not char vectors), use strings ourselves
            use_string_vals = all(cellfun(@isstring, varargin(2:2:end)));
            if ~use_string_vals && ~all(cellfun(@ischar, varargin(2:2:end)))
                error('OSFI:CLP:InvalidArgTypes', 'Values must be either all char vectors or all string scalars');
            end
            if use_string_vals
                conf = {missing, missing};
            else
                conf = {'', ''};
            end
            inp = {};
            out = {};
            for kv=reshape(varargin,2,[])
                if any(strcmpi(kv{1}, {'-g', '--global', 'global'}))
                    conf{1} = kv{2};
                elseif any(strcmpi(kv{1}, {'-l', '--local', 'local'}))
                    conf{2} = kv{2};
                elseif any(strcmpi(kv{1}, {'-i', '--input', 'input'}))
                    inp{end+1} = kv{2};
                elseif any(strcmpi(kv{1}, {'-o', '--output', 'output'}))
                    out{end+1} = kv{2};
                else
                    error('OSFI:CLP:InvalidKey', 'Unrecognized key %s. Must be one of "global", "local", "input" or "output"', kv{1});
                end
            end
            if use_string_vals
                conf = string([conf{:}]);
                inp = string([inp{:}]);
                out = string([out{:}]);
            end
        end

        function res = looksLikeCLIv2(varargin)
            keys = varargin(1:2:end);
            if ~isempty(keys) && strncmp(keys{1}, '-', 1)
                % Definitely yes, 1st arg starts with "-" so parse as CLPv2
                res = true;
                return
            end
            res = false;
            % TODO: Recognition of keywords for a function call like module('Global', 'X', 'Input', 'Y.in') ?
            % if mod(nargin, 2) ~= 0 && ...
        end
    end
    
end


