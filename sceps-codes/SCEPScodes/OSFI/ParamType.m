classdef ParamType
%
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
    enumeration
        INTEGER, FLOAT, BOOLEAN, STRING, FILE, FOLDER, TIME
    end

    methods
        function val = parseValues(obj, txt)
            % PARSEVALUES receives either a string or a cell array of strings
            % and returns the parsed value according to the parameter type.
            % If given a cell array, the output will be an array or cell, with
            % the same shape as the input.
            switch obj % Each branch must define val _and_ parseErrors (or return on its own)
            case {ParamType.INTEGER, ParamType.FLOAT}
                val = str2double(txt);
                parseErrors = isnan(val) & ~strcmpi(txt, 'nan');
                if obj == ParamType.INTEGER
                    % The range mandated in the E2E-ICD is that of Java int = Matlab int32. Type double
                    % has enough precision to represent all int32 values exactly, so if any fpVal ~= val
                    % (conversion doesn't roundtrip), it was not a valid input. This catches both 
                    % out-of-range integers and floating point values.
                    fpVal = val;
                    val = int32(val);
                    parseErrors = parseErrors | fpVal ~=val;
                end

            case ParamType.BOOLEAN % NOTE: per the E2E-ICD, only capitals!
                val = strcmp(txt, 'TRUE');
                % Check that the false values in val are actually FALSE
                shouldBeNotVal = strcmp(txt, 'FALSE');
                parseErrors = (val == shouldBeNotVal); % Fields that are neither TRUE nor FALSE
            
            case ParamType.TIME
                if iscellstr(txt)
                    % Parse into a cell, storing empty values in any unparseable inputs
                    val = cellfun(@(s) TimeValue.parse(s), txt, 'UniformOutput', false, ...
                            'ErrorHandler', @(id, msg, idx) []);
                    parseErrors = cellfun(@isempty, val);
                    if ~any(parseErrors(:)) % Recast into a TimeValue object array
                        val = reshape([val{:}], size(val));
                    end
                else
                    val = TimeValue.parse(txt);
                    return % Error checking performed by the function itself
                end

            case ParamType.STRING % No-op
                val = txt;
                return

            case {ParamType.FILE, ParamType.FOLDER}
                val = ParamType.getFileFormattedValues(txt);
                return % Error checking is performed by the function itself

            otherwise
                error('Unknown parameter type %s', char(obj));
            end

            %% Check for errors and report any with the offending values
            if any(parseErrors(:)) % Logical matrix with the same dimensions as txt/val
                error('Invalid value(s) for %s: %s', char(obj), ...
                        ParamType.parseErrorBadValues(txt, parseErrors));
            end
        end

        function split = splitValues(obj, txt, expect_empty)
            % SPLITVALUES splits a string depending on the parameter type. Numeric types
            % and BOOLEAN are split on whitespace, while the string-like types expect a
            % set of quoted values separated by whitespace. The output is always a cell
            % array of the tokens.
            txt = strtrim(txt);
            if isempty(txt) && nargin > 2 && expect_empty
                split = {}; % For empty nodes in ARRAY-typed parameters, since regex split would return {''}
            elseif obj.nonScalarsAreQuoted()
                % In this case we don't want to collapse delimiters, since '' is a valid element
                split = strsplit(txt, '''\s+''', 'CollapseDelimiters', false, 'DelimiterType', 'RegularExpression');
                if ~isempty(split)
                    % Since we have split on "' '", we need to remove the begin-quote from
                    % the first item and the end-quote from the last one.
                    split{1}   = regexprep(split{1},   '^''', '');
                    split{end} = regexprep(split{end}, '''$', '');
                end
            else
                % Split on whitespace, collapse multiple delimiters
                split = strsplit(txt);
            end
        end

        function isquoted = nonScalarsAreQuoted(obj)
            switch obj
            case {ParamType.INTEGER, ParamType.FLOAT, ParamType.BOOLEAN, ParamType.TIME}
                isquoted = false;

            case {ParamType.STRING, ParamType.FILE, ParamType.FOLDER}
                isquoted = true;
            otherwise
                error('Unknown parameter type %s', char(obj));
            end
        end
    end

    methods (Static=true)
        function str = parseErrorBadValues(txt, badMask)
            if ~iscellstr(txt) % Single value
                str = ['"' txt '"'];
                return
            end
            % Non scalar: build a string with (indices)="value" with all the invalid values
            badVals = txt(badMask);
            badLinIdx = find(badMask);
            dims = size(badMask);
            if isvector(badMask) % Avoid matrix indices for vectors
                dims = length(badMask); 
            end
            formattedBadVals = cell(1,numel(badVals));
            for iBad = 1:numel(badVals)
                subs = cell(1,length(dims));
                [subs{:}] = ind2sub(dims, badLinIdx(iBad));
                formattedBadVals{iBad} = sprintf('(%s)="%s"', ...
                        strjoin(cellfun(@num2str, subs, 'UniformOutput', false), ','), badVals{iBad});
            end
            str = [ 'at ' strjoin(formattedBadVals, ', ')];
        end

        function ffv = getFileFormattedValues(val)
            ffv = val;
            if ispc % Windows accepts both separators so normalize to \
                ffv = strrep(ffv, '/', '\');
            end
            
            % If ffv is a string, this function converts from string to char in the input (to then call cellfun)
            inputStr = isstring(val);
            if inputStr
                val = convertStringsToChars(val);
            end

            % Check if any of the path(s) is not absolute (no such function in Matlab)
            is_abs_fun = @(f) logical(javaMethod('isAbsolute', javaObject('java.io.File', f)));
            if iscell(val)
                isabs = cellfun(is_abs_fun, val);
            else
                isabs = is_abs_fun(val);
            end
            if all(isabs)
                return % XXX we might want to perform separator conversion here?
            end
            % If any relative path, apply the base dir from environment variables, to them (only!)
            base = ParamType.getFilesBaseDir(OsfiEnv.getInstance);
            if ~isempty(base)
                base_path = javaMethod('toPath', javaObject('java.io.File', base));
                make_abs = @(f) char(base_path.resolve(f));
                if iscell(val)
                    ffv(~isabs) = cellfun(make_abs, ffv(~isabs), 'UniformOutput', false);
                elseif inputStr
                    % If ffv is a string, convert from char to string
                    ffv = string(make_abs(val));
                else
                    ffv = make_abs(val);
                end
            end
        end
        
        function base = getFilesBaseDir(env)
            % GETFILESBASEDIR returns the "base directory" used to resolve
            % relative paths of FILE and FOLDER types. Per the the E2E-ICD,
            % it is read from the enviroment variable E2E_HOME. As a
            % deprecated option, if that variable is not present, the
            % OSFI_HOME variable is consulted next. If both are absent from
            % the environment, then the current working directory is used.
            ENVKEY_BASE_DIR = 'Parameter.base_dir';
            if nargin == 0
                env = OsfiEnv.getInstance;
            end
            storedVal = env(ENVKEY_BASE_DIR);
            if iscellstr(storedVal) && length(storedVal) == 1
                base = storedVal{1};
                return
            end
            
            e2e_home = env.getExternalVar('E2E_HOME');
            legacy_home = env.getExternalVar('OSFI_HOME');
            if ~isempty(e2e_home)
                base = e2e_home;
            elseif ~isempty(legacy_home)
                warning('OSFI:ConFM:Deprecated', ['E2E_HOME is not set, but OSFI_HOME is. '...
                        'Note that OSFI_HOME is deprecated and may be removed.']);
                base = legacy_home;
            else
                warning('OSFI:ConFM:HOME_UNDEF', ...
                        'E2E_HOME is not set, using the current folder for relative paths.');
                base = '';
            end
            env(ENVKEY_BASE_DIR) = {base}; % Prevent '' from being swallowed
        end
    end
end
