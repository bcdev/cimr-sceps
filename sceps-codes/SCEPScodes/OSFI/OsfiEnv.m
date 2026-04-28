classdef (Sealed=true) OsfiEnv < handle
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
properties (Access=private)
        mockedEnv % Empty in normal operation, otherwise it replaces getenv calls
        settings  % Stores the internal configuration data from OSFI classes
    end

    methods (Static)
        function obj = getInstance()
            % GETINSTANCE retrieves the current instance of the singleton.
            obj = OsfiEnv.instance();
        end
    end

    methods
        function var = getExternalVar(self, name)
            % GETEXTERNALVAR retrieves a value from the external
            % configuration, which is normally the environment variables.
            if ~isempty(self.mockedEnv)
                if self.mockedEnv.isKey(name)
                    var = self.mockedEnv(name);
                else
                    var = [];
                end
            else
                var = getenv(name);
            end
        end

        %% Delegate assignments to the settings map
        function ret = subsref(self, S)
            if ~isscalar(S) || ~strcmp(S.type, '()')
                ret = builtin('subsref', self, S);
            elseif ~iscellstr(S.subs) || length(S.subs) ~= 1
                error('OSFI:BadArgument', 'The settings keys must be a single string');
            else
                key = S.subs{1};
                if self.settings.isKey(key)
                    ret = self.settings(key);
                else % Avoid throwing like Map does, just return an empty value
                    ret = [];
                end
            end
        end

        function ret = subsasgn(self, S, val)
            if ~isscalar(S) || ~strcmp(S.type, '()')
                ret = builtin('subsasgn', self, S, val);
            elseif ~iscellstr(S.subs) || length(S.subs) ~= 1
                error('OSFI:BadArgument', 'The settings keys must be a single string');
            else
                key = S.subs{1};
                if ~isempty(val)
                    self.settings(key) = val;
                elseif isa(val,'double') && self.settings.isKey(key)
                    % We got a literal [] (or as close as we can detect), so delete the key
                    self.settings.remove(key);
                end
                ret = self;
            end
        end
    end

    %% Functions only accessible internally and by the UTs
    methods (Static, Access=?EnvMocker)
        function ret = instance(replacement, onlyIfCurrent)
            % INSTANCE stores the persistent shared instance.
            %
            % env = INSTANCE() retrieves the current singleton
            % env = INSTANCE(new) replaces the current singleton by new and returns the prior value
            % env = INSTANCE(new,cur) works as the previous overload if and only if the current
            %       singleton is "cur". Otherwise, it is not replaced and [] is returned.
            persistent singleton
            if isempty(singleton)
                singleton = OsfiEnv;
            end

            if nargin > 1 && singleton ~= onlyIfCurrent
                % Somebody else replaced the current env, reject change and return []
                ret = [];
                return
            end

            ret = singleton;
            if nargin > 0 % Replace the current value, for UTs
                singleton = replacement;
            end
        end
    end

    methods(Access=?EnvMocker)
        function self = OsfiEnv(varargin)
            % OSFIENV creates an instance of the settings storage used by OSFI.
            %
            % env = OSFIENV() creates a normal environment that delegates to getenv()
            % env = OSFIENV('key', 'value'...) creates an enviroment with mock external vars

            self.settings = containers.Map('KeyType', 'char', 'ValueType', 'any');
            if nargin == 0
                self.mockedEnv = [];
            else % Create an object that will not consult getenv, but rather the given "env vars"
                mockedEnv = convertContainedStringsToChars(varargin);
                if ~iscellstr(mockedEnv) || mod(numel(mockedEnv),2) ~= 0
                    error('OSFI:BadArgument', 'arguments must be key-value pairs');
                elseif ~isempty(mockedEnv)
                    self.mockedEnv = containers.Map(mockedEnv(1:2:end), mockedEnv(2:2:end));
                else % Fully empty mocked env
                    self.mockedEnv = containers.Map('KeyType', 'char', 'ValueType', 'char');
                end
            end
        end
    end
end
