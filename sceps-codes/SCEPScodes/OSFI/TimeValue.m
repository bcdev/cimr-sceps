classdef (Sealed=true) TimeValue
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
    properties (SetAccess = private)
        year    %int16
        month   %int8
        dom     %int8
        hour    %int8
        minute  %int8
        sec     %int8
        nanosec %int32
    end
    % Regexes for parse. Note that unlike other languages, Matlab does not accept nested groups,
    % so the optional fraction-of-second part cannot be (?:\.(\d+))? capturing the digits. Instead,
    % it must be (\.\d+)?, capturing also the period.
    properties (Access=private, Constant)
        % CCSDS time code A format, as YYYY-MM-DD"T"HH:mm:SS.f+"Z"
        TIME_CCSDS_ASCII_A = '^(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2})(\.\d+)?Z?$'
        % CCSDS time code B format, as YYYY-ddd"T"HH:mm:SS.f+"Z"
        TIME_CCSDS_ASCII_B = '^(\d{4})-(\d{3})T(\d{2}):(\d{2}):(\d{2})(\.(\d+)?Z?$'
    end

    methods
        function self = TimeValue(year, month, dom, hour, minute, sec, nanosec)
            % TIMEVALUE constructs a new instance of the object.
            %
            %   obj = TIMEVALUE(y, mo, dom, h, mi, s, ns)
            %
            % Each field is range-checked according to the definition: years in 1-9999, seconds 
            % in 0-60 (to allow for leap second representation) and other fields in their natural
            % ranges. Note that the date component is interpreted as proleptic Gregorian for the
            % purpose of determining the number of days in February.
            % Error 'OSFI:TimeValue:BadFieldValue' is raised if a value has an out-of-range value.
            %
            % See also PARSE
            self.year    = TimeValue.checkFieldRange('year', year, 1, 9999);
            self.month   = TimeValue.checkFieldRange('month', month, 1, 12);
            self.dom     = TimeValue.checkFieldRange('day-of-month', dom, ...
                    1, TimeValue.daysInProlepticGregorianMonth(year, month));
            self.hour    = TimeValue.checkFieldRange('hour', hour, 0, 23);
            self.minute  = TimeValue.checkFieldRange('minute', minute, 0, 59);
            self.sec     = TimeValue.checkFieldRange('sec', sec, 0, 60); % Allow '60' in leap seconds, which we don't check
            self.nanosec = TimeValue.checkFieldRange('nanosecond', nanosec, 0, 999999999);
        end

        %% Conversion functions: direct to string, then to datevec so that datestr can be used
        function str = char(self)
            % CHAR Formats this object for output
            %   str = CHAR(scalar) Returns a char array
            %   C   = CHAR(nonscalar) Returns a cellstr of the same shape
            %
            % See also PARSE
            
            if isscalar(self)
                str = self.formatSingleValue(true);
            else
                str = arrayfun(@char, self, 'UniformOutput', false);
            end
        end

        function str = string(self)
            % STRING Formats this object for output
            %   str = STRING(scalar) Returns a scalar string
            %   sa  = STRING(nonscalar) Returns a string array of the same shape
            %
            % See also PARSE
            
            if isscalar(self)
                str = self.formatSingleValue(false);
            else
                str = arrayfun(@string, self);
            end
        end

        function vm = datevec(self)
            % DATEVEC converts this object into datevec format
            %   v = DATEVEC(scalar) returns a 6-element vector [Y,M,D,H,MN,S]
            %   m = DATEVEC(nonscalar) returns a n-by-6 matrix where n=numel(nonscalar)
            %
            % All columns in the vector except for the last (seconds) are
            % integral. The following caveats apply:
            % - The nanoseconds from the object are added to the seconds
            %   element (6th column) seconds as a fractional part, which 
            %   may cause precision to be lost.
            % - No conversions are applied, so any out-of-range values will
            %   be passed along, INCLUDING the possible leap-second (60).
            [y,m,d,h,mn,s,ns] = deal([self(:).year]', [self(:).month]', [self(:).dom]', ...
                    [self(:).hour]', [self(:).minute]', [self(:).sec]', [self(:).nanosec]');
            vm = [y,m,d,h,mn,s+ns*1e-9];
        end

        %% Equality and relational operators: EQ and LT are implemented directly, then the rest in terms of them
        function res = eq(A,B)
            % EQ compares two TimeValue objects or arrays for equality. No conversions or
            % adaptations are performed, so month=4 and dom=0 is different from month=3 and dom=31
            %
            % See also LT
            res = arrayfun(@(a,b) (a.year==b.year) && (a.month==b.month) && (a.dom==b.dom) ...
                        && (a.hour==b.hour) && (a.minute==b.minute) && (a.sec==b.sec) && (a.nanosec==b.nanosec), ...
                    A, B);
        end
        function res = lt(A,B)
            % LT compares two TimeValue objects or arrays for and returns true if (or where) A<B.
            % No conversions or adaptations are performed, so month=3 and dom=31 is less-than month=4 and dom=0
            %
            % See also EQ
            res = arrayfun(@(a,b) a.year<b.year || ...
                        (a.year==b.year && (a.month<b.month || ...
                            (a.month==b.month && (a.dom<b.dom || ...
                                (a.dom==b.dom && (a.hour<b.hour ||...
                                    (a.hour==b.hour && (a.minute<b.minute || ...
                                        (a.minute==b.minute && (a.sec<b.sec || ...
                                            (a.sec==b.sec && a.nanosec<b.nanosec))))))))))), ...
                    A, B);
        end
        function res = ne(A,B)
            res = ~eq(A,B);
        end

        function res = le(A,B)
            res = ~(B<A);
        end
        function res = ge(A,B)
            res = ~(A<B);
        end
        function res = gt(A,B)
            res = B<A;
        end
    end

    methods (Access = private)
        function str = formatSingleValue(self, isChar)
            % Compute the value and width of the "fraction-of-second" field
            secFracVal = self.nanosec;
            secFracWidth = 9; % Prevent trailing zeros while keeping leading ones
            while secFracWidth > 3 && mod(secFracVal, 10) == 0
                secFracWidth = secFracWidth - 1;
                secFracVal = secFracVal / 10;
            end
            % Write, always including the fraction-of-second even if it is ".0"
            date = "%04d-%02d-%02dT%02d:%02d:%02d.%0*d";
            if isChar
                date = char(date);
            end
            str = sprintf(date, ...
                    self.year, self.month, self.dom, ...
                    self.hour, self.minute, self.sec, secFracWidth, secFracVal);
        end
    end

    methods (Static)
        function obj = parse(s)
            % PARSE Convert a *single* TIME-formatted string into a TimeValue instance.
            % The input must be in a recognized time format (CCSDS ASCII time code), or an
            % error is raised.
            %   obj = PARSE(s) processes a single string
            %
            % See also CHAR
            m = regexp(s, TimeValue.TIME_CCSDS_ASCII_A, 'once', 'tokens');
            matchMaD = ~isempty(m); % Day-of-month groups
            if ~matchMaD % Try with the day-of-year form
                m = regexp(s, TimeValue.TIME_CCSDS_ASCII_B, 'once', 'tokens');
            end
            matchDoY = ~matchMaD && ~isempty(m); % Day-of-year groups
            if ~matchMaD && ~matchDoY
                error('OSFI:TimeValue:BadFormat','Bad format for TIME value: %s', s)
            end
    
            % Common fields
            iM_year = 1;
            iM_hour = 4 * matchMaD + 3 * matchDoY;
            % The fit into each data type is forced by the number of digits (4, 2, etc.)
            year = str2double(m{iM_year});
            hour = str2double(m{iM_hour});
            min  = str2double(m{iM_hour+1});
            sec  = str2double(m{iM_hour+2});
    
            secfrac = m{iM_hour+3};
            if ~isempty(secfrac)
                secfrac = secfrac(2:end); % Don't consider the decimal period
                sfrac_digits = length(secfrac);
                if sfrac_digits > 9 % We only support up to nanos, so '000000000'..'999999999'
                    error('OSFI:TimeValue:BadSecFrac', 'Second fractions supported only up to ns, read: %s', secfrac)
                end
                % Multiply the fraction by 10^(9-number of digits) to convert to nanos
                nanos = str2double(secfrac) * 10^(9 - sfrac_digits);
            else % Fraction-of-second part unmatched, consider as '.0'
                nanos = 0;
            end
    
            % Parse either -month-day or -dayofyear
            if matchMaD
                month = str2double(m{iM_year+1});
                dom   = str2double(m{iM_year+2});
            else
                % FIXME: not mentioned in the examples of the E2E-ICD v2.4,
                % but it refers to "a CCSDS ASCII time format", so we should
                % support the day-of-year format, which is time code "type B".
                error('OSFI:NotImplemented', 'Day-of-year format for TIME parameters not implemented');
            end
    
            obj = TimeValue(year, month, dom, hour, min, sec, nanos); % Checks values
        end
    end

    methods (Static, Access=private)
        function val = checkFieldRange(name, val, min, max)
            % Private helper for checking the individual values in the TimeValue ctor
            if val < min || val > max
                error('OSFI:TimeValue:BadFieldValue', 'Invalid %s: %g', name, val)
            end
        end
        
        function days = daysInProlepticGregorianMonth(year, month)
            % Helper for the days in a month on an "eternally" extended Gregorian calendar (in years 0000-9999)
            if any(month == [1, 3, 5, 7, 8, 10, 12])
                days = 31;
            elseif any(month == [4, 6, 9, 11])
                days = 30;
            elseif month == 2
                isGrLeap = (mod(year, 4) == 0) && ((mod(year, 100) ~= 0) || (mod(year, 400) == 0));
                days = 28 + isGrLeap;
            else
                assert(false, 'Invalid month in internal function call: %02g/%04g', month, year)
            end
        end  
    end
end
