classdef Parameter < handle
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
    properties (Access = private)
        path;
        name;
        description;
        type;
        dims;
        value;
        min;
        max;
        %AN-010: Matrix  representation 
        complexType;
    end
    
    methods
        function obj = Parameter (path, name, description, type, dims, value, min, max, complexType)
            obj.path = path;
            obj.name = name;
            obj.description = description;
            obj.type = type;
            if ~isreal(dims)
                error('Dims argument should be a vector of dimensions');
            end
			obj.dims = dims;
			
            obj.value = value;
            obj.min = str2double(min);
            obj.max = str2double(max);
            obj.complexType = complexType;
        end
        function [path] = getPath (obj)
            path = obj.path;
        end       
        function [name] = getName (obj)
            warning('OSFI:Parameter:DeprecatedAPI', 'Ambiguous function, instead use getLocalName or getPath.');
            name = obj.name;
        end
        function [name] = getLocalName (obj)
            name = obj.name;
        end
        function [description] = getDescription (obj)
            description = obj.description;
        end
        function [type] = getType (obj)
            type = obj.type;
        end
        function type = getElementType(obj)
            % GETELEMENTTYPE returns the type of individual components of this parameter, as an
            % enum object.
            % See also PARAMTYPE.
            type = ParamType.(obj.type);
        end
        function [ndims] = getNdims (obj)
			ndims = length(obj.dims);
        end
        function [dims] = getDims (obj)
			dims = obj.dims;
        end
        function val = getRawValue(obj)
            val = obj.value;
        end
        function [value] = getValue (obj, pt)
            % GETVALUE returns the parsed value of the parameter.
            % val = p.GETVALUE() parses the value as its declared type
            % val = p.GETVALUE(pt) parses the value as the given type
            % See also GETARRAYVALUE.
            if nargin < 2
                pt = obj.getElementType();
            end
			if strcmp(obj.complexType,'ARRAY')
                value = obj.parseArrayNode(pt, obj.value, true, []);
			else
                value = obj.parseNonArrayValues(pt, obj.value);
			end
        end
        function [min] = getMin (obj)
            min = obj.min;
        end
        function [max] = getMax (obj)
            max = obj.max;
        end         
        function [complexType] = getComplexType (obj)
            complexType = obj.complexType;
        end  
        function print (obj)
            disp('');
            disp('Path:'); disp(obj.path); disp('');
            disp('Name:'); disp(obj.name); disp('');
            disp('Description:'); disp(obj.description); disp('');
            
            if(isempty(obj.complexType))
            	disp('Type:'); disp(obj.type); disp('');
            else
            	disp('Type:'); disp(obj.complexType); disp(''); 
            	disp('Element Type:'); disp(obj.type); disp(''); 
            end
            disp('Dims:'); disp(obj.dims); disp('');
            disp('Value:'); disp(obj.value); disp('');
            disp('Min:'); disp(obj.min); disp('');
            disp('Max:'); disp(obj.max); disp('');
        end

        function result = getArrayValue(obj, varargin)
            % GETARRAYVALUE parses the tree structure defined by this ARRAY parameter.
            % The return value may be an array of values, or a cell array containing structured
            % data. If varargin contains indices, they are applied before parsing, so the call
            %     v = p.GETARRAYVALUE(1,2);
            % is equivalent to:
            %     a = p.GETARRAYVALUE();
            % followed by either
            %     v = a{1}{2}; or v = a{1}(2);
            % depending on whether a{1} is a leaf node or not; except that the latter parses
            % the full tree while the former only parses the node selected by the given indices.
            if ~strcmp(obj.complexType,'ARRAY')
        		error('OSFI:InvalidState', 'Parameter: "%s" is not an ARRAY', obj.name);
            end
            pt = ParamType.(obj.type);
            % If any indices are given (in varargin), find that subnode before parsing
            indices = cell2mat(varargin);
            node = obj.value;
            idxIntoData = [];
            for iDim = 1:length(indices)
                if ParamParserComplex.isRawArrayValueNode(node)
                    % Trying to index into a data node: allowable if it is
                    % the last index (size is checked later)
                    if iDim == length(indices)
                        idxIntoData = indices(iDim);
                        break
                    else
                        error('OSFI:InvalidRank', 'Parameter "%s{%s}" is a data item, but there are extra indices (%s)', ...
                                obj.name, num2str(indices(1:iDim)), num2str(indices(iDim+1:end)));
                    end
                elseif indices(iDim) > length(node)
                    error('OSFI:OutOfBounds', 'Parameter "%s{%s}" has length %d, requested index %d', ...
                            obj.name, num2str(indices(1:iDim-1)), length(node), indices(iDim));
                end
                node = node{indices(iDim)};
            end
            % Parse tree without flattening it afterwards
            result = obj.parseArrayNode(pt, node, false, []);
            if ~isempty(idxIntoData)
                if indices(iDim) > length(node)
                    error('OSFI:OutOfBounds', 'Parameter "%s{%s}" has length %d, requested index %d', ...
                            obj.name, num2str(indices(1:end-1)), length(node), indices(iDim));
                end
                result = result(idxIntoData);
            end
        end
    end
    
    methods (Access = private)
        function parsed = parseNonArrayValues (obj, pt, raw)
            if prod(obj.dims) > 1
                vals = pt.splitValues(raw);
                dims = obj.dims;
                if prod(dims) ~= length(vals)
                    error('OSFI:Parameter:getValues', '%s: Expected %d values (%s) but found %d.', obj.name, ...
                            prod(dims), strjoin(arrayfun(@num2str, dims, 'UniformOutput', false), ','), length(vals));
                end
                if ~isscalar(dims)
                    if length(dims) == 2 % E2E-ICD gives columns first, flip for the reshape
                        dims = fliplr(dims);
                    end
                    % Matlab works in Fortran element order (1st dimension contiguous) but the
                    % serialized strings are in C element order (last dimension contiguous), permute
                    vals = reshape(vals, fliplr(dims));
                    vals = permute(vals, ndims(vals):-1:1);
                end
            else
                vals = strtrim(raw);
            end
            parsed = pt.parseValues(vals);
        end % parseNonArrayValues
        
        function [nodeVals] = parseArrayNode(obj, pt, node, flatten, curPath)
            % PARSEARRAYNODE creates a new cell array tree from node, with the same structure
            % but with parsed data nodes replacing the unparsed holding {dim, textdata}.
            % If flatten is true, the structure is flattened into a vector in depth-first order.
            % An error is raised if the number of values generated by splitting textdata according
            % to pt does not match the expected number (in dim), or if parsing fails.       
            makePath = @() strjoin(arrayfun(@(idx) ['{' num2str(idx) '}'], curPath(:)', 'UniformOutput', false), '');
            if ParamParserComplex.isRawArrayValueNode(node) % Data node
                vals = pt.splitValues(node{2}, node{1} == 0);
                if node{1} ~= length(vals)
                    error('OSFI:Parameter:getValues', '%s%s: Expected %d values but found %d.', ...
                            obj.name, makePath(), node{1}, length(vals));
                end
                try
                    nodeVals = pt.parseValues(vals);
                catch e
                    error('OSFI:Parameter:getValues', '%s%s: %s.', obj.name, makePath(), e.message);
                end
            else % Internal node, flatten all data nodes into single array
                nodeVals = arrayfun(@(iS, sub) obj.parseArrayNode(pt, sub{1}, flatten, [curPath iS]), ...
                        1:length(node), node(:)', 'UniformOutput', false);
                if flatten % Flatten nested cellstrs/arrays
                    nodeVals = [ nodeVals{:} ];
                end
            end
        end  %parseArrayNode
    end
end

