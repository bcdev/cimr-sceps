classdef ParamParserComplex < handle
% PARAMPARSERCOMPLEX is an OSFI-internal class that deals with parameters of type ARRAY and MATRIX
% XML constraints of new format: 
% - New types: MATRIX and ARRAY 
% - Must have 'dims' attribute 
% - all parameter child must have 'parameter' name 
% - all 'elementType' attribute must be the same 
% 
% ARRAY: 
% - 'dims' in each xml line must have only one value 
% - Maximum of 3 dimensions 
% - number of child/values must match the number defined in 'dims' 
% - intermediate child don't need 'elementType' attribute 
% 
% MATRIX: 
% - 'dims' must have two values 
% - all childs of parameter must have 'elementType' attribute 
% - number of columns and lines must match the number defined in 'dims'

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
		fqName;
		type;
        value = {''};
        listDims;
        maxDims;
		complexType;
		convertToString;
	end    
    
    methods
		function obj = ParamParserComplex (fqName, parameterNode, convertToString) %parameterNode is an Element
			% convertToString is used to decide whether to return the parameter name, complexType, type and value
			% as char vector or string for both, matrix and array
			obj.convertToString = convertToString;
			obj.fqName = obj.convertToString(fqName);
			obj.complexType = obj.convertToString(parameterNode.getAttribute('type'));
			if ~any(strcmp(obj.complexType, {'ARRAY', 'MATRIX'}))
				error('OSFI:ConFM:BadStructuredType', 'Unknown structured type "%s"', obj.complexType);
			end

			%% Determine elementType by checking _all_ descendent "parameter" elements
			allSubParams = parameterNode.getElementsByTagName('parameter');
			childTypes = cell(1, allSubParams.getLength()+1); % L+1 so we can read it from the root too
			for iSub = 1:allSubParams.getLength()
				childTypes{iSub} = char(allSubParams.item(iSub-1).getAttribute('elementType'));
			end
			childTypes{end} = char(parameterNode.getAttribute('elementType')); % Allow root element to provide elType too
			childTypes = childTypes(~cellfun(@isempty, childTypes));
			if isempty(childTypes) || ~all(strcmp(childTypes{1}, childTypes(2:end))) % Not empty, all present values match
				error('OSFI:ConFM:BadElementType', 'Could not determine the element type (missing, or conflicting values)');
			end
			obj.type = obj.convertToString(childTypes{1});

			if strcmp(obj.complexType, 'MATRIX')
				[obj.value, obj.listDims] = obj.readMatrix(parameterNode);
				return
			elseif strcmp(obj.complexType, 'ARRAY')
				obj.value = obj.readArray(parameterNode, []);
				obj.listDims = ParamParserComplex.parseArrayDims(obj.value);
			end
        end  %ParamParserComplex

        function [type] = getType (obj)
            type = obj.convertToString(obj.type);
        end
		function [value] = getValue (obj)
			% Value may be either a string (for matrices) containing the data in row-major order,
			% or a structure of nested cell arrays (for ARRAY parameters).
			value = obj.value;
        end
		function [dims] = getListDims (obj)
            dims = obj.listDims;
        end
        function [complexType] = getComplexType (obj)
            complexType = obj.complexType;
        end
    end  %methods
    
    methods (Access = private)
		function [value, dims] = readMatrix(self, paramEl)
			% Reads a MATRIX-typed parameter and returns a string corresponding to the
			% concatenation of the contents of the rows, separated by spaces. A MATRIX param
			% thus ends up with the same "value" as if it had been written in the old format.
			%fprintf("reading matrix %s", self.fqName)
			dims = ConFM.parseDims(char(paramEl.getAttribute('dims')), 2);
			[nC, nR] = deal(dims(1), dims(2)); % E2E-ICD: columns go first in a matrix
			rows = ParamParserComplex.getNamedChildNodes(paramEl, 'parameter'); % Not recursive, only direct children
			if length(rows) ~= nR
				error('OSFI:ConFM:BadSize', '%s: Number of rows is not the expected (%d vs %d)', self.fqName, length(rows), nR);
			end
			
			value = '';
			for iR = 1:length(rows)
				r = rows{iR};
				makeRowPath = @() [ self.fqName '(' num2str(iR) ')' ]; % Create if needed
				if ~strcmp(char(r.getAttribute('type')), 'ARRAY')
					error('OSFI:ConFM:BadType', '%s: Row does not have the expected type ARRAY', makeRowPath())
				end
				% Rows can omit the dims (assumed to be nC), but if they do specify it, it must match nC
				if r.hasAttribute('dims') && ConFM.parseDims(char(r.getAttribute('dims')), 1) ~= nC
					error('OSFI:ConFM:BadSize', '%s: Number of columns changed in row vs matrix', makeRowPath())
				end
				subElems = r.getElementsByTagName('parameter');
				if subElems.getLength() ~= 0 % Contains sub-elements!
					error('OSFI:ConFM:BadStructure', '%s: Row contains sub-parameter elements instead of data', makeRowPath())
				end
				% XXX using getTextContent allows us to ignore comment nodes, but it also means that if the
				% user has added any subelements we will get their text too.
				data = char(r.getTextContent());
				% TODO split data as the target type to actually verify len against nC? However, this
				% behaviour (fail-fast, preventing whole file parsing) does not match that of old vectors.
				value = [ value, ' ', data ];
			end
			value = self.convertToString(value);
		end

		function value = readArray(self, curElem, curPath)
			% Reads an ARRAY-typed parameter and sets value to a structure of cell arrays representing
			% the *unparsed* parameter values. Leaf nodes are represented by {dim, raw-text}.
			curFq = [ self.fqName, strjoin(arrayfun(@(idx) ['{' num2str(idx) '}'], curPath(:)', 'UniformOutput', false), '') ];
			%fprintf("reading array element %s", curFq)
			try
				dim = ConFM.parseDims(char(curElem.getAttribute('dims')), 1);
				subElems = ParamParserComplex.getNamedChildNodes(curElem, 'parameter'); % Not recursive, only direct children

				if isempty(subElems) % Leaf node, read values and split
					% value = ParamType.valueOf(self.type).splitValues(curElem.text)
					% TODO split data as the target type to actually verify len against dim? However, this
					% behaviour (fail-fast, preventing whole file parsing) does not match that of old vectors.
					value = { dim, self.convertToString(curElem.getTextContent())};
				else % Non-leaf, parse subelements
					if length(subElems) ~= dim
						error('OSFI:ConFM:BadSize', '%s: Number of sub-elements is not the expected (%d vs %d)', ...
								curFq, length(subElems), dim);
					end
					value = arrayfun(@(iS, s) self.readArray(s{1}, [curPath iS]), ...
							1:length(subElems), subElems(:)', 'UniformOutput', false);
				end
				% disp(value)
			catch e
				if strncmp(e.identifier, 'OSFI:ConFM:', length('OSFI:ConFM:'))
					rethrow(e) % Path already noted, pass along
				else
					error('OSFI:ConFM:ParamReader', '%s: Parse failure. At %s@%d: %s', curFq, e.stack(1).name, e.stack(1).line, e.message)
				end
			end
		end
	end

	methods(Static, Access = private)
		function [childNodes] = getNamedChildNodes(node, name)
			nl = node.getChildNodes();
			childNodes = cell(1,nl.getLength());
			for iN = 1:length(childNodes)
				cur = nl.item(iN-1); % API is 0-based (Java)
				if strcmp(char(cur.getNodeName()), name)
					childNodes{iN} = cur;
				end
			end
			childNodes = childNodes(~cellfun(@isempty, childNodes));
		end

		function [rectEnvDims] = parseArrayDims(curVal)

			if ParamParserComplex.isRawArrayValueNode(curVal)
				% Data node, dimension size is the first element
				rectEnvDims = curVal{1};
			else % Internal node, current dimension size is the length of subelements (recurse into them)
				allSubDims = cellfun(@(sub) ParamParserComplex.parseArrayDims(sub), curVal(:), 'UniformOutput', false);
				% Before stacking each list of sizes in the cell array we have to extend all items up 
				% to the deepest dimension, so we can call max across them all.
				deepestDim = max(cellfun(@length, allSubDims));
				allSubDims = cell2mat(cellfun(@(d) [d zeros(1, deepestDim-length(d))], allSubDims, 'UniformOutput', false));
				furtherDims = max(allSubDims, [], 1); % Keep the largest value for each dimension
				rectEnvDims = [length(curVal), furtherDims];
			end
		end
	end

	methods (Static = true)
		function isRawArray = isRawArrayValueNode(curVal) 
			isRawArray = length(curVal) == 2 && isreal(curVal{1}) && (ischar(curVal{2}) || isstring(curVal{2}));
		end
	end
end

