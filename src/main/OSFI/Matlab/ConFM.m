classdef ConFM < handle
% CONFM parses the XML configuration files provided by the end-to-end simulation framework,
% with the following elements:
% - Element <group_name>. This element can create nested groups to enclose parameters and
% provide scope information. Optional element.
% - Element parameter. This element can define the following attributes:
%     - Name. This is the parameter identifier. Names cannot contain spaces
%     - Description. Short definition or meaning for the parameter
%     - Type. Possible values for this attribute are: INTEGER, DOUBLE, BOOLEAN, STRING,
%         FILE or FOLDER. These are system-independent type of values, only intended for
%         GUI formatting validations.
%     - Value. This is the numerical, string or file location value of the parameter
%     - Units. Physical units of measurements if applicable. This attribute is optional
%     - Dims. Size of the dimensions. For example 1, 3 is a vector of three elements and
%         3, 3 is a square matrix of 3x3 elements. The first number refers to columns
%         and the second one to rows when describing a matrix. This attribute is optional
%         for scalar variables.
%         For variables of vector/matrix types, the value attribute must contain a
%         comma-separated (or blank-separated) list of values by rows, following this example:
%             A 3 columns by 2 rows matrix (dims = 3,2)
%             1    2    3
%             4    5    6
%         Will be represented as: 1,2,3,4,5,6.
%         String vectors must enclose each element in single quotes. For example
%         'a string' 'second string' 'third and last string'.

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
        DOMnode;
        paramMapList;
        returnStrings;
        convertToString;
    end
 
    methods
        function obj = ConFM (fileName, textType, textValue)
            % Accept only one option of configuration where the last two arguments are optional
            % T = ConFM(fileName, 'TextType', textValue) specifies the text type, where TextValue is either
            % 'char' or 'string'.

            if nargin == 0
                error('ConFM must be called with the path to the XML file');
            elseif nargin == 1
                value = 'char';
            elseif nargin == 3
                if strcmpi(textType, 'TextType')
                    value = textValue;
                else
                    error('OSFI:BadArgument', 'Optional argument not recognised')
                end
            else
                error('OSFI:ConFM:InvalidNargin', 'ConFM: either 1 or 3 arguments expected, got %d', nargin');
            end

            % Decide whether to use char or strings based on TextValue
            if strcmpi(value, 'string')
                obj.returnStrings = true;
                obj.convertToString = @string;
            elseif strcmpi (value, 'char')
                obj.returnStrings = false;
                obj.convertToString = @char;
            else
                error('OSFI:ConFM:InvalidArgTypes', 'Input value must be char vectors or strings arrays')
            end

            obj.parseFile(fileName)
        end
        function parseFile (obj, fileName)
            % PARSEFILE reads the given file and replaces the contents of this ParamReader.

            %% Read the XML into a DOM document
            % Note that the Matlab function xmlread does the same but appears NOT to set the
            % "namespaceAware" feature, so XSD validation may fail. Thus, manually call the
            % Java-based XML parsing.
            dbf = javaMethod('newInstance', 'javax.xml.parsers.DocumentBuilderFactory');
            dbf.setCoalescing(true);
            dbf.setNamespaceAware(true);
            dbf.setIgnoringComments(true);
            db = dbf.newDocumentBuilder();
            obj.DOMnode = db.parse(ConFM.javaFile(fileName));
            %% Parse and store the map of path => Parameter objects
            obj.paramMapList = obj.getParameters();
        end
        function [parameter] = getParameter (obj, strPath)
            parameter = obj.paramMapList(strPath);
        end
        function [paramMap] = getParameters (obj)
            % Create map
            paramMap = containers.Map();
            
            % Find all <parameter ... />
            xRoot   = obj.DOMnode.getDocumentElement;
            paramElements = xRoot.getElementsByTagName ('parameter');
            nParameters = paramElements.getLength;
            
            rootNodeName = char(xRoot.getNodeName);
            
            % Iterate over all the parameters (DOM objects are Java, indices are 0-based)
            for iP=1:nParameters
                pElement = paramElements.item(iP-1);
                parent = pElement;
                
                % getElementsByTagName recursed into all parameter sub-elements, so skip
                % those elements whose parent is also a "parameter" element.
                if strcmp(pElement.getParentNode.getNodeName,'parameter')
                    continue
                end
                
                % Get path
                ePath = obj.convertToString('');
                nodeName = obj.convertToString('');
                while ~strcmp(nodeName, rootNodeName)
                    if strlength(nodeName) ~= 0
                        ePath = strcat(nodeName, '.', ePath);                            
                    end
                    parent = parent.getParentNode;
                    nodeName = obj.convertToString(parent.getNodeName);
                end

                % Get the required attributes and place them as fields in a struct
                paramAttrs = {'name', 'type', 'description', 'dims', 'min', 'max', 'units'};
                paramAttrs(2,:) = cellfun(@(attr) { obj.convertToString(pElement.getAttribute(attr)) }, paramAttrs(1,:));
                paramAttrs = struct(paramAttrs{:});

                ePath = strcat(ePath, paramAttrs.name);
                try
                    %If parameter has MATRIX or ARRAY type it needs a different processing
                    if any(strcmp(paramAttrs.type, {'ARRAY', 'MATRIX'}))
                        paramInfo = ParamParserComplex(ePath, pElement, obj.convertToString);
                        [paramAttrs.type, paramAttrs.dims, paramAttrs.value, paramAttrs.complexType] = ...
                            deal(paramInfo.getType(), paramInfo.getListDims(), paramInfo.getValue(), paramInfo.getComplexType());
                    else
                        if pElement.hasAttribute('value') && pElement.hasChildNodes
                            error('A parameter cannot have both the value attribute and a non-empty content');
                        elseif pElement.hasChildNodes % getTextContent clears comments and joins text nodes
                            paramAttrs.value = obj.convertToString(pElement.getTextContent());
                        else
                            paramAttrs.value = obj.convertToString(pElement.getAttribute('value'));
                        end
                        paramAttrs.dims = ConFM.parseDims(paramAttrs.dims);
                        paramAttrs.complexType = obj.convertToString('');
                    end
                    % TODO: implement the "units" field
                    parameter = Parameter (ePath, paramAttrs.name, paramAttrs.description,...
                            paramAttrs.type, paramAttrs.dims, paramAttrs.value, ...
                            paramAttrs.min, paramAttrs.max, paramAttrs.complexType);
                    % Add to map
                    paramMap (ePath) = parameter;
                catch e
                    newErr = MException('OSFI:ConFM:ParamReadFailure', ...
                            'Failed to parse parameter "%s". %s', ePath, e.message);
                    newErr.addCause(e);
                    throw(newErr);
                end
            end
        end

        function valOk = validateAgainst(obj, xsdFile)
            % VALIDATEAGAINST validates the XML file read against the given XSD schema, 
            % logging the errors found (if any) using the Logger facilities.
            %
            % valPasses = obj.VALIDATEAGAINST(path_to_xsd)
            %
	        % The result value is true if the validation succeeds, false otherwise.
            % If the schema cannot be loaded, an error is thrown.
            %
            % See also VALIDATEAGAINSTINTERNALSCHEMA.
            valOk = ConFM.doValidate(obj.DOMnode, xsdFile);
        end
        function valOk = validateAgainstInternalSchema(obj)
            % VALIDATEAGAINSTINTERNALSCHEMA validates the XML file read against the XSD schema
            % that is referenced in the document, logging the errors found (if any) using the Logger facilities.
            %
            % valPasses = obj.VALIDATEAGAINSTINTERNALSCHEMA()
            %
	        % The result value is true if the validation succeeds, false otherwise.
            %
            % See also VALIDATEAGAINST.
            valOk = ConFM.doValidate(obj.DOMnode);
        end
    end
    
    methods (Access = private)
        function [child] = getChildNode (~, father, string)
           children = father.getChildNodes;
           for i=0:children.getLength-1
               childNode = children.item(i);
               if childNode.getNodeType == childNode.ELEMENT_NODE
                   nodeName = childNode.getNodeName;
                   if (nodeName == string)
                       child = childNode;
                       return;
                   end
               end
           end

           % If not found... return the father
           child = father;
        end   
        function [parameter] = getParameterByName (~, father, name)
           children = father.getChildNodes;
           for i=0:children.getLength-1
               childNode = children.item(i);
               if childNode.getNodeType == childNode.ELEMENT_NODE
                   parameterName = childNode.getAttribute('name');
                   if (parameterName == name)
                       parameter = childNode;
                       return;
                   end
               end
           end

           % If not found... return the father
           parameter = [];
        end
    end

    methods (Static = true)
        function dims = parseDims(aDims, expectedRank)
            % PARSEDIMS Parse a "dims" attribute into an array of dimensions. If the expectedRank argument is,
            % provided, the function also asserts that the number of dimensions is equal to that argument.
            if isempty(strtrim(aDims))
                dims = [];
            else
                dims = cellfun(@str2double, strsplit(aDims));
            end
            if nargin > 1 && length(dims) ~= expectedRank
                error('Expected number of dimensions is %d, but dims is ''%s\''', expectedRank, aDims)
            end
        end
    end
    methods (Static = true, Access = private)
        function file = javaFile(f)
            % JAVAFILE builds a java.io.File for path f, taking into account the fact that Java does NOT
            % share the CWD with matlab, so if the path is not absolute, Matlab's CWD is used.
            file = javaObject('java.io.File', f);
            if ~file.isAbsolute()
                file = javaObject('java.io.File', fullfile(pwd, f));
            end
        end

        function schema = loadSchema(schemaFile)
            % LOADSCHEMA builds a javax.xml.validation.Schema object from the given file path, or, if not provided,
            % a special XSD schema object that will use the schema referenced by the XML document.
            W3C_XML_SCHEMA_NS_URI = 'http://www.w3.org/2001/XMLSchema';
            schemaFactory = javaMethod('newInstance', 'javax.xml.validation.SchemaFactory', W3C_XML_SCHEMA_NS_URI);
            if nargin == 0
                schema = schemaFactory.newSchema(); % Load the internal schema in the XML document.
            else
                xsdFile = ConFM.javaFile(schemaFile);
                schema = schemaFactory.newSchema(xsdFile);
            end
        end

        function valOk = doValidate(docObj, varargin)
            % DOVALIDATE implements the validation functions. Returns true if validation succeeds and false otherwise.
            % Throws an error if the schema cannot be loaded, or if validation throws anything other than SAXParseException.
            try
                schema = ConFM.loadSchema(varargin{:});
            catch e
                error('OSFI:ConFM:validate:loadSchema', 'Could not load the schema: %s', e.message)
            end
            val = schema.newValidator();
            val.setErrorHandler([]) % With this, we lose warnings (we cannot implement ErrorHandler from Matlab)

            % FIXME since we cannot implement LSResourceResolver either, we at least _try_ to make the
            % "Java CWD" match Matlab's for the duration. HOWEVER, Java actually does NOT allow users
            % to reliably modify the CWD, so this is a big hack which works as long as the validator
            % accesses files in a way that reads 'user.dir' dynamically instead of caching it.
            oldCWD = javaMethod('getProperty', 'java.lang.System', 'user.dir');
            try
                javaMethod('setProperty', 'java.lang.System', 'user.dir', pwd);
                val.validate(javaObject('javax.xml.transform.dom.DOMSource', docObj));
                javaMethod('setProperty', 'java.lang.System', 'user.dir', oldCWD); % Matlab has no finally!
                valOk = true;
            catch e % Should be a SAXParseException
                javaMethod('setProperty', 'java.lang.System', 'user.dir', oldCWD); % Matlab has no finally!
                if ~isa(e,'matlab.exception.JavaException') || ~isa(e.ExceptionObject, 'org.xml.sax.SAXParseException')
                    rethrow(e);
                end
                javaExc = e.ExceptionObject;
                log = Logger(true); % OSFI-Matlab EHLog has instance, not static functions
                log.error(['OSFI::ConFM::validate. Error ', char(javaExc.getMessage())])
                valOk = false;
            end
        end
    end
     
end
