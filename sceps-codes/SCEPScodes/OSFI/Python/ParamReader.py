#
# openSF Integration Libraries (OSFI)
# Deimos Space, S.L.U.
# 
# This file is part of OSFI. OSFI is free software; you can redistribute it
# and/or modify it under the terms of the 'ESA Software Community Licence Permissive' as
# published by the European Space Agency; either version 2.4 of the License,
# or (at your option) any later version. You should have received a
# copy of the 'ESA Software Community Licence Permissive - v2.4' along with this program
# or one can be found at <http://eop-cfi.esa.int/index.php/docs-and-mission-data/licensing-documents>.
#

try:  # The Python-bundled etree does not implement XMLSchema, so try to use the lxml version
    from lxml import etree as ET

    XML_PARSER_LIB = 'lxml'
except ImportError:  # Validation will be unavailable
    import xml.etree.ElementTree as ET

    XML_PARSER_LIB = 'xml.etree'

try:
    from OSFI import Logger
    from OSFI.Parameter import Parameter
except ImportError:
    import Logger
    from Parameter import Parameter


class ParamParsingError(ValueError):
    def __init__(self, path, message):
        super(ParamParsingError, self).__init__(path + ": " + message)


class ParamReader(object):
    """
    This ConFM class parses the XML configuration files provided by the openSF HMI formatted
	with the following elements:
	- Element <group_name>. This element can create nested groups to enclose parameters and
	provide scope information. Optional element.
	- Element parameter. This element can define the following attributes:
	        - Name. This is the parameter identifier. Names cannot contain spaces
	        - Description. Short definition or meaning for the parameter
	        - Type. Possible values for this attribute are: INTEGER, DOUBLE, BOOLEAN, STRING,
	            FILE or FOLDER. These are system-independent type of values, only intended for
	            GUI formatting validations.
	        - Value. This is the numerical, string or file location value of the parameter
	        - Units. Physical units of measurements if applicable. This attribute is optional
	        - Dims. Size of the dimensions. For example 1, 3 is a vector of three elements and
	            3, 3 is a square matrix of 3x3 elements. The first number refers to columns
	            and the second one to rows when describing a matrix. This attribute is optional
	            for scalar variables.
	            For variables of vector/matrix types, the value attribute must contain a
	            comma-separated (or blank-separated) list of values by rows, following this example:
	                A 3 columns by 2 rows matrix (dims = 3,2)
	                1    2    3
	                4    5    6
	            Will be represented as: 1,2,3,4,5,6.
	            String vectors must enclose each element in single quotes. For example
	            'a string' 'second string' 'third and last string'.
    """

    def __init__(self, xmlFile, xsdFile=None):
        """
        Default class constructor.
        Creates an instance of the ParamReader class, parses a given configuration XML file and
        stores every valid parameter.</p>
        <p>Validating against an XSD schema is optional.</p>
        @param xmlFile - Configuration XML file
        @param xsdFile - IGNORED, will generate a deprecation warning
        @throw exception if error while parsing
        """
        if xsdFile is not None:
            import warnings
            msg = DeprecationWarning(
                "OSFI::ParamReader. The xsdFile argument has always been ignored; if schema validation is " \
                + "actually desired then build ParamReader without it and use the validateAgainst methods.")
            Logger.warning(msg.args[0])
            warnings.warn(msg)
        self.__params = {}
        try:
            # Using a custom parser is required by lxml, but xml.etree removes comments by default and actually
            # does not even accept the remove_comments argument to the parser, so just use None (default)
            try:
                parser = ET.XMLParser(remove_comments=True)
            except Exception as e:
                parser = None  # If using the CPython implementation, just use the default parser
            self.__tree = ET.parse(xmlFile, parser)
            self.__root = self.__tree.getroot()
            self.__traverseNodes(None, self.__root)
        except Exception as e:
            Logger.error("OSFI::ParamReader::constructor. Error parsing file '{0}': {1}".format(xmlFile, str(e)))
            raise e

    def getParameter(self, paramName):
        """
         <p>Gets a parameter object corresponding with the first occurrence of a parameter name into a configuration file.</p>
         @param paramName Parameter name (full path and name)
         @return .Parameter instance.
        """
        return self.__params.get(paramName)

    def getParameters(self, groupName):
        """
         Gets the list of parameters under a certain group.
         @param groupName - full name of the group.
        """
        parameters = []
        for key, value in self.__params.items():
            if groupName in key:
                parameters.append(value)

        return parameters

    def getAllParameters(self):
        """
          Gets a map containing all stored parameters.
        """
        return self.__params

    def existParameter(self, paramName):
        """
         <p>Checks the existence of a parameter within a configuration file.</p>
        @param paramName Parameter name (full path and name)
        @return boolean True if the parameter exists False otherwise.
        """
        result = paramName in self.__params.keys()
        return result

    def setParameter(self, paramName, value):
        """
         Changes the stored value of a certain parameter.
         Shows an error if cannot find parameter.
         @param paramName - full name of the parameter
         @param value - new value
         @throw exception if value is invalid
        """
        self.__params[paramName] = value

    def write(self):
        """
         <p>Prints a textual representation of the list of Parameters on screen.</p>
        """
        for key, value in self.__params.items():
            Logger.info("Complete name : " + key)
            value.write()

    def __traverseNodes(self, parentPath, element):
        """
         <p>Traverses the DOM tree nodes and stores the parameters of a configuration file.</p>
         @param parentPath tuple with the names of parent elements, or None if this is the root element
         @param element DOM document element, must not be a "parameter" element.
         @throw exception if error while parsing.
        """
        # If this is the root element, its tag should _not_ be included in the path
        myPath = parentPath + (element.tag,) if parentPath is not None else ()
        myName = ".".join(myPath) if myPath else "(root)"
        try:
            # Find all "parameter" elements that are direct children of this element
            for param in element.iterfind('parameter'):
                pName = param.get('name')
                if pName is None:
                    raise ValueError("Parameter without a name attribute under " + myName)
                completeName = ".".join(myPath + (pName,))
                pType = param.get('type')
                if pType in ('MATRIX', 'ARRAY'):
                    try:
                        from OSFI.ParamParserComplex import ParamParserComplex
                    except ImportError:
                        from ParamParserComplex import ParamParserComplex
                    p = ParamParserComplex(param, completeName)
                    pComplexType = pType  # Store for ARRAY types
                    pType, pValue = p.type, p.value
                else:
                    if param.text is not None and "value" in param.attrib:
                        raise ValueError("A parameter cannot have both the value attribute and a non-empty content")
                    # Try the node content (new convention), then the attribute value
                    pValue = param.text if param.text is not None else param.get('value')
                    pComplexType = 'NONE'
                self.__readParameter(param, pType, pComplexType, pValue, completeName)

            # Find non-parameter elements
            for subElem in element.iterfind('*'):
                if subElem.tag == 'parameter':
                    continue  # Skip
                self.__traverseNodes(myPath, subElem)
        except ParamParsingError:
            raise  # Path already noted
        except Exception as e:
            Logger.error("OSFI::ParamReader::traverseNodes. Error at/under {0}: {1}".format(myName, str(e)))
            raise ParamParsingError(myName, str(e))

    @staticmethod
    def _parseDims(aDims, expectedRank=None):
        """
        Parse a "dims" attribute into a tuple of dimensions. If the expectedRank argument is,
        provided, the function also asserts that the number of dimensions is equal to that argument.
        """
        dims = tuple(int(d) for d in aDims.split())
        if expectedRank is not None and len(dims) != expectedRank:
            raise ValueError("Expected number of dimensions is {0}, but dims is '{1}'".format(expectedRank, aDims))
        return dims if expectedRank != 1 else dims[0]  # For a single _expected_ dim, return the number directly

    def __readParameter(self, currentElement, aType, aComplexType, pValue, completeName):
        """
         Reads the XML element containing the parameter information.
         @param currentElement - XML element to read, already known to be a parameter element with a name
         @throw exception if element is not valid.
        """
        try:
            if pValue is None:
                pValue = ""
            aName = currentElement.get("name")  # Name cannot be empty or traverseNodes would fail
            aDescription = currentElement.get("description", "")
            aUnits = currentElement.get("units", "")
            aMin = currentElement.get("min", "")
            aMax = currentElement.get("max", "")
            aDims = ParamReader._parseDims(currentElement.get("dims", ""))  # list[int], may be empty
            param = Parameter(name=aName, description=aDescription, elType=aType, value=pValue, aDims=aDims,
                              aUnits=aUnits, aMin=aMin, aMax=aMax, aParamComplexType=aComplexType)
            self.__params[completeName] = param

        except ParamParsingError as e:
            raise  # Path of error has already been noted
        except Exception as e:
            Logger.error("OSFI::ParamReader::readParameter. Error at {0}: {1}".format(completeName, str(e)))
            raise ParamParsingError(completeName, str(e))

    def __doValidate(self, xsdPathOrUrl):
        """
        Helper function that implements the validateAgainst methods. It parses the schema source given
        (which may throw ET.ParseError), creates the schema object (which may raise NotImplementedError)
        and perform validation, logging all problems found. Returns True if validation succeeds without errors.
        """
        schemaTree = ET.parse(xsdPathOrUrl)
        try:  # Check that we are using lxml, or fail since the base Python ElementTree does not implement XMLSchema
            XMLSchema = ET.XMLSchema
        except AttributeError:
            raise NotImplementedError("XML ElementTree library does not implement XSD validation - try installing lxml")
        schema = XMLSchema(schemaTree)
        result = schema.validate(self.__tree)
        if not result:  # Log all problems found
            for err in schema.error_log:
                logfun = Logger.warning if err.level == ET.ErrorLevels.WARNING else Logger.error
                logfun("OSFI::ParamReader::validate. {0} at {1}:{2}:{3}. {4}".format(
                    err.level_name, err.filename, err.line, err.column, err.message))
        return result

    def validateAgainst(self, xsdFile):
        """
        Validate the XML document read against the given schema, logging problems found.
        @return True if validation succeds, False otherwise.
        @throw Exception If the named schema cannot be loaded
        @throw NotImplementedError If the ElementTree library does not implement XSD validation
        """
        return self.__doValidate(xsdFile)

    def validateAgainstInternalSchema(self):
        """
        Validate the XML document read against the schema linked in the document itself, logging problems found.
        @return True if validation succeds, False otherwise.
        @throw ValueError If the document does not link to a XSD schema using the schemaLocation/noNamespaceSchemaLocation
            attributes from the "http://www.w3.org/2001/XMLSchema-instance" namespace.
        @throw Exception If the linked schema cannot be loaded
        @throw NotImplementedError If the ElementTree library does not implement XSD validation
        """
        # Unlike in Java, there is no ET.XMLSchema() or ET.XMLSchema(None) to build a schema
        # that auto-locates the actual XSD from hints in the file, we have to find it ourselves
        XMLSchemaNamespace = '{http://www.w3.org/2001/XMLSchema-instance}'
        schemaLink = self.__root.get(XMLSchemaNamespace + 'schemaLocation')
        if schemaLink is None:
            schemaLink = self.__root.get(XMLSchemaNamespace + 'noNamespaceSchemaLocation')
        if schemaLink is None:
            raise ValueError("XML file does not link to a schema")
        return self.__doValidate(schemaLink)
