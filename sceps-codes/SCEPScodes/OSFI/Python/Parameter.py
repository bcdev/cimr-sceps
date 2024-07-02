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

from __future__ import division
import os

try:
    from OSFI import Logger
    from OSFI.ParamType import ParamType, parseSingleFile
except ImportError:
    import Logger
    from ParamType import ParamType, parseSingleFile


class Parameter(object):
    """
     This class encapsulates every element of information from the configuration file.
     <p>A parameter is a self-describing element that couples a name and a value, describing
     its contents and adding syntactic and semantic constraints.</p>
    """
    
    """
     Characters for delimiting elements in a parameter
    """
    __delimiter = '[ ,]+'
    
    """
     Characters for delimiting elements in a string parameter
    """
    __stringDelimiter = '[ ,\']*'
    
    def __init__(self, name, description, elType, value, aDims=(), aUnits='', aMin='', aMax='',
                 aParamComplexType='NONE'):
        """
        <p>Complete constructor.</p>
         @param name Name.
         @param description Description.
         @param type Type. See ParamType. Either a ParamType value or a string can be used
         @param value Value. May be a string or, for ARRAY parameters, the format returned by ParamParserComplex.__readArray.
         @param aDims Dimensions sizes. List or tuple of integers.
         @param aUnits Units. Blank or "-" for unitless parameters.
         @param aMin Minimum valid value.
         @param aMax Maximum valid value.
         @param aParamComplexType ParamComplexType. One of the following valid options::  <ul complexType="disc">  <li>NONE;</li> <li>ARRAY;</li>  <li>MATRIX;</li> </ul>
        """
        super(Parameter, self).__init__()
        self.__name = name
        self.__description = description
        # XXX consider that this causes fail-fast if a single parameter has an invalid type
        self.__type = elType if isinstance(elType, ParamType) else ParamType.valueOf(elType)
        self.__paramComplexType = ParamComplexType.valueOf(aParamComplexType)
        self.__value = value
        self.__units = aUnits
        self.__min = aMin
        self.__max = aMax
        self.__dims = aDims

        if self.__paramComplexType == 'ARRAY' and not isinstance(self.__value, (tuple, list)):
            raise ValueError("Invalid value for ARRAY-valued parameter, must be structure")
        

    def write(self):
        """
          <p>Prints out to standard output a textual description of the parameter.</p>
        """
        output = "Parameter ("
        output += self.__name + ", "
        output += self.__description + ", "
        output += ParamType.getParamTypeName(self.__type) + ", "
        if self.__paramComplexType != ParamComplexType.NONE:
            output += ParamComplexType.getParamComplexTypeName(self.__paramComplexType) + ", "
        output += str(self.__value) + ", "
        output += self.__units + ", "
        output += str(self.__min) + ", "
        output += str(self.__max) + ", "
        for i in range(0, len(self.__dims)):
            output +=  " " + str(self.__dims[i])
        output += ")" + '\n'
        Logger.getInfoStream().write(output)

    def getRawValue(self):
        """
        Returns the unparsed parameter value, as a string.
        """
        return str(self.__value)
    
    def getValue(self, asType=None):
        """
        Returns the parsed value of the parameter.
        """
        if isinstance(self.__value, ArrayNode):
            return self.getArrayValue(asType=asType)
        else:
            nd = len(self.__dims)
            pt = asType if asType is not None else self.__type
            if nd == 0:
                return self.__parseValAsScalar(pt)
            elif nd == 1:
                return self.__parseValAsVector(pt)
            elif nd == 2:
                return self.__parseValAsMatrix(pt)
            else:
                raise ValueError("Invalid number of dimensions for non-ARRAY")
    
    def __parseSingleValAs(self, paramType, value):
        """
        Helper that wraps paramType.parseSingleValue in a try block with adequate logging of errors.
        """
        try:
            return paramType.parseSingleValue(value)
        except Exception as e:
            Logger.error("OSFI::Parameter. Parameter: '{0}'. Value '{1}' is not valid for type {2}: {3}"
                        .format(self.__name, value, str(paramType), str(e)))
            return None
    
    def __parseValAsScalar(self, paramType):
        """
        Just calls __parseSingleValAs which already logs an error on failure.
        """
        return self.__parseSingleValAs(paramType, self.__value)

    def getIntValue(self):
        """
         <p>Gets the value of the parameter as an integer type.</p>
         <p>Valid when there is only one value.</p>
         @return Parameter value as integer, or None if the parse fails
        """
        return self.__parseValAsScalar(ParamType.INTEGER)

    def getDoubleValue(self):
        """
         <p>Gets the value of the parameter as a double type.</p>
         <p>Valid when there is only one value.</p>
         @return Parameter value as double, or None if the parse fails
        """
        return self.__parseValAsScalar(ParamType.FLOAT)

    def getBooleanValue(self):
        """
         <p>Gets the value of the parameter as a boolean type.</p>
         <p>Valid when there is only one value.</p>
         @return TRUE if value = "TRUE", or None if the parse fails
        """
        return self.__parseValAsScalar(ParamType.BOOLEAN)
    
    def getStringValue(self):
        """
        <p>Gets the value of the parameter as a string type.</p>
         @return Parameter value as string, or None if the parse fails
        """
        return self.__parseValAsScalar(ParamType.STRING)
    
    def getFileValue(self):
        """
        <p>Gets the value of the parameter as a file.</p>
        <p>Checks the file validity</p>
        @return Parameter value as a string, or None if the parse fails
        """
        return self.__parseValAsScalar(ParamType.FILE)

    def getTimeValue(self):
        """
         <p>Gets the value of the parameter as a timecode.</p>
         <p>Valid when there is only one value.</p>
         @return Parameter value as a timecode, or None if the parse fails
        """
        return self.__parseValAsScalar(ParamType.TIME)


    def __parseValAsVector(self, paramType):
        try:
            if isinstance(self.__value, ArrayNode):
                return self.__value.parseAsFlattened(paramType)
            else: # XXX check length here too?
                return paramType.parseMultipleValues(self.__value)
        except Exception as e:
            Logger.error("OSFI::Parameter. Error parsing value as a vector of {0}. {1}"
                        .format(str(paramType), str(e)))
            Logger.warning("OSFI::Parameter. Returning empty vector.")
            return []

    def getVectorBoolean(self):
        """
         <p>Gets the value of the parameter as a vector of booleans.</p>
        @return Parameter value as a vector of booleans.
        """
        return self.__parseValAsVector(ParamType.BOOLEAN)

    def getVectorString(self):
        """
         <p>Gets the value of the parameter as a vector of strings.</p>
         @return Parameter value as a vector of strings.
        """
        return self.__parseValAsVector(ParamType.STRING)

    def getVectorFile(self):
        """
         <p>Gets the value of the parameter as a vector of file names.</p>
         @return Parameter value as a vector of file names.
        """
        return self.__parseValAsVector(ParamType.FILE)

    def getVectorDouble(self):
        """
         <p>Gets the value of the parameter as a vector of double values.</p>
        @return Parameter value as a vector of double values.
        """
        return self.__parseValAsVector(ParamType.FLOAT)

    def getVectorInt(self):
        """
         <p>Gets the value of the parameter as a vector of integer values.</p>
        @return Parameter value as a vector of integer values.
        """
        return self.__parseValAsVector(ParamType.INTEGER)

    def getVectorTime(self):
        """
         <p>Gets the value of the parameter as a vector of timecodes.</p>
        @return Parameter value as a vector of timecode values.
        """
        return self.__parseValAsVector(ParamType.TIME)
    

    def __parseValAsMatrix(self, paramType):
        """
        Helper function that does the heavy lifting to parse a value as a matrix of some type.
        """
        try:
            if self.__paramComplexType == ParamComplexType.ARRAY:
                raise TypeError("Not valid for Parameter with ARRAY type")
            
            dims = self.getDims()
            if (len(dims) != 2):
                raise TypeError("Parameter is not a matrix")
            columns, rows = dims # NOTE: "swapped" (cols goes first) due to the E2E-ICD spec!
            if rows == 0:
                raise ValueError("Parameter dims are not valid")
            matrix = [([None] * columns) for i in range(rows)]

            tokens = paramType.splitValues(self.__value)
            if len(tokens) != rows * columns:
                raise ValueError("Expected {0}*{1}={2} elements, got {3}"
                                 .format(rows, columns, rows * columns, len(tokens)))
            for i in range(0, rows):
                for j in range(0, columns):
                    n = i * columns + j
                    val = self.__parseSingleValAs(paramType, tokens[n])
                    if val is None: # Error with the wrong value logged by __parseSingleValAs
                        raise ValueError()
                    matrix[i][j] = val
            return matrix
        except Exception as e:
            Logger.error("OSFI::Parameter. Error parsing value as a matrix of {0}. {1}"
                        .format(str(paramType), str(e)))
            Logger.warning("OSFI::Parameter. Returning empty matrix.")
            return []

    def getMatrixDouble(self):
        """
         <p>Gets the value of the parameter as a dynamic matrix of double values.</p>
         <p>Not valid for ARRAY type.</p>
         @return Value as matrix of doubles.
        """
        return self.__parseValAsMatrix(ParamType.FLOAT)

    def getMatrixInt(self):
        """
         <p>Gets the value of the parameter as a dynamic matrix of integer values.</p>
         <p>Not valid for ARRAY type.</p>
         @return Value as matrix of integers.
        """
        return self.__parseValAsMatrix(ParamType.INTEGER)

    def getMatrixString(self):
        """
         <p>Gets the value of the parameter as a dynamic matrix of string values.</p>
         <p>Not valid for ARRAY type.</p>
         @return Value as matrix of string.
        """
        return self.__parseValAsMatrix(ParamType.STRING)

    def getMatrixBoolean(self):
        """
         <p>Gets the value of the parameter as a dynamic matrix of boolean values.</p>
         <p>Not valid for ARRAY type.</p>
         @return Value as matrix of booleans.
        """
        return self.__parseValAsMatrix(ParamType.BOOLEAN)

    def getMatrixFile(self):
        """
         <p>Gets the value of the parameter as a dynamic matrix of file name values.</p>
         <p>Not valid for ARRAY type.</p>
         @return Value as matrix of file names.
        """
        return self.__parseValAsMatrix(ParamType.FILE)
    
    def getMatrixTime(self):
        """
         <p>Gets the value of the parameter as a dynamic matrix of timecode values.</p>
         <p>Not valid for ARRAY type.</p>
         @return Value as matrix of timecodes.
        """
        return self.__parseValAsMatrix(ParamType.TIME)


    def getArrayValue(self, *indices, **kwargs):
        """Returns the parsed value of an ARRAY parameter as a structure of nested lists.
        Each dimension corresponds to a level of nesting, and lists may have different
        sizes due to the inhomogeneous nature of ARRAY parameters. If positional arguments
        are provided, they are indices into the parameter, so the following:
            v = p.getArrayValue(1,0)
        is equivalent to:
            v = p.getArrayValue()[1][0]
        with the exception that in the former call only the relevant part of the tree is parsed.
        """
        pt = kwargs.pop('asType', None)
        pt = pt if pt is not None else self.__type # Non-present or None: use declared type
        if kwargs:
            raise ValueError('Unrecognized kwargs: ' + ", ".join(str(k) for k in kwargs.keys()))
        if not isinstance(self.__value, ArrayNode):
            raise TypeError("Not an ARRAY parameter: " + self.__name)
        node = self.__value
        idxIntoData = None
        for iDim, idx in enumerate(indices):
            if node.children is not None: # Indexing into internal node: replace node to be parsed
                node = node.children[idx]
            elif iDim == len(indices) - 1:  # Indexing into the data node: allow if last index
                idxIntoData = idx # Save for later, since before parsing we cannot index into the string
            else:
                raise ValueError("Invalid indices: '{0}'[{1}]".format(self.__name, ",".join(str(d) for d in indices)))
        
        val = node.parseAs(pt)
        if idxIntoData is not None: # If necessary, apply the last index after parsing
            val = val[idxIntoData]
        return val

    def fileExists(self):
        """
         <p>Check the existence of the files specified within a FILE parameter</p>
         <p>Not valid for ARRAY type.</p>
        Note: 
        When the parameter value is a single file name, the boolean is in the first pos [0][0]
        When the parameter value is a vector of file names, the array is a single row one in horizontal 
        [0][0:vector.size]
        @return A matrix with a boolean for each file TRUE -> file found, FALSE ->file not found
                A file is considered found when it is correctly located and ready to be read.
        """
        if self.__paramComplexType == ParamComplexType.ARRAY:
                Logger.error("OSFI::Parameter::fileExists. Not valid for Parameter with ARRAY type")
                return None
            
        dims = self.getDims()
        if not dims: # scalar
            filenames = [[self.getFileValue()]]
        elif len(dims) == 1:
            filenames = [self.getVectorFile()]
        else:
            filenames = self.getMatrixFile()
        result = []
        for i in range(0, len(filenames)):  # number of rows
            result.append([])
            for j in range(0, len(filenames[0])):  # number of columns
                result[i].append(Parameter.FileExists(filenames[i][j]))
        return result


    @staticmethod
    def getFileFormattedValue(fileName):
        """
         <p>Gets the value of the parameter as an array of characters.</p>
         <p>Checks the fileName validity</p>
         @param fileName String where fileName path name is stored
         @return Parameter value as array of characters.
        """
        return parseSingleFile(fileName)
        

    def setValue(self, aValue):
        """
         <p>Sets the parameter value as a string.</p>
         @param aValue Parameter value.
        """
        self.__value = aValue

    def getNdims(self):
        """
         * <p>Gets the number of dimensions.</p>
         * @return Dimensions number.
        """
        return len(self.__dims)


    def getDims(self):
        """
         <p>Gets a vector with the size of the dimensions.</p>
         @return Vector with dimension sizes.
        """
        return self.__dims

    def getName(self):
        """
         <p>Gets the name of the parameter.</p>
         @return string with the parameter name.
        """
        import warnings

        warnings.warn(" Ambiguous function, instead use getLocalName or getPath", DeprecationWarning)
        return self.getPath()

    def getLocalName(self):
        """
         <p>Gets the name of the parameter.</p>
         @return string with the parameter name.
        """
        return self.__name[self.__name.rindex(".") + 1: ]

    def getPath(self):
        """
         <p>Gets the full path of the parameter.</p>
         @return string with the parameter full path and name.
        """
        return self.__name

    def getDescription(self):
        """
         <p>Gets a brief description of the parameters.</p>
         @return string with the description.
        """
        return self.__description


    def getUnits(self):
        """
         <p>Gets the units of the parameters.</p>
         @return string with the minimum value.
        """
        return self.__units


    def getMin(self):
        """
         <p>Gets the minimum value of the parameters.</p>
         @return string with the minimum value.
        """
        return self.__min


    def geMax(self):
        """
        @deprecated see #getMax()
        """
        return self.getMax()
    
    def getMax(self):
        """
         <p>Gets the maximum value of the parameters.</p>
         @return string with the maximum value.
        """
        return self.__max


    def getType(self):
        """
        <p>Gets the type of the parameter.</p>
         @return String with the type of the parameter.
        """
        return ParamType.getParamTypeName(self.__type)

    def getParamType(self):
        """Deprecated in favour of the uniform name getElementType"""
        return self.getElementType()
    def getElementType(self):
        """
        Returns the type of individual components of this parameter, as a ParamType instance.
        See also getType.
        """
        return self.__type


    def getParamComplexType(self):
        """
        <p>Gets the complex type of the parameter.</p>
         @return Integer with the type of the parameter.
         """
        return self.__paramComplexType
    
    def cleanString(self, string):
        """
        <p>Replaces "," and multiple empty spaces with one space.</p>
        @return string.
        """
        return ' '.join(string.replace(",", " ").split())
    
    def splitValues(self, string):
        """
        Split values by their delimiter ("'" for String types and spaces for the other types)
        @return String
        """

        string = self.cleanString(string)
        return self.__type.splitValues(string)
    

    @staticmethod
    def FileExists(filename):
        """
        <p> Checks if a file is found within the file system </p>
         @param filename
         @return True if file exists 
        """
        return os.path.exists(filename) 
    
    # OPENSF-AN-010: Matrix  representation
 
 
class ParamComplexType():
    """
     Internal definition of parameter complex types
    """
    NONE = 0
    ARRAY = 1
    MATRIX = 2
    
    @staticmethod
    def getParamComplexTypeName(paramComplexType):
        lookFor = None

        for member in dir(ParamComplexType):
            if (getattr(ParamComplexType, member) == paramComplexType):
                lookFor = member
        return lookFor
    
    @staticmethod
    def valueOf(name):
        return getattr(ParamComplexType, name)

class ArrayNode(object):
    """Represents a structure of unparsed values in a parameter of type ARRAY.
    It may contain either a list of children ArrayNodes; or a tuple of the
    declared size and unparsed text of the leaf node, which may be None for
    zero-sized nodes."""

    __slots__ = ['elems', 'children']
    
    def __init__(self, val):
        if isinstance(val, tuple) and len(val) == 2 and isinstance(val[0], int): # Leaf node
            self.elems = val # Tuple of (declared length, unparsed string)
            self.children = None
        elif isinstance(val, list): # Non-leaf node, elements should be ArrayNode too
            self.elems = None
            self.children = val
        else:
            raise ValueError("Invalid structure input for ARRAY node")
    
    def parseAs(self, paramType, indices=()):
        makePathId = lambda: "[{0}]".format(",".join(str(i) for i in indices))
        try:
            if self.elems: # Leaf node, parse values
                return paramType.parseMultipleValues(self.elems[1], expectedLen=self.elems[0])
            else: # Internal node, parse sub-elements
                return [s.parseAs(paramType, indices + (iS,)) for iS, s in enumerate(self.children)]
        except Exception as e:
            Logger.error("OSFI::Parameter. Error parsing value as an array node of {0}, at node {1}. {2}"
                        .format(str(paramType), makePathId(), str(e)))
            Logger.warning("OSFI::Parameter. Returning empty node.")
            return []
    
    def parseAsFlattened(self, paramType):
        vals = self.parseAs(paramType)
        # Flatten the array structure, losing the information about which 
        # specific elements are "missing" from the rectangular envelope
        def flatiter(parsedNode):
            if not isinstance(parsedNode, list):
                yield parsedNode # Leaf value
            else: # Depth-first flatten
                for n in parsedNode:
                    for v in flatiter(n): # yield from is Py3
                        yield v
        return list(flatiter(vals))
