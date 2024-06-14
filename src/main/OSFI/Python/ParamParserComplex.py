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

from __future__ import print_function

try:
    from OSFI.Parameter import ArrayNode
    from OSFI.ParamReader import ParamParsingError, ParamReader
except ImportError:
    from Parameter import ArrayNode
    from ParamReader import ParamReader, ParamParsingError


class ParamParserComplex():
    """
     This class parses XML elements of the type 'parameter' that have a type ARRAY or MATRIX
    """

    def __init__(self, elem, completeName):
        try:
            self.__param = elem
            self.__fqName = completeName
            self.__complexType = elem.attrib["type"]
            self.type = self.__determineConsistentType()

            if self.__complexType == 'ARRAY':
                self.value = self.__readArray(self.__param, ())
            elif self.__complexType == 'MATRIX':
                self.value = self.__readMatrix()
            else:  # This shouldn't happen unless ParamReader and ParamParserComplex become desynced
                raise ParamParsingError(self.__fqName, "Invalid structured type " + self.__complexType)
        except ParamParsingError:
            raise  # Path already noted
        except Exception as e:
            raise ParamParsingError(self.__fqName, str(e))

    def __determineConsistentType(self):
        elType = None
        for elem in self.__param.iter('parameter'):  # iter() goes up to any depth
            curElType = elem.get('elementType')
            if curElType is None:  # No elementType attribute
                continue
            elif elType is None:  # First found, store
                elType = curElType
            elif elType != curElType:  # Then ensure they are all the same
                raise ParamParsingError(self.__fqName, "Inconsistent element types")
        if elType is None:
            raise ParamParsingError(self.__fqName, "Could not determine element type")
        return elType

    def __readMatrix(self):
        """Reads a MATRIX-typed parameter and returns a string corresponding to the
        concatenation of the contents of the rows, separated by spaces. A MATRIX param
        thus ends up with the same "value" as if it had been written in the old format.
        """
        # print("reading matrix " + self.__fqName)
        nC, nR = ParamReader._parseDims(self.__param.get('dims', ''),
                                        expectedRank=2)  # E2E-ICD: columns go first in a matrix
        rows = self.__param.findall('parameter')  # Not recursive, only direct children
        if len(rows) != nR:
            raise ParamParsingError(self.__fqName, "Number of rows is not the expected ({0} vs {1})"
                                    .format(len(rows), nR))

        value = ''
        for iR, r in enumerate(rows):
            makeRowPath = lambda: self.__fqName + '[' + str(iR) + ']'  # Create if needed
            if r.get('type') != 'ARRAY':
                raise ParamParsingError(makeRowPath(), "Row does not have the expected type ARRAY")
            rowDim = r.get('dims')
            # Rows can omit the dims (assumed to be nC), but if they do specify it, it must match nC
            if rowDim is not None and ParamReader._parseDims(rowDim, expectedRank=1) != nC:
                raise ParamParsingError(makeRowPath(), "Number of columns changed in row vs matrix")
            if r.findall('parameter'):  # Contains sub-elements!
                raise ParamParsingError(makeRowPath(), "Row contains sub-parameter elements instead of data")
            data = r.text
            # TODO split data as the target type to actually verify len against nC? However, this
            # behaviour (fail-fast, preventing whole file parsing) does not match that of old vectors.
            value += ' ' + data
        return value

    def __readArray(self, curElem, curPath):
        """
        Reads an ARRAY-typed parameter and sets value to a structure of nested lists representing
        the *unparsed* parameter values. Leaf nodes are represented by tuples (dim, raw-text).
        """
        curFq = self.__fqName + "".join('[{0}]'.format(idx) for idx in curPath)
        # print("reading array element " + curFq)
        try:
            dim = ParamReader._parseDims(curElem.get('dims', ''), expectedRank=1)  # Returns number, not tuple!
            subElems = curElem.findall('parameter')  # Not recursive, only direct children

            if not subElems:  # Leaf node, read values and split
                # value = ParamType.valueOf(self.type).splitValues(curElem.text)
                # TODO split data as the target type to actually verify len against dim? However, this
                # behaviour (fail-fast, preventing whole file parsing) does not match that of old vectors.
                value = ArrayNode((dim, curElem.text))
            else:  # Non-leaf, parse subelements
                if len(subElems) != dim:
                    raise ParamParsingError(curFq, "Number of sub-elements is not the expected ({0} vs {1})"
                                            .format(len(subElems), dim))
                value = ArrayNode([self.__readArray(s, curPath + (iS,)) for iS, s in enumerate(subElems)])
            return value
        except ParamParsingError:
            raise  # Path already noted
        except Exception as e:
            raise ParamParsingError(curFq, str(e))
