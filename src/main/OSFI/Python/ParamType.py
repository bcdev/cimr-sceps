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

import os
import re
import pathlib

try:
    from . import Logger
    from .TimeValue import TimeValue
    from .OsfiEnv import get_env
except ImportError: # XXX Deprecation of non-package mode
    import Logger
    from TimeValue import TimeValue
    from OsfiEnv import get_env


class ParamType(object):
    """
     Internal definition of parameter types
     """
    __values = {} # Keep a central DB of all param types
    
    def __new__(cls, name, parseFun, is_str):
        if name in ParamType.__values:
            raise ValueError("ParamType with that name already exists")
        instance = super(ParamType, cls).__new__(cls)
        ParamType.__values[name] = instance
        return instance

    def __init__(self, name, parseFun, is_str):
        self.__name = name
        self.__parseSingleVal = parseFun
        self.__isStringBased  = is_str
    
    def __str__(self):
        return self.__name
    
    def parseSingleValue(self, x):
        return self.__parseSingleVal(x)
    
    def parseMultipleValues(self, value, expectedLen=None):
        """
        Parse a sequence of values in a string into a vector. If expectedLen is given,
        also assert that the number of parsed values equals that argument.
        """
        vector = []
        tokens = self.splitValues(value, expectedLen is not None and expectedLen == 0)
        if expectedLen is not None and len(tokens) != expectedLen:
            raise ValueError("Number of elements is different than expected: {0} vs {1}"
                    .format(len(tokens), expectedLen))
        for iT, token in enumerate(tokens):
            try:
                vector.append(self.__parseSingleVal(token))
            except Exception as e:
                raise ValueError("Invalid value at element {0} for type {1}: '{2}'"
                        .format(iT, self.__name, str(e)))
        return vector
    
    # TODO: this was moved from Parameter.splitValues, it needs to be cleaned up
    def splitValues(self, string, expected_empty=False):
        """
        Split values into the tokens to be parsed
        @param string unparsed str or None (only for empty dimensions in ARRAY parameters)
        @param expected_empty True if the data corresponds to an empty dimension in an ARRAY parameter
        @return List[str] with the unparsed but split values
        """

        if (string is None or not string.strip()) and expected_empty:
            return []
        elif self.__isStringBased:
            # String should be in the format ""'a' 'b c' 'd'". The old algorithm (which is kept with
            # the substitution of str methods for re methods) is to split by "' '" and then to
            # remove the initial and final quotes.
            splitValues = re.split(r"'\s+'", string) # Note that "''" is not a valid delimiter
            firstQuote = re.match(r"\s*'", splitValues[0]) # "Match" is like an automatic ^ in the RE...
            if firstQuote:
                splitValues[0] = splitValues[0][firstQuote.end():]
            else: # A single unquoted value, remove leading and trailing space
                splitValues[0] = splitValues[0].strip()
            lastQuote = re.search(r"'\s*$", splitValues[-1]) # ... so here we use search instead
            if lastQuote:
                splitValues[-1] = splitValues[-1][:lastQuote.start()]
            # XXX given that "''" is not a valid separator, should it be considered as escaping "'"
            # and thus replaced in the split tokens?
            return splitValues
        else:
            # String should be in the format "1 2 3". Here, direct usage of str.split(sep=None)
            # takes care of ignoring multiple separating spaces.
            return string.split()

    @staticmethod
    def getParamTypeName(paramType):
        return paramType.__name
    
    @staticmethod
    def valueOf(typeName):
        return ParamType.__values[typeName]

def parseSingleBool(x):
    """
    Parse a single value of type BOOLEAN, according to the E2E-ICD rules.
    """
    x = x.strip()
    if x == "TRUE":
        return True
    elif x == "FALSE":
        return False
    else:
        raise ValueError('Accepted values are TRUE or FALSE (case-sensitive), got "{0}"'.format(x))

def parseSingleFile(fileName):
    """
    Parse the value of a single item in a FILE or FOLDER parameter. If the value is an absolute path
    under the current platform, it is returned as-is, whereas a relative path is composed with the
    "base folder" determined from environment variables per the E2E-ICD.
    See filesBaseDir.
    """
    # Check whether the fileName has absolute path or relative
    f = pathlib.PurePath(fileName)
    if f.is_absolute():
        return str(f) # Also performs / -> \ slash conversion in Windows
    base = filesBaseDir()
    if not f.parts and not base.parts:
        return '' # Empty base and fileName: return '' instead of '.', keeping previous behaviour
    return str(base / fileName)

def filesBaseDir(env=None):
    """
    Returns the "base directory" used to resolve relative paths of ParamType.FILE and
    ParamType.FOLDER types. Following the E2E-ICD, it is read from the enviroment variable
    "E2E_HOME".	As a deprecated option, if that variable is not present, the "OSFI_HOME" variable
    is consulted next. If both are absent from the environment, then the current working directory is used.
    
    :param env: The default (None) uses the default environment, but it can be replaced (e.g. for
        unit tests) by passing it in.
    :return: Base directory used to resolve relative paths in parameters.
    :rtype: pathlib.PurePath
    """
    if env is None:
        env=get_env()
    ENVKEY_BASE_DIR = "Parameter.base_dir"
    stored = env[ENVKEY_BASE_DIR]
    if stored is not None:
        return stored
    
    e2e_home = env.external_var("E2E_HOME")
    legacy_home = env.external_var("OSFI_HOME")
    if e2e_home is not None:
        newBdir = pathlib.PurePath(e2e_home)
    elif legacy_home is not None:
        Logger.warning("OSFI: E2E_HOME is not set, but OSFI_HOME is. Note that OSFI_HOME is deprecated and may be removed.")
        newBdir = pathlib.PurePath(legacy_home)
    else:
        Logger.warning("OSFI: E2E_HOME is not set, using the current directory for relative paths.")
        newBdir = pathlib.PurePath('')
    env[ENVKEY_BASE_DIR] = newBdir
    return newBdir

#
# Define the user-visible parameter types
#
ParamType.INTEGER = ParamType("INTEGER", parseFun=int, is_str=False)
ParamType.FLOAT   = ParamType("FLOAT", parseFun=float, is_str=False)
ParamType.BOOLEAN = ParamType("BOOLEAN", parseFun=parseSingleBool, is_str=False)
ParamType.STRING  = ParamType("STRING", parseFun=lambda x: x, is_str=True) # No parsing required for strings
ParamType.FILE    = ParamType("FILE", parseFun=parseSingleFile, is_str=True)
ParamType.FOLDER  = ParamType("FOLDER", parseFun=parseSingleFile, is_str=True)
ParamType.TIME    = ParamType("TIME", parseFun=TimeValue.parse, is_str=False)
