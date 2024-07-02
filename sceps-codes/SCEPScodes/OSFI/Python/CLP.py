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

import sys
import argparse

try:
    from . import Logger
    from .UsageReader import UsageReader
except ImportError:
    import Logger
    from UsageReader import UsageReader


class CLP():
    """
    CLP class provide routines to parse the command line arguments, generating
    a list of tokens grouped as:
    - Configuration file
    - Input files
    - Output files

    Command line must follow the following format (expressed in Extended Backus-Naur form):

        <command_line> ::= <executable_name> <whitespaces> <configuration_files> <whitespaces>
                           <input_files> <whitespaces> <output_files> <EOL>
        <executable_name> ::= <file_name>
        <whitespace>  ::= (“ “)
        <whitespaces> ::= <whitespace>+
        <file_name>   ::= (<alphanumeric>)+
        <configuration_files> ::= <list_of_files>
        <input_files>   ::= <list_of_files>
        <output_files>  ::= <list_of_files>
        <list_of_files> ::= <file_name> (“,” <file_name>)*

    In order to increase the flexibility of the OpenSF/Modules CLI, the following new CLI should be adopted (CLP-v2):

        <command_line> ::= <executable> <options>
        <options> ::= (<global_configuration> | <local_configuration> | <input_file> | <output_file>)*
        <executable>> ::= <file_name>
        <global_configuration> ::= (--global|-g) <file_name>
        <local_configuration> ::= (--local|-l) <file_name>
        <input_file> ::= (--input|-i) <file_name>
        <output_file> ::= (--output|-o) <file_name>
        <file_name> ::= (<alphanumeric>)+

    File names must be valid OS-dependant file locations. Per the E2E-ICD, the configuration, input
    and output file names can be either a full-path name or a relative path name. Relative path names
    are interpreted by the module as relative to the current working directory. In particular, no
    reference to the possible E2E_HOME variable is made.
    """

    """
     Character for delimiting file names in a list
    """
    __delimiter = ","


    def __init__(self, argv=None):
        """
        Constructor, may be given the command line to parse or take it from sys.argv.
        Note that the command line is expected to have a first element (argv[0]) representing the
        application name, and it is dropped while parsing. Thus, a list of length 3 or 4 is needed,
        where the last element (argv[3], representing the output files) is optional.
        """

        # Initialise attributes
        self.__confFile = "" # Legacy value: location of the _single_ configuration file.
        self.__confFiles   = [] # List of locations of the configuration files.
        self.__inputFiles  = [] # List of locations of the input files.
        self.__outputFiles = [] # List of locations of the output files.
            
        # Command line options. If not provided specifically, use the command line
        if argv is None:
            argv = sys.argv

        if len(argv) >= 2 and argv[1].startswith('-'):  # Assuming CLP-v2
            for arg in argv[1::2]:
                if not arg.startswith('-'):
                    self.__raiseErrorMixingOptions()
            self.__parseArgsV2(argv)

        else:  # Assuming CLP-v1
            for arg in argv:
                if arg.startswith('-'):
                    self.__raiseErrorMixingOptions()

            if len(argv) == 4: # Three arguments: conf input output
                self.__confFiles, self.__inputFiles, self.__outputFiles = \
                    (self.__parseFilesV1(x) for x in argv[1:])
                self.__confFiles = [None if file == "" else file for file in self.__confFiles]  # Converting to None if
                #                                        configuration file is empty string, to be in line with v2 CLP
            elif len(argv) == 3: # Two arguments: input output (no configuration files, LEGACY for pre-E2E-ICD modules)
                Logger.warning("Calling OSFI modules with only input and output files is deprecated");
                self.__inputFiles, self.__outputFiles = (self.__parseFilesV1(x) for x in argv[1:])
            else: # Otherwise -> error
                Logger.error("Program needs two or three arguments")
                try:
                    appName = argv[0]
                    usageReader = UsageReader(appName + ".usage.xml")
                    usage = usageReader.getUsage()
                    print(usage, file=sys.stderr)
                except Exception:
                    pass # Ignore errors showing usage
                raise ValueError("Wrong number of elements to CLP argv: expected ['appname', 'confs', 'inputs', 'outputs']")
            self.__confFile = argv[1] # Legacy, single conf file support


    def __parseArgsV2(self, argv):
        """
         Parses the arguments received on the CLI, for the v2 CLP version. This version allows to give the arguments
         without required order, preceded by a string stating the argument option that is being given.
        """

        # Defining arguments
        parser = argparse.ArgumentParser(description="CLI parser of module parameters")
        parser.add_argument("module", action="store", default=None, help="Module to be run with the next arguments as "
            "inputs", type=str)
        parser.add_argument("-g", "--global", action="store", default=None, help="Global configuration file path. Only "
            "one needed. If several are given, only the last one will be considered", type=str)
        parser.add_argument("-l", "--local", action="store", default=None, help="Local configuration file path. Only "
            "one needed. If several are given, only the last one will be considered", type=str)
        parser.add_argument("-i", "--input", action="append", default=[], help="Input(s) file path(s)", type=str)
        parser.add_argument("-o", "--output", action="append", default=[], help="Output(s) file path(s)", type=str)

        # Parsing arguments
        arguments = parser.parse_args(argv)
        self.__confFiles = [getattr(arguments, 'global'), arguments.local]
        self.__inputFiles = arguments.input
        self.__outputFiles = arguments.output

    def __raiseErrorMixingOptions(self):
        """
         Used to raise error when arguments given to CLP follow patterns mixing the v1 and v2 CLP versions.
        """
        Logger.error("Mixing options and positional arguments in CLP")
        raise ValueError(
            "Wrong value of elements to CLP argv: expected positional arguments "
            "['appname', 'confs', 'inputs', 'outputs'] or options ['--global', 'g', '--local', 'l', "
            "'--input', 'i1', '--input', 'i2', '--input', 'i3', '--output', 'o1', '--output', 'o2']")

    def __parseFilesV1(self, argv):
        """
         For the CLP v1 version. Parses a string to obtain a list of tokenized strings, separated by a delimiter.
         Throws an exception in case of invalid strings.
        """
        Logger.debug("OSFI::CLP::parseFiles. Parsing " + argv)
        if not argv: # Fully empty string
            return []
        result = argv.split(CLP.__delimiter)
        for arg in result:
            self.__checkValidFile(arg)
        
        return result
            
    def __checkValidFile(self, fileName):
        """
         Checks the OS validity of a file name.
         @param file - file name to validate
         @throw exception in case of invalid file.
         """
        import re
        # Test to see if the file is ok.
        invalid = ","
        if re.search("[" + re.escape(invalid) + "]", fileName):
            msg = "OSFI::CLP.checkValidFile. File '" + fileName + "' is not a valid OS file location"
            Logger.error(msg)
            raise ValueError(msg)

    def getConfFile(self):
        """
         Gets the name of the configuration file.
         Configuration file can be "" is no configuration file is passed.
        """
        return self.__confFile


    def getConfFiles(self):
        """
         Gets a copy of the list of the configuration files.
         Update: for CLP-v2, this list will always have two elements, with the first one representing the global
         configuration file and the second one the local configuration file. If any of them is not given in the command
         line, the corresponding slot will be None.
        """
        return self.__confFiles[:]

    def getInputFiles(self):
        """
         Gets a copy of the list of the input files.
        """
        return self.__inputFiles[:]

    def getOutputFiles(self):
        """
         Gets a copy of the list of the output files.
        """
        return self.__outputFiles[:]
