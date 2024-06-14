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

"""
Logger module provides routines for presenting messages in the standard output.
Any E2E-ICD-complieant framework can later intercept these messages.

Logger formats the logging messages as follows:

    <message>   ::= (<progress> | <log>) <EOL>
    <progress>  ::= "Progress" <whitespaces> <delimiter> <whitespace> <progress_body>
    <delimiter> ::= "|"
    <progress_body> ::= <integer> " of " <integer>
    <log> ::= <type> <whitespaces> <delimiter> <whitespaces> <text> [<whitespaces>
              <delimiter> <whitespaces> <version>]
    <type>    ::= "Error" | "Warning" | "Info" | "Debug"
    <version> ::= <digit>("." <digit>)*
    <digit>   ::= "0" | "1" | "2" | "3" | "4" | "5" | "6" | "7" | "8" | "9"
    <whitespace>  ::= (" ")
    <whitespaces> ::= <whitespace>+

This format defines five different types of messages:
- Information. This is a message raised by the model describing a harmless event,
    giving knowledge to the user. Model execution should continue with no interruptions.
- Warning. The model has detected a non-fatal error or situation that may cause a
    fatal error. This is a harmless event, thus, the execution should continue with no
    interruption.
- Error. A fatal error has happened in the model execution and the model itself
    informs the user about it, so the model has time to gracefully close the execution.
- Debug. Detailed information of the model execution given to the user. Information
    is intended to lead the user (or model developer) upon fixing a problem. This is a
    harmless event so model execution should continue with no interruptions. Debug messages
    are only shown if an environment variable named DEBUG_MODE is defined and set as On in
    the model execution context.
- Progress. Numerical information on the amount of model execution performed.

Logger also provides a way to finish the model execution (after a fatal error or at
the expected end of the execution) and return a non-zero code to the operating system.
"""

import os
import sys

try:
    from OSFI.vt100 import setColors, finalize, VT_DEFAULT, VT_RED, VT_YELLOW, VT_GREEN, VT_CYAN
except ImportError:
    from vt100 import setColors, finalize, VT_DEFAULT, VT_RED, VT_YELLOW, VT_GREEN, VT_CYAN

_delimiter = "|"  # Character delimiting the type and the body of the message.

# OPENSF-AN-005: Log messages in OSFI should include the library version
def getVersion():
    if getVersion.osfiVer is None:
        try: from OSFI._osfi_version import __version__
        except ImportError:
            try: from _osfi_version import __version__
            except ImportError:
                import warnings
                warnings.warn('No _osfi_version module found to be imported! Running OSFI from the development tree?')
                __version__ = '(unknown)'
        getVersion.osfiVer = __version__
    return getVersion.osfiVer
getVersion.osfiVer = None # Stores the OSFI version from the CMake-substituted module


def isDebugging():
    """
    Checks if the system is in debug mode. Environment variable
    DEBUG_MODE controls this setting, with "On" enabling debug output.
    """
    if isDebugging.debugOutput is None:
        if 'DEBUG_MODE' in os.environ:
            isDebugging.debugOutput = (os.environ['DEBUG_MODE'] == "On")
        else:
            warning("OSFI::EHLog. DEBUG_MODE variable is not initialized. Setting to Off")
            isDebugging.debugOutput = False
    return isDebugging.debugOutput
isDebugging.debugOutput = None # Stores the debug mode.


def isColored():
    """
     Checks if the system is in color mode.
      OSFI_LOG_COLOR if ON the OSFI logs are colored.
     - Error : Red
     - Warning: Yellow
     - Progress: Cyan
     - Info: Green
    """
    if isColored.colorEnabled is None:
        # Read value from environment variable, defaulting to no colour
        isColored.colorEnabled = (os.environ['OSFI_LOG_COLOR'] == "On") \
            if 'OSFI_LOG_COLOR' in os.environ else False
    return isColored.colorEnabled


isColored.colorEnabled = None  # Stores the OSFI log colour mode.


def _writeLogFlair(flair, fg_color):
    """Internal utility function that writes only the log flair, padded to the right size."""
    padded = "{: <9s}".format(flair)
    if fg_color is not None and isColored():
        setColors(fg_color, VT_DEFAULT)
        sys.stdout.write(padded)
        finalize()
        sys.stdout.write(_delimiter + " ")
    else:
        sys.stdout.write(padded + _delimiter + " ")
    return sys.stdout

def _writeLogMsgWithFlair(flair, fg_color, message, withVer = False):
    """Internal utility function that writes a full log message, with the given flair."""
    stream = _writeLogFlair(flair, fg_color)
    stream.write(message)
    if withVer:
        stream.write(_delimiter + getVersion())
    stream.write('\n')


def error(message):
    """
      Shows a formatted error message.
      @param message - text of the message
    """
    _writeLogMsgWithFlair("Error", VT_RED, message)

def warning(message):
    """
      Shows a formatted warning message.
      @param message - text of the message
    """
    _writeLogMsgWithFlair("Warning", VT_YELLOW, message)

def info(message):
    """
      Shows a formatted info message.
      @param message - text of the message
    """
    _writeLogMsgWithFlair("Info", VT_GREEN, message)

def debug(message):
    """
     Shows a formatted debug message, if the environment variable "DEBUG_MODE" is
     equal to "On".
     First time this function is called, checks the declaration of that variable. If
     it is not declared, presents a warning message and assumes it as "Off".
      @param message - text of the message
    """
    if isDebugging():
        _writeLogMsgWithFlair("Debug", None, message)


def getErrorStream():
    """
     Returns the output stream used to show error messages. Developers can
     use this function as a convenient method for showing complex data types.
     """
    sys.stdout.flush()
    return _writeLogFlair("Error", VT_RED)

def getDebugStream():
    """
     Returns the output stream used to show debugging messages. Developers can
     use this function as a convenient method for showing complex data types.
     Returns a "null stream" is not in debug mode, implying that nothing is written
     on output stream.     
    """
    if isDebugging():
        sys.stdout.flush()
        return _writeLogFlair("Debug", None)
    else:
        return os.devnull

def getInfoStream():
    """
     Returns the output stream used to show info messages. Developers can
     use this function as a convenient method for showing complex data types.
     """
    sys.stdout.flush()
    return _writeLogFlair("Info", VT_GREEN)

def getWarningStream():
    """
     Returns the output stream used to show warning messages. Developers can
     use this function as a convenient method for showing complex data types.
     """
    sys.stdout.flush()
    return _writeLogFlair("Warning", VT_YELLOW)


def progress(step, nSteps):
    """
     Shows a formatted progress message.
     @param step - current step number
     @param nSteps - maximum number of steps.
    """
    _writeLogMsgWithFlair("Progress", VT_CYAN, str(step) + " of " + str(nSteps))


def finishExecution(errorCode):
    """
    Shows an information message and exits the program execution with
    an specific error code.
    @param errorCode - Code of error to exit with.
    """
    info("Finishing module execution")
    exit(errorCode)


def qualityReport(name, value):
    """
     Shows a formatted message reporting a quality indicator
     @param name - Name of the quality indicator
     @param value - Double associated to the quality indicator pointed by name
    """
    return _writeLogMsgWithFlair("Quality", None, name + ": " + str(value))


def getQualityStream():
    """
     Returns the output stream used to show quality messages. Developers can
     use this function as a convenient method for showing complex data types.
     """
    sys.stdout.flush()
    return _writeLogFlair("Quality", None)
