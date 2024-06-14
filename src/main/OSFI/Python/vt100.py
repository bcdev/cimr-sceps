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
----------------------------------------------
 Terminal functions to use in console programs
----------------------------------------------
"""

import sys

VT_BLACK = "0"
VT_RED = "1"
VT_GREEN = "2"
VT_YELLOW = "3"
VT_BLUE = "4"
VT_MAGENTA = "5"
VT_CYAN = "6"
VT_WHITE = "7"
VT_DEFAULT = "9"


def setColors(fg, bg):
    sys.stdout.write("\33[3" + fg + ";4" + bg + "m")

def finalize():
    """
    Sets the default attributes. Do this before your program ends.
    """
    default_attributes = "\33[0m"
    sys.stdout.write(default_attributes)
