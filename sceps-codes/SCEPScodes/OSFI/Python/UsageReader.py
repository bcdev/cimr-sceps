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

import xml.etree.ElementTree as ET


class UsageReader:
    """
     This class is responsible to read a XML file showing the usage instructions of a
     certain executable in case the Command Line Arguments Parser (CLP) detects a
     failure in the executable invocation.
    """

    def __init__(self, xmlFile):
        """
         Class constructor.
         @param xmlFile Location of the XML file with the usage description.
         @throws exception if an error occurred while parsing the file
        """
        self.__usage = ""
        
        tree = ET.parse(xmlFile)
        root = tree.getroot()
        
        self.__usage = root[0].text
        
    def getUsage(self):
        """
         Returns a string with the usage description.
        """
        return self.__usage
