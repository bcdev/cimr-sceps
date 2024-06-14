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

import re


def checkFieldRange(name, val, min, max):
    """
    Private helper for checking the individual values in the TimeValue ctor
    """
    if val < min or val > max:
        raise ValueError("Invalid {0}: {1}".format(name, val))
    return val



def daysInProlepticGregorianMonth(year, month):
    """
    Helper for the days in a month on an "eternally" extended Gregorian calendar (in years 0000-9999)
    """
    if month in (1, 3, 5, 7, 8, 10, 12):
        return 31
    elif month in (4, 6, 9, 11):
        return 30
    elif month == 2:
        isGrLeap = (year % 4 == 0) and ((year % 100 != 0) or (year % 400 == 0))
        return 29 if isGrLeap else 28
    raise AssertionError("Invalid month in internal function call!")

class TimeValue(object):
    """
    Class that represents the parsed value of a TIME parameter. It holds attributes for the
    fields of a CCSDS ASCII "A" variant time code.
    """
    __slots__ = ('year', 'month', 'dom', 'hour', 'minute', 'sec', 'nanosec')

    def __init__(self, year, month, dom, hour, minute, sec, nanosec):
        # Limits specified by the CCSDS definition of time code "A"
        self.year    = checkFieldRange("year", year, 1, 9999)
        self.month   = checkFieldRange("month", month, 1, 12)
        self.dom     = checkFieldRange("day-of-month", dom, 1, daysInProlepticGregorianMonth(year, month))
        self.hour    = checkFieldRange("hour", hour, 0, 23)
        self.minute  = checkFieldRange("minute", minute, 0, 59)
        self.sec     = checkFieldRange("sec", sec, 0, 60) # Allow "60" in leap seconds, which we don't check
        self.nanosec = checkFieldRange("nanosecond", nanosec, 0, 999999999)

    def __format__(self, format):
        """
        Format this object for output. The only accepted format is 'a', the CCSDS ASCII "A" time code variant A.
        """
        if format == 'a': # CCSDS ASCII time code "A": YYYY-MM-DD"T"HH:mm:SS.f
            # Compute the value and width of the "fraction-of-second" field
            secFracVal = self.nanosec
            secFracWidth = 9 # Prevent trailing zeros while keeping leading ones
            while secFracWidth > 3 and secFracVal % 10 == 0:
                secFracWidth -= 1
                secFracVal //= 10
            return "{0.year:04d}-{0.month:02d}-{0.dom:02d}T{0.hour:02d}:{0.minute:02d}:{0.sec:02d}.{1:0{2:d}d}" \
                .format(self, secFracVal, secFracWidth)
        raise ValueError("Invalid format specifier '{0}'".format(format))
    
    def __repr__(self):
        return self.__format__('a')
    
    def __astuple(self):
        """Returns this object's fields as a tuple. May disappear if TimeValue is changed to extend NamedTuple."""
        return (self.year, self.month, self.dom, self.hour, self.minute, self.sec, self.nanosec)
    
    def __eq__(self, value): # __ne__ is autogenerated as "not __eq__"
        """Two instances are equal if all their field values are exactly equal.
        No conversions or adaptations are performed so month=4 and dom=0 is different from month=3 and dom=31."""
        if not isinstance(value, TimeValue):
            return NotImplemented
        return self.__astuple() == value.__astuple()
    
    def __lt__(self, value):
        """Two instances are ordered according to their fields, from most to least significant.
        No conversions or adaptations are performed so month=3 and dom=31 is less-than month=4 and dom=0."""
        if not isinstance(value, TimeValue):
            return NotImplemented
        return self.__astuple() < value.__astuple()
    
    # NOTE: in the future we will be able to remove all ops except for __eq__ and __lt__ and use functools.total_ordering
    def __le__(self, value):
        if not isinstance(value, TimeValue):
            return NotImplemented
        return self.__astuple() <= value.__astuple()
    def __ge__(self, value):
        if not isinstance(value, TimeValue):
            return NotImplemented
        return self.__astuple() >= value.__astuple()
    def __gt__(self, value):
        if not isinstance(value, TimeValue):
            return NotImplemented
        return self.__astuple() > value.__astuple()

    @staticmethod
    def parse(s):
        """
        Parse the given string in a recognized time format (CCSDS ASCII time code) into a TimeValue instance.
        Raises ValueError if the format is unrecognized or the field values are invalid
        """
        m = TIME_CCSDS_ASCII_A.match(s)
        matchMaD = m is not None # Day-of-month groups
        if not matchMaD: # Try with the day-of-year form
            m = TIME_CCSDS_ASCII_B.match(s)
        matchDoY = not matchMaD and m is not None # Day-of-year groups
        if not (matchMaD or matchDoY):
            raise ValueError("Bad format for TIME value: " + s)

        # Common fields
        iM_year = 1
        iM_hour = 3 if matchDoY else 4
        # The fit into each data type is forced by the number of digits (4, 2, etc.)
        year = int(m.group(iM_year))
        hour = int(m.group(iM_hour))
        min = int(m.group(iM_hour + 1))
        sec = int(m.group(iM_hour + 2))

        secfrac = m.group(iM_hour + 3)
        if secfrac is not None:
            sfrac_digits = len(secfrac)
            if sfrac_digits > 9: # We only support up to nanos, so "000000000".."999999999"
                raise ValueError("Second fractions supported only up to ns, read: " + secfrac)
            # Multiply the fraction by 10^(9-number of digits) to convert to nanos
            nanos = int(secfrac) * 10 ** (9 - sfrac_digits)
        else: # Fraction-of-second part unmatched, consider as ".0"
            nanos = 0

        # Parse either -month-day or -dayofyear
        if matchMaD:
            month = int(m.group(iM_year + 1))
            dom = int(m.group(iM_year + 2))
        else:
            # FIXME: not mentioned in the examples of the E2E-ICD v2.4,
            # but it refers to "a CCSDS ASCII time format", so we should
            # support the day-of-year format, which is time code "type B".
            raise NotImplementedError("Day-of-year format for TIME parameters not implemented")

        return TimeValue(year, month, dom, hour, min, sec, nanos) # Checks values

# CCSDS time code A format, as YYYY-MM-DD"T"HH:mm:SS.f+"Z"
TIME_CCSDS_ASCII_A = re.compile(r"(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2})(?:\.(\d+))?Z?")
# CCSDS time code B format, as YYYY-ddd"T"HH:mm:SS.f+"Z"
TIME_CCSDS_ASCII_B = re.compile(r"(\d{4})-(\d{3})T(\d{2}):(\d{2}):(\d{2})(?:\.(\d+))?Z?")
