function [year, month, day, hour, minute, second] = unixsecs2date(secs)
%UNIXSECS2DATE Number of seconds since 00:00:00 1 January 1970 to date.
%
%   [YEAR, MONTH, DAY, HOUR, MINUTE, SECOND] = UNIXSECS2DATE(SECS) returns
%   the Gregorian calendar date (year, month, day, hour, minute, and second)
%   corresponding to given number of seconds since 00:00:00 1 January 1970.
%
%   Any missing MONTH or DAY will be replaced by ones.  Any missing HOUR,
%   MINUTE or SECOND will be replaced by zeros.  If no date is specified,
%   the current date and time is used.
%
%   In UNIX, the smallest time unit is a signed 32-bit integer counting the
%   number of seconds since 00:00:00 1 January 1970.  The range is from
%   1901-12-13 20:45:52, when the number of seconds is 2^31-1, to 2038-01-19
%   03:14:07, when the number of seconds is 2^31.
%
%   This function is compatible but the number of seconds is not limited to
%   a 32-bit integer, any MATLAB double precision number may be used.  Also,
%   fractional seconds are allowed.
%
%   See also DATE2UNIXSECS.
%   $Id: unixsecs2date.m 8345 2013-04-17 18:16:40Z gerrit $

%   Author:      Peter J. Acklam
%   Time-stamp:  2003-01-14 21:32:11 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   narginchk(1, 1);

if isscalar(secs)
   dv = datevec(addtodate(datenum(1970, 1, 1, 0, 0, 0), double(secs), 'second'));
   C = num2cell(dv);
   [year, month, day, hour, minute, second] = deal(C{:});
else
   [year, month, day, hour, minute, second] ... 
	 	      = jd2date(double(secs) / 86400 + date2jd(1970, 1, 1)); 
end

return

function [year, month, day, hour, minute, second] = jd2date(jd)
%JD2DATE Gregorian calendar date from modified Julian day number.
%
%   [YEAR, MONTH, DAY, HOUR, MINUTE, SECOND] = JD2DATE(JD) returns the
%   Gregorian calendar date (year, month, day, hour, minute, and second)
%   corresponding to the Julian day number JD.
%
%   Start of the JD (Julian day) count is from 0 at 12 noon 1 JAN -4712
%   (4713 BC), Julian proleptic calendar.  Note that this day count conforms
%   with the astronomical convention starting the day at noon, in contrast
%   with the civil practice where the day starts with midnight.
%
%   Astronomers have used the Julian period to assign a unique number to
%   every day since 1 January 4713 BC.  This is the so-called Julian Day
%   (JD). JD 0 designates the 24 hours from noon UTC on 1 January 4713 BC
%   (Julian calendar) to noon UTC on 2 January 4713 BC.

%   Sources:  - http://tycho.usno.navy.mil/mjd.html
%             - The Calendar FAQ (http://www.faqs.org)

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-05-24 15:24:45 +0200
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   nargsin = nargin;
   error(nargchk(1, 1, nargsin));

   % Adding 0.5 to JD and taking FLOOR ensures that the date is correct.
   % Here are some sample values:
   %
   %  MJD     Date       Time
   %  -1.00 = 1858-11-16 00:00 (not 1858-11-15 24:00!)
   %  -0.75 = 1858-11-16 06:00
   %  -0.50 = 1858-11-16 12:00
   %  -0.25 = 1858-11-16 18:00
   %   0.00 = 1858-11-17 00:00 (not 1858-11-16 24:00!)
   %  +0.25 = 1858-11-17 06:00
   %  +0.50 = 1858-11-17 12:00
   %  +0.75 = 1858-11-17 18:00
   %  +1.00 = 1858-11-18 00:00 (not 1858-11-17 24:00!)

   ijd = floor(jd + 0.5);               % integer part

   if nargout > 3
      fjd = jd - ijd + 0.5;             % fraction part
      [hour, minute, second] = days2hms(fjd);
   end

   % The following algorithm is from the Calendar FAQ.

   a = ijd + 32044;
   b = floor((4 * a + 3) / 146097);
   c = a - floor((b * 146097) / 4);

   d = floor((4 * c + 3) / 1461);
   e = c - floor((1461 * d) / 4);
   m = floor((5 * e + 2) / 153);

   day   = e - floor((153 * m + 2) / 5) + 1;
   month = m + 3 - 12 * floor(m / 10);
   year  = b * 100 + d - 4800 + floor(m / 10);

return
