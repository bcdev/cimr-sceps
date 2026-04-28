%------------------------------------------------------------------------

function [utc, utcdaydiff] = loctime2utc(lt,lon)
%TOOLS_UTC2LOCTIME Conversion from local time to UTC
%   Detailed explanation goes here
if isdatetime(lt)
    lt = hour(lt) + minute(lt)/60;
end
utc = lt - lon*12/180;
utcdaydiff = zeros(size(utc));

ind = find( utc >= 24);
utc(ind)      = utc(ind) - 24;
utcdaydiff(ind) = 1;


ind = find( utc < 0);
utc(ind)      = utc(ind) + 24;
utcdaydiff(ind) = -1;

return
