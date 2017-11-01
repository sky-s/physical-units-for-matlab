function d = u2duration(v)
% u2duration  Convert DimVar time to duration data type.
%   D = u2duration(V) converts time value contained in DimVar V to a
%   duration.
% 
%   u2duration inspects the name of the time dimension and recognizes most
%   common string expressions of seconds, minutes, hours, or days. If you
%   have redefined the time dimension name to something that could be
%   ambiguous, use caution.
% 
%   Note that units uses the Julian year (which defines the light-year),
%   rather than the Gregorian year that duration uses.
% 
%   See also duration, seconds, units, u2num.

if isa(v/u.s,'DimVar')
    error('A pure time DimVar (with exponent of one) is required.')
else
    switch lower(v.names{3})
        case {'s' 'ss' 'sec' 'secs' 'second' 'seconds'}
            d = seconds(v.value);
        case {'m' 'mm' 'min' 'mins' 'minute' 'minutes'}
            d = minutes(v.value);
        case {'h' 'hh' 'hr' 'hrs' 'hour' 'hours'}
            d = hours(v.value);
        case {'day' 'days'}
            d = days(v.value);
        otherwise
            error('Unknown time unit name.')
    end        
end