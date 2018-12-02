function v = duration2u(d)
% u2duration Convert duration DimVar.
%   V = duration2u(D) simply executes v = u.s*seconds(D) after checking to
%   ensure that D is in fact a duration type.
% 
%   See also duration, seconds, units, u.

if isduration(d)
    v = u.s*seconds(d);
    switch d.Format
        % Approximately mimic display formats used by duration.
        case 'y'
            v = scd(v,'yr');
        case 'd'
            v = scd(v,'day');
        case 'h'
            v = scd(v,'hr');
        case 'm'
            v = scd(v,'min');
        case 's'
            v = scd(v,'sec');
        otherwise
            v = scd(v);
    end
else
    error('Input must be duration type.')
end