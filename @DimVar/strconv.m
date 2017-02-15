function [cTo,cInverse] = strconv(sFrom,sTo)
% [cTo,cInverse] = STRCONV(sFrom,sTo) finds conversion factor, cTo, that
% converts from the unit indicated by the string sFrom to the unit indicated by
% the string sTo. cInverse is the inverse.
% 
%   Examples:
%       DimVar.STRCONV('MPa/min','(lbf/cm^2)/hr')
%       DimVar.STRCONV('deg/min','1/hr')
% 
%   See also u, str2u.

sTo = str2u(sTo);
sFrom = str2u(sFrom);

cTo = sFrom/sTo;

if isa(cTo, 'DimVar')
    error('Incompatible units.')
end

cInverse = 1/cTo;

end

% 2017-01-21 replaced subfunction with new str2u, enabling simplification
% (removing subfunction, removing scrubbing inputs for empty strings, removing
% eval call, ...).