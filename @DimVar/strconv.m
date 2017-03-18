function [cTo,cInverse] = strconv(sFrom,sTo,varargin)
% strconv  Find conversion factor between units defined by strings.
% 
%   [c, cInv] = STRCONV(sFrom, sTo) finds the conversion factor that converts
%   from the unit indicated by the string sFrom to the unit indicated by the
%   string sTo such that sFrom = c*sTo and sTo = cInv*sFrom. If the units are
%   not compatible, an error is thrown.
%
%   [c, cInv] = unitConversionFactor(sFrom, sTo, 'Force', true) allows for
%   incompatible units, in which case the function may return DimVars instead of
%   throwing an error.
% 
%   Examples:
%       DimVar.STRCONV('MPa/min','(lbf/cm^2)/hr')
%       DimVar.STRCONV('deg/min','1/hr')
% 
%   See also u, str2u, unitConversionFactor.

if ischar(sFrom)
    sFrom = str2u(sFrom);
end
if ischar(sTo)
    sTo = str2u(sTo);
end

[cTo,cInverse] = unitConversionFactor(sFrom,sTo,varargin{:});