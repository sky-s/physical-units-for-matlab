function [cTo,cInverse] = unitconversionfactor(from,to,varargin)
% unitconversionfactor  Find conversion factor between units.
% 
%   [c, cInv] = unitconversionfactor(FROM, TO) returns a factor c such that FROM
%   is equal to c*TO and TO = cInv*FROM. If the units of FROM and TO are not
%   compatible, an error is thrown.
%
%   [c, cInv] = unitconversionfactor(FROM, TO, 'Force', true) allows for
%   incompatible units, in which case the function may return DimVars instead of
%   throwing an error. This is the same as calling c = FROM./TO.
% 
%   unitconversionfactor can work with a heterogeneous mix of inputs that can be
%   DimVar, char, numeric, or sym produced by symunit. 
% 
%   Examples:
%       unitconversionfactor(u.radian,u.arcminute)
%       unitconversionfactor('MPa/min','(lbf/cm^2)/hr')
%       unitconversionfactor('deg/min',1/u.hour)
%       unitconversionfactor(u.meter,u.kilogram,'Force',true)
%       unitconversionfactor(symunit('HP'),u.W)
% 
%   See also u, strconv, str2u, compatible, symbolic/unitConversionFactor.

% This implementation using symunit2str my cause issues for some situations
% where the unit names might not line up between u and symunit. This is accepted
% due to the capabilities introduced.
if isa(from,'sym')
    from = symunit2str(from);
end
if isa(to,'sym')
    to = symunit2str(to);
end

if ischar(from)
    from = str2u(from);
end
if ischar(to)
    to = str2u(to);
end

p = inputParser;
p.addParameter('Force',false);
p.parse(varargin{:});

if p.Results.Force || (isfloat(from) && isfloat(to))
    % Do nothing if anything was provided for Force (string, double, etc.), or
    % if both inputs are normal numbers, in which case they are compatible (but
    % are invalid inputs for the compatible function below).
else
    % Throw error if incompatible.
    compatible(from,to); % Not defined for normal inputs.
end


cTo = from./to;

cInverse = 1./cTo;

end