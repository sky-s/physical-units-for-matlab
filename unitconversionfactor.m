function [cTo,cInverse] = unitConversionFactor(from,to,varargin)
% unitConversionFactor  Find conversion factor between units.
% 
%   [c, cInv] = unitConversionFactor(FROM, TO) returns a factor c such that FROM
%   is equal to c*TO and TO = cInv*FROM. If the units of FROM and TO are not
%   compatible, an error is thrown.
%
%   [c, cInv] = unitConversionFactor(FROM, TO, 'Force', true) allows for
%   incompatible units, in which case the function may return DimVars instead of
%   throwing an error. This is the same as calling c = FROM./TO.
% 
%   See also u, strconv, str2u, compatible, symunit/unitConversionFactor.

if ischar(from)
    from = str2u(from);
end
if ischar(to)
    to = str2u(to);
end

p = inputParser;
p.addParameter('Force',false);
p.parse(varargin{:});

if p.Results.Force
    % Do nothing if anything was provided for Force (string, double, etc.).
else
    % Throw error if incompatible.
    compatible(from,to);
end


cTo = from./to;

cInverse = 1./cTo;

end