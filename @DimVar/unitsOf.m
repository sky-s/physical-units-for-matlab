function [uv,us] = unitsOf(v)
% uv = UNITSOF(V) returns a value = 1 scalar dimensioned variable with same
% units as V.
% 
%   unitsOf(V) == V/u2num(V) or V/double(V).
% 
%   [uv,us] = UNITSOF(V) returns the the same uv as with s single output
%   argument but in addition returns a string, us, describing the units of V.
% 
%   Example: Use unitsOf to use MATLAB functions that are undefined for
%   dimensioned variables. (Note: MESHGRID is actually defined for dimensioned
%   variables - this is just an example.)
%       V = u.kph*(0:5:80); P = u.kW*(15:10:65);
%       tempV = u2num(V); tempP = u2num(P);
%       [tempV,tempP] = meshgrid(tempV,tempP);
%       V = tempV*unitsOf(V); P = tempP*unitsOf(P); 
% 
%   See also U, U2NUM. DISPLAYINGVALUE.

uv = DimVar(v.exponents,1);

if nargout > 1
    [~,us] = displayparser(v);
end