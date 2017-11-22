function [varargout] = unitsOf(in)
% Called when unitsOf is used on a variable without physical units.
% 
%   See also DimVar.unitsOf.

if ~isnumeric(in)
    error('Non-numeric input.')
end

varargout = {1,''};