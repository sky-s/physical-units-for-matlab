function [varargout] = unitsOf(in)
% Called when unitsOf is used on a variable without physical units.
% 
%   See also DimVar.unitsOf.

if ~isnumeric(in)
    error('Non-numeric input.')
end

switch nargout
    case 0
        disp('-no units-');
    case 1
        varargout = {1};
    case 2
        varargout = {1,''};
end
