function [varargout] = unitsOf(in)
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
% This function is called when unitsOf is used on a non-DimVar.
