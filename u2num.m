function [out] = u2num( varargin )
out = varargin{1};
if ~isnumeric(out)
    error('Non-numeric input.')
end
end
% This function is called when u2num is used on a non-DimVar.
