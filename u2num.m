function [out] = u2num( varargin )
% This function is called when u2num is used on a non-DimVar.
% See also DimVar.u2num.
out = varargin{1};
if ~isnumeric(out)
    error('Non-numeric input.')
end
end
