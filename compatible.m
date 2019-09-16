function compatible(varargin)
% compatible(v1, v2, ...) throws an error unless all inputs have the same units.
%
%   If throwing an error is not desired, use iscompatible.
%
%   See also u, iscompatible.

for i = 1:nargin
    if ~isnumeric(varargin{i})
        ME = MException('DimVar:incompatibleUnits',...
            ['Incompatible units. Cannot perform operation on '...
            'variables with different units.']);
        throwAsCaller(ME);
    end
end