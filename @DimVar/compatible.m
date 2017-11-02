function compatible(varargin)
% compatible(v1, v2, ...) returns TRUE if all inputs are DimVar with the same
% units and throws an error otherwise.
% 
%   If throwing an error is not desired, use incompatible.
% 
%   See also u, incompatible.

if nargin == 2
    % Simplified and faster version for only two inputs.
    if ~isa(varargin{1},'DimVar') || ~isa(varargin{2},'DimVar') || ...
            ~isequal(varargin{1}.exponents,varargin{2}.exponents)
        ME = MException('DimVar:incompatibleUnits',...
            ['Incompatible units. Cannot perform operation on '...
            'variables with different units.']);
        throwAsCaller(ME);
    end
    return
end

try
    c = cellfun(@(x) round(x.exponents,5),varargin,'UniformOutput',false);
catch ME
    if strcmp(ME.identifier,'MATLAB:structRefFromNonStruct')
        throwAsCaller(MException('DimVar:incompatibleUnits',...
            ['Incompatible units. Cannot perform operation on '...
            'variables with different units.']));
    else
        rethrow(ME)
    end
end

if nargin == 1 || isequal(c{:})
    % Single input is always compatible with itself.
else
    throwAsCaller(MException('DimVar:incompatibleUnits',...
        ['Incompatible units. Cannot perform operation on '...
        'variables with different units.']));
end