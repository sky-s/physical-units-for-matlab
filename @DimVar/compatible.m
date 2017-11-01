function compatible = compatible(varargin)
% Checks if two or more inputs are all DimVars with the same units (within
% tolerance).
% 
%   compatible(v1, v2, ...) returns TRUE if all inputs are DimVar with the same
%   units and throws an error otherwise.
% 
%   If throwing an error is not desired, use incompatible.
% 
%   See also u, incompatible.

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

if isequal(c{:})
    compatible = true;
else
    throwAsCaller(MException('DimVar:incompatibleUnits',...
        ['Incompatible units. Cannot perform operation on '...
        'variables with different units.']));
end