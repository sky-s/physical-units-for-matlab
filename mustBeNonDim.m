function mustBeNonDim(x)
% mustBeNonDim  Validate input is numeric and non-dimensional.
%
%   mustBeNonDim(X) throws an error if X is not a numeric array or if X is a
%   dimensioned variable (a `DimVar` produced by the Physical Units Toolbox).
%   This function is intended for use in input validation and assert-style
%   checks inside functions that require plain numeric data.
% 
%   See also mustBeNumeric, isnumeric, u.

if ~isnumeric(x) && ~islogical(x)
    throwAsCaller(matlab.internal.validation.util.createValidatorException('MATLAB:validators:mustBeNumericOrLogical'));
end
if isa(x,'DimVar')
    ME = MException('MATLAB:validators:mustBeNonDim', ...
        'Value must be a non-dimensional numeric type without units.');
    ME.throwAsCaller();
end
