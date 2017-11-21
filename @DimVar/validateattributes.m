function validateattributes(a,classes,varargin)
% Categorize DimVar as numeric for purposes of validating attributes.

if any(strcmp('numeric',classes)) || any(strcmp('DimVar',classes))
    % DimVar is allowed. Evaluate against value.
    
    classes = strrep(classes,'DimVar','double');
    % Does not need a 'unique' call since validateattributes is okay
    % with doubling up the classes input (e.g. {'double' 'double'}.
    
    validateattributes(a.value,classes,varargin{:})
else
    builtin('validateattributes',a,classes,varargin{:})
    % Use built-in function to throw appropriate error for
    % non-DimVar input.
end
end