function compatible = compatible(v1, v2, throwErrorFlag)
% Checks if two inputs are both DimVars with the same units (within tolerance).
% 
%   compatible(v1, v2) or compatible(v1, v2, 'throwError') returns TRUE if both
%   inputs are DimVar with the same units and throws an error otherwise.
% 
%   TF = compatible(v1, v2, 'noError') will return FALSE instead of throwing an
%   error for incompatible units.
% 
%   See also u.

if isa(v1,'DimVar') && isa(v2,'DimVar') && ...
    all(abs(v1.exponents - v2.exponents) <= v1.exponentsZeroTolerance)
    
    compatible = true;
else
    compatible = false;
    
    % Checking error flag here to reduce the overhead for the nominal case,
    % since this function is called a lot by many DimVar methods.
    if nargin < 3
        throwErrorFlag = 'throwError';
    end
    throwErrorFlag = validatestring(throwErrorFlag,{'noError','throwError'});
    
    if strcmp(throwErrorFlag,'throwError')
        throwAsCaller(MException('DimVar:incompatibleUnits',...
            ['Incompatible units. Cannot perform operation on '...
            'variables with different units.']));
    end

end
