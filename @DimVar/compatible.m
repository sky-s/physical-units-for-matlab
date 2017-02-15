function compatible = compatible(v1, v2)
% Checks if two inputs are both DVs with the same units (within tolerance).

if isa(v1,'DimVar') && isa(v2,'DimVar') && ...
    all(abs(v1.exponents - v2.exponents) <= v1.exponentsZeroTolerance);
    
    compatible = true;
else
    throwAsCaller(MException('DimVar:incompatibleUnits',...
        ['Incompatible units. Cannot perform operation on '...
        'variables with different units.']));
end