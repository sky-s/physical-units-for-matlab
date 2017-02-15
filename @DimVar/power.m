function v = power(v,y)

if ~isa(y,'DimVar') && isscalar(y)
    v.value = v.value.^y;
    v.exponents = y*v.exponents;
    
    v = clearcanceledunits(v);
else
    error('For DV.^b, b must be dimensionless scalar.')
end

% 2014-05-16/Sartorius: reworked.

