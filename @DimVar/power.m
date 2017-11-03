function v = power(v,y)

if isa(y,'DimVar')
    error('For Z = X.^Y, Y may not be a DimVar.')
end
v.value = v.value.^y;
v.exponents = y*v.exponents;
    
v = clearcanceledunits(v);