function v1 = times(v1,v2)

if ~isa(v2,'DimVar') % v1 is the only DimVar.
    v1.value = v1.value .* v2;
elseif ~isa(v1,'DimVar') % v2 is the only DimVar.
    v2.value = v1 .* v2.value;
    v1 = v2;
else % BOTH v1 and v2 are DimVars.
    v1.value = v1.value .* v2.value;
    v1.exponents = v1.exponents + v2.exponents;
    v1 = clearcanceledunits(v1);
end