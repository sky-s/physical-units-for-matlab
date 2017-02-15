function vOut = mldivide(v1,v2)


if ~isa(v2,'DimVar') % v1 is the DimVar.
    vOut = v1;
    vOut.value = v1.value \ v2;
    vOut.exponents = - v1.exponents;
    
elseif ~isa(v1,'DimVar') % v2 is the DimVar.
    vOut = v2;
    vOut.value = v1 \ v2.value;


else % BOTH v1 and v2 are DimVars.
    vOut = v1;
    vOut.value = v1.value \ v2.value;
    vOut.exponents = v2.exponents - v1.exponents;
    
    vOut = clearcanceledunits(vOut);
end

% 2014-05-16/Sartorius: Modestly simplified.