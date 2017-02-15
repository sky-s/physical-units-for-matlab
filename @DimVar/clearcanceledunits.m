function v = clearcanceledunits(v)
% If all DV unit exponents are within tolerance, return normal, non-DV
% variable.

if all( abs(v.exponents) <= v.exponentsZeroTolerance )
    v = v.value;
end