function v = clearcanceledunits(v)
% If all DimVar unit exponents are zero, return normal (double) variable.
% Exponent tolerance is to fifth decimal.

if ~nnz(round(v.exponents,5))
    v = v.value;
end