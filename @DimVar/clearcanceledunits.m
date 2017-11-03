function v = clearcanceledunits(v)
% If all DimVar unit exponents are zero, return normal (double) variable.
% Exponent tolerance is to fifth decimal.

if ~any(round(1e5*v.exponents)) % Seems to be faster than round(x,5).
    v = v.value;
end