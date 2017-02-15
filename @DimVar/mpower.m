function v = mpower(v,y)
% 

if isa(y,'DimVar')
    error('Exponent may not be a DimVar.');
else
    v.value = v.value^y; % ^ will catch if y isn't scalar.
    v.exponents = y*v.exponents;
    
    v = clearcanceledunits(v);
end

% 2014-05-16/Sartorius: reworked.