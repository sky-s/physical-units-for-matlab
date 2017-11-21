function a = colon(a,b,c)

if nargin == 2
    error('DimVar vector creation with colon is only defined for three inputs.')
end

compatible(a,b,c)

a.value = a.value:b.value:c.value;