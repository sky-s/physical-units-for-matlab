function v1 = minus(v1,v2)

if compatible(v1,v2)
    v1.value = v1.value - v2.value;
end

% 2014-05-15/Sartorius: new, majorly simplified.