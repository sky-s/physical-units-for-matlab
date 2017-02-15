function v1 = plus(v1,v2)

if compatible(v1,v2)
    v1.value = v1.value + v2.value;
end

% 2014-05-16/Sartorius: new, simpler, more elegant.