function result = ne(v1,v2)

if compatible(v1,v2)
    result = v1.value ~= v2.value;
end

% 2014-05-16/Sartorius: new, simpler.