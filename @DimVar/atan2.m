function out = atan2(v1,v2)

if compatible(v1,v2)
    out = atan2(v1.value, v2.value);
end

% Added Oct 2015.