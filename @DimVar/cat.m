function v = cat(dim,v,varargin)

for i = 1:length(varargin)
    vi = varargin{i};
    if compatible(v,vi)
        v.value = cat(dim,v.value,vi.value);
    end
end