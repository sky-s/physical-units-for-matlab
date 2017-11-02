function v = cat(dim,varargin)

compatible(varargin{:})
val = cellfun(@(x)x.value,varargin,'UniformOutput',false);

v = varargin{1};
v.value = cat(dim,val{:});
