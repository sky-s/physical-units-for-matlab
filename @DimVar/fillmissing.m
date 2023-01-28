function [v,I] = fillmissing(v,varargin)

if strcmpi(varargin{1},'constant')
    c = varargin{2};
    compatible(v,c)
    varargin{2} = c.value;
end
[v.value,I] = fillmissing(v.value,varargin{:});
