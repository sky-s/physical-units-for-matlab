function v = diff(v,varargin)

v.value = diff(v.value, varargin{:});

% 2014-05-14/Sartorius: simplified; allow full inputs for diff.