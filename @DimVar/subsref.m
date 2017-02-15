function v = subsref(v,varargin)

v.value = subsref(v.value, varargin{:});

% 2014-05-16/Sartorius: new, simpler, way faster (no more eval).