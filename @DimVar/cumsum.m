function v = cumsum(v,varargin)

v.value = cumsum(v.value,varargin{:});

% 2014-05-14/Sartorius: updated to allow more inputs, simplified syntax.