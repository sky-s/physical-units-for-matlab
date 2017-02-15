function v = diag(v,varargin)

v.value = diag(v.value,varargin{:});

% 2014-05-14/Sartorius: new code to allow more inputs, simplified syntax.