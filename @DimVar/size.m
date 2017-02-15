function varargout = size(x,varargin)

[varargout{1:nargout}] = size(x.value,varargin{:});

% Updated 2014-05-13 at suggestion of Emanuele Ruffaldi to include full
% functionality of size function.