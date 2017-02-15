function [v1,I] = min(v1,varargin)

if nargin == 1
    [v1.value, I] = min(v1.value);
elseif nargin == 2 
    v2 = varargin{1};
    if compatible(v1,v2)
        v1.value = min(v1.value, v2.value);
    end
else
    if isempty(varargin{1})
        % As of 2014-05-15, there is no isempty method for DV.
        [v1.value, I] = min(v1.value,varargin{:});
    else
        error(['MIN with two matrices to compare and a working'...
            'dimension is not supported.']);
    end
end

% 2014-05-15/Sartorius: Expanded to full functionality of function.