function Vout = interp1(varargin)
% See also interp1.

if nargin == 2
    % interp1(V,Xq)
    outUnits = unitsOf(varargin{1});
else
    outUnits = unitsOf(varargin{2});
end

if nargin >= 3 && isnumeric(varargin{3})
    % interp1(X,V,Xq,...)
    compatible(varargin{[1,3]});
end
if nargin == 5 && isnumeric(varargin{5})
    % interp1(X,V,Xq,METHOD,EXTRAPVAL)
    compatible(varargin{[2,5]});
end

for i = 1:nargin
    if isnumeric(varargin{i})
        varargin{i} = double(varargin{i});
    end
end

Vout = outUnits*interp1(varargin{:});