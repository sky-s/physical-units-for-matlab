function Vq = interp2(varargin)
% See also interp2.

nArgs = nargin;
if isnumeric(varargin{end}) && isscalar(varargin{end}) ...
        && (ischar(varargin{end-1}) ...
        || (isstring(varargin{end-1}) && isscalar(varargin{end-1})))
    % User supplied a method and extrap val.
    nArgs = nArgs-2;
elseif ischar(varargin{end}) ...
        || (isstring(varargin{end}) && isscalar(varargin{end}))
    nArgs = nArgs-1;
end


if nArgs <= 3
    % interp2(V,...)
    outUnits = unitsOf(varargin{1});
else
    % interp2(X,Y,V,Xq,Yq,...)
    outUnits = unitsOf(varargin{3});
    compatible(varargin{[1,4]});
    compatible(varargin{[2,5]});
end

if nArgs == nargin - 2
    % interp2(...,METHOD,EXTRAPVAL)
    compatible(outUnits,varargin{end});
end

for i = 1:nargin
    if isnumeric(varargin{i})
        varargin{i} = double(varargin{i});
    end
end

Vq = outUnits*interp2(varargin{:});