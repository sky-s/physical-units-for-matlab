function Vq = interp3(varargin)
% See also interp3.

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


if nArgs <= 4
    %  interp3(V,...)
    outUnits = unitsOf(varargin{1});
elseif nArgs == 7
    %  interp3(X,Y,Z,V,Xq,Yq,Zq,...)
    outUnits = unitsOf(varargin{4});
    for i = 1:3
        compatible(varargin{[i,i+4]});
    end
end

if nArgs == nargin - 2
    % interp3(...,METHOD,EXTRAPVAL)
    compatible(outUnits,varargin{end});
end

for i = 1:nargin
    if isnumeric(varargin{i})
        varargin{i} = double(varargin{i});
    end
end

Vq = outUnits*interp3(varargin{:});