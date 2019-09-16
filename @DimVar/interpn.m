function Vq = interpn(varargin)
% See also interpn.

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


if nArgs <= 2
   %  interpn(V,...)   
   outUnits = unitsOf(varargin{1});
else
   if ~isvector(varargin{1}) && nArgs == (ndims(varargin{1})+1)
       % interpn(V,X1q,X2q,X3q,...)
       outUnits = unitsOf(varargin{1});
   elseif rem(nArgs,2) == 1
       % interpn(X1,X2,X3, ... V,X1q,X2q,X3q, ...)
       vPos = (nArgs + 1)/2;
       outUnits = unitsOf(varargin{vPos});
       for i = 1:vPos-1
           compatible(varargin{[i,i+vPos]});
       end
   end
end

if nArgs == nargin - 2
    % interpn(...,METHOD,EXTRAPVAL)
    compatible(outUnits,varargin{end});
end

for i = 1:nargin
    if isnumeric(varargin{i})
        varargin{i} = double(varargin{i});
    end
end

Vq = outUnits*interpn(varargin{:});