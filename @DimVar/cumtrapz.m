function vOut = cumtrapz(v1,varargin)

if nargin == 1 || (nargin > 1 && ~isa(varargin{1},'DimVar'))
    % v1 is the only DimVar.
    vOut = v1;
    vOut.value = cumtrapz(v1.value,varargin{:});

elseif ~isa(v1,'DimVar')
    % v2 is the only DimVar.
    v2 = varargin{1};
    vOut = v2;
    vOut.value = cumtrapz(v1,v2.value,varargin{2:end});

else 
    % BOTH v1 and v2 are DimVars.
    v2 = varargin{1};
    vOut = v1;
    vOut.value = cumtrapz(v1.value,v2.value,varargin{2:end});
    vOut.exponents = v1.exponents + v2.exponents;
    
    vOut = clearcanceledunits(vOut);
end
