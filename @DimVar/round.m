function v = round(v,varargin)

warning('DimVar:round','Using round with DimVars may yield unexpected results.')

dispVal = displayparser(v);
delta = round(dispVal,varargin{:})./dispVal;

v.value = v.value.*delta;

% v.value = round(v.value,varargin{:});