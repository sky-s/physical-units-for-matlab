function v = round(v,varargin)
warning('DimVar:round','Using round with DimVars may yield unexpected results.')

v.value = round(v.value,varargin{:});