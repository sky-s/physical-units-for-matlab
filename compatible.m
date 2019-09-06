function compatible(varargin)
% See also DimVar.compatible.

for i = 1:nargin
    if ~isnumeric(varargin{i})
        ME = MException('DimVar:incompatibleUnits',...
            'Incompatible units.');
        throwAsCaller(ME);        
    end
end