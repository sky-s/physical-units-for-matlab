function v = cat(dim,v,varargin)

if ~isa(v,'DimVar')
    
    ME = MException('DimVar:incompatibleUnits',...
        'Incompatible units. All inputs must be DimVar.');
    throwAsCaller(ME);
    
end

vExpos = v.exponents;

for i = 1:numel(varargin)
    vi = varargin{i};
    if ~isa(vi,'DimVar') || ~isequal(vExpos,vi.exponents)
        
        ME = MException('DimVar:incompatibleUnits',...
            ['Incompatible units. Cannot perform operation on '...
            'variables with different units.']);
        throwAsCaller(ME);
        
    end
    v.value = cat(dim,v.value,vi.value);
end

% Functionality of compatible method is integrated for the sake of speed.