function compatible(v,varargin)
% compatible(v1, v2, ...) returns TRUE if all inputs are DimVar with the same
% units and throws an error otherwise.
% 
%   If throwing an error is not desired, use iscompatible.
% 
%   See also u, iscompatible.

if ~isa(v,'DimVar')
    
    ME = MException('DimVar:incompatibleUnits',...
        'Incompatible units. All inputs must be DimVar.');
    throwAsCaller(ME);
    
end

vExpos = v.exponents;

for i = 1:numel(varargin)
    
    if ~isa(varargin{i},'DimVar') || ~isequal(vExpos,varargin{i}.exponents)
        
        ME = MException('DimVar:incompatibleUnits',...
            ['Incompatible units. Cannot perform operation on '...
            'variables with different units.']);
        throwAsCaller(ME);
        
    end
end