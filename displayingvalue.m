function v = displayingvalue(v)
% DISPLAYINGVALUE(V)  Returns value of DimVar used in displaying and plotting
% (which varies with selected base unit systems and preferred display units,
% either from displayUnits or per variable).
% 
%   See also scd, makehgtform, DimVar.display, DimVar.u2num, DimVar.plot.

if isa(v,'DimVar')
    v = displayparser(v);
end