function v = displayingvalue(v)
% DISPLAYINGVALUE(V)  Returns value of DimVar used in displaying (which varies
% with selected base unit systems and preferred display units, either from
% displayUnits or per variable).
% 
%   See also scd, plottingvalue, makehgtform, DimVar.display, DimVar.u2num.

if isa(v,'DimVar')
    v = displayparser(v);
end