function v = plottingvalue(v)
% PLOTTINGVALUE(V)  Returns value of DimVar used in plotting (which varies with
% selected base unit systems and preferred display units).
% 
%   See also scd, displayingvalue, makehgtform, DimVar.u2num, DimVar.plot.

if isa(v,'DimVar')
    v = displayparser(scd(v));
end