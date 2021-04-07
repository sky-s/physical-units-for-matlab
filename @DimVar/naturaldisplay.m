function v = naturaldisplay(v,sets)
% naturaldisplay  Sets custom display of v to the units that bring its magnitude
% close but not below 1.
% 
%   See also DimVar.scd.

if nargin < 2
    sets = {
        ["um" "mm" "cm" "m" "km" "au" "lightYear"]
        ["mil" "in" "ft" "mi"]
        ["mil" "in" "ft" "nmi"]
        ["ug" "mg" "g" "kg" "tonne"]
        ["oz" "lbm" "lb" "ton"]
        ["ms" "s" "min" "hr" "day" "yr"]
        ["sqmil" "sqin" "sqft" "acre" "sqmi"]
        ["sqmm" "sqcm" "sqm" "ha" "sqkm"]
        ["ul" "ml" "dl" "L" "m3"]
        ...["oz" "gal"]
        ...["uA
        ...["uK'
        ...["bit" "byte" "
        };
end

%%
ln = 0;
for i = numel(sets):-1:1
    if any(sets{i} == v.customDisplay)
        ln = i;
        break
    elseif iscompatible(v,u.(sets{i}(1)))
        ln = i;
    end
end
if ~ln
    return
end

s = sets{ln};

m = median(v(:));
for i = numel(s):-1:1
    k = m/u.(s{i});
    if k>1
        break
    end
end

v = scd(v,s{i});
