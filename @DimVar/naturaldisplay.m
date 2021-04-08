function v = naturaldisplay(v,lists)
% v = naturaldisplay(v) sets custom display units of v such that its magnitude
% close to but not below 1. For arrays, the median value is used to pick a
% magnitude.
% 
%   v = naturaldisplay(v,lists), instead of using the default lists of "normal"
%   units, searches for a match in the provided lists. Each element of the lists
%   cell array is a string array of compatible units defined by u in ascending
%   order of magnitude.
% 
%   lists can contain different string arrays representing same unit type. First
%   priority when selecting a list is always for one matching the input's
%   customDisplay property, otherwise it is the first elements of lists whose
%   first element represents a compatible unit with the input.
% 
%   Example: See light travel times in units that make the most sense.
%     v = u.lightSpeed;
%     naturaldisplay(10*u.meter/v)
%     naturaldisplay(1000*u.mile/v)
%     naturaldisplay(1*u.astronomicalUnit/v)
%     naturaldisplay(1*u.parsec/v)
% 
%   Example: Scale down a metric recipe.
%     lists = {["oz" "lb"];["tsp" "Tbls" "cup" "quart" "gal"]};
%     naturaldisplay(30*u.mL,lists)
%     naturaldisplay(250*u.g,lists)
%     naturaldisplay(500*u.mL,lists)
% 
%   See also DimVar.scd, u, iscompatible.

if nargin < 2
    lists = {
        ["um" "mm" "cm" "m" "km" "au" "lightYear"]
        ["mil" "in" "ft" "mi"]
        ["mil" "in" "ft" "nmi"]
        ["ug" "mg" "g" "kg" "tonne" "Mt"]
        ["oz" "lbm" "lb" "ton"]
        ["ms" "s" "min" "hr" "day" "yr"]
        ["sqmil" "sqin" "sqft" "acre" "sqmi"]
        ["sqmm" "sqcm" "sqm" "ha" "sqkm"]
        ["ul" "ml" "dl" "L" "m3"]
        ["floz" "cup" "quart" "gal" "acft"]
        ["uN" "dyne" "mN" "N" "kN" "MN"]
        ["lbf" "kip"]
        ["" "k" "M" "G" "T"] + "Hz"
        ["n" "u" "m" "" "k" "M" "G" "T"] + "J"
        ["u" "m" "" "k" "M" "G"] + "Pa"
        ["psi" "ksi" "Msi"]
        ["u" "m" "" "k" "M" "G" "T"] + "W"
        ["u" "m" "" "k"] + "A"
        ["u" "m" "" "k" "M"] + "V"
        ["n" "u" "m" "" "k" "M" "G"] + "Ohm"
        ["p" "n" "u" "m" ""] + "F"
        ["n" "u" "m" "" "k" "M" "G"] + "H"
        ["bit" ["" "k" "M" "G" "T" "P" "E"] + "B"]
        ["" "k" "M" "G" "T"] + "bps"
        };
end

%% Find a matching list
ln = 0;


for i = numel(lists):-1:1
    if any(lists{i} == v.customDisplay)
        % Stop search if finding an exact match.
        ln = i;
        break
    elseif iscompatible(v,u.(lists{i}(1)))
        % Then search for any compatible match, prioritizing higher-up entries
        % (hence the reason for searching backwards through the list).
        ln = i;
    end
end
if ~ln
    return
end

s = lists{ln};

%% Choose a unit from the list
m = median(v(:));

% Loop from largest to smallest, stopping when the input divided by that unit is
% above 1.
for i = numel(s):-1:1
    k = m/u.(s{i});
    if k>1
        break
    end
end

v = scd(v,s{i});
