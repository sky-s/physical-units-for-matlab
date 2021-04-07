function s = multiunitstr(in, sets, options)
% multiunitstr  Return a string representing a DimVar that is a good guess at
% 'natural', readable, understandable units in both SI and imperial units.
% 
%   The returned string is automatically copied to the clipboard.
% 
%   Supports the following (expand when desired):
%       Speeds of magnitude kph.
%       Lengths of magnitude km, m, or cm.
%       Masses of magnitude kg.
% 
%   See also u, scd, displayUnits.

% This method will not scale and so should not be included in the toolbox.
% However, there are a couple of potential approaches to fulfilling the
% function. 
% 
% One would be to allow redundant entries in the displayUnits list such that the
% tool picks the units closest to but greater than 1 (or half?, or log10 closes
% to 0) or so of those (assuming no customdisplay is set). An example would be
% cm, m, and km on the list, where it would choose cm for below 1 m, m for below
% 1 km, etc. This wouldn't allow for elegant handling of cases where multiple
% unit systems want to be carried, though (e.g. km and mi), but could be a cool
% feature to add, which would require a modification of loops in displayparser
% to go through everything and pick the best of multiple matches. Displayparser
% could also get much simpler if dispunits are always enforced to be the
% two-column format, in which case the loops would no longer be necessary (since
% figure out if things are pre-defined vs needing str2u would only happen once).

% Another way would be to have a separate list for project-specific multi-units.
% This is what would handle multi-units, potentially choosing which set (e.g.,
% ft|m versus mi|km) based solely on if the scd contains a member.

% Copyright Sky Sartorius. All rights reserved.
% Contact: www.mathworks.com/matlabcentral/fileexchange/authors/101715 

options.parens = ''; %'()'; % Must be char vector of length 2 (or empty).
options.separator = ' | '; %', ';
options.format = '%g';

% This would be project-specific. Primary unit goes first (maybe: customDisplay
% supersedes as primary?)
sets = {
    ["kph" "ktas" "mph"]
    ["km" "nmi" "mi"]
    ["cm" "in"]
    ["m" "ft"]
    ["kg" "lb"]
    };

base = in.customDisplay;
if isempty(base)
    error('multiunitstr requires a defined custom display unit to anchor.')
    % Consider pulling in the option of using a matching dispUnit. 
end

for i = 1:numel(sets)
    set = sets{i};
    if ismember(base,set)
        s = sprintf(options.format + " " + set(1),in/u.(set(1)));
        n = numel(set);
        if n > 1 && ~isempty(options.parens)
            s = s + " " + options.parens(1) + ...
                sprintf(options.format + " " + set(2),in/u.(set(2)));
            startInd = 3;
        else
            startInd = 2;
        end
        for j = startInd:numel(set)
            s = sprintf("%s%s" + options.format + " " + set(j),...
                s,options.separator,in/u.(set(j)));
        end
        if ~isempty(options.parens)
            s = s + options.parens(end);
        end
        
        break
    end
end


% 
% if iscompatible(in,u.kts)
%     s = sprintf('%.0f kph | %.0f ktas | %.0f mph',...
%                  in/u.kph,  in/u.KTAS,  in/u.mph);
% elseif iscompatible(in,u.ft)
%     if in > 10*u.km
%         s = sprintf('%.0f km | %.0f nmi | %.0f mi',...
%                      in/u.km,  in/u.nmi,  in/u.mi);
%     elseif in < 10*u.m
%         s = sprintf('%.0f cm | %.0f in',...
%                      in/u.cm,  in/u.in);
%     else
%         s = sprintf('%.0f m | %.0f ft',...
%                      in/u.m,  in/u.ft);
%     end
% elseif iscompatible(in,u.kg)
%     s = sprintf('%.0f kg | %.0f lb',...
%                  in/u.kg,  in/u.lb);
%     
% end

copy(char(s))
