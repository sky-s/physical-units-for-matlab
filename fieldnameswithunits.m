function [names,S,T,C] = fieldnameswithunits(S,exceptionNames,...
    exceptionUnitStr,exceptionUnits)
% fieldnameswithunits  Return field names with (units) appended.
% 
%   names = fieldnameswithunits(S) returns a cell array of character vectors
%   containing the names of the fields in structure S. If a field contained
%   DimVars, the name with be appended with parentheses containing the string
%   representing the unit.
% 
%   [names,S] = fieldnameswithunits(S) returns a NEW version of S that contains
%   only the value of the variables without the corresponding unit.
% 
%   [names,S,T] = fieldnameswithunits(S) returns struct2table(S).
% 
%   [names,S,T,C] = fieldnameswithunits(S) returns [names';table2cell(T)]
% 
%   ... = fieldnameswithunits(S,exceptionNames,exceptionUnitStr) will override
%   the currently set display units for the fields matching any in the cell
%   array of strings exceptionNames and instead will use exceptionUnitStr for
%   the returned names and str2u(exceptionUnitStr) for de-uniting returned S.
% 
%   ... = fieldnameswithunits(S,exceptionNames,exceptionUnitStr,exceptionUnits)
%   uses exceptionUnits instead of str2u(exceptionUnitStr).
% 
%   [names,S,T,C] = fieldnameswithunits(S,exceptionNames,...
%       exceptionUnitStr,exceptionUnits)
% 
%   See also u, displayparser str2u. 


% Copyright 2018 Sky Sartorius
% Contact: www.mathworks.com/matlabcentral/fileexchange/authors/101715

if nargin < 2
    exceptionNames = {};
end
if nargin == 3
    exceptionUnits = str2u(exceptionUnitStr);
%     exceptionUnitStr = cellstr(exceptionUnitStr);
elseif nargin < 3
    exceptionUnits = {};
end
if numel(exceptionNames) ~= numel(exceptionUnits)
    error('Exception inputs must have the same number of elements.')
end
if ~(iscellstr(exceptionNames) || isstring(exceptionNames))
    error('exceptionNames must be cellstr or string array.')
end

names = fieldnames(S);
T = struct2table(S);

uN = names;

for i = 1:numel(names)
    fn = names{i};
    var = T.(fn);
    exceptionInd = strcmp(fn,exceptionNames);
    if ~isa(var,'DimVar')
        uN{i} = '';
    elseif any(exceptionInd)
        uN{i} = exceptionUnitStr{exceptionInd};
        if ~isempty(uN{i})
            names{i} = strcat(fn,' (',uN{i},')');
        end
        T.(fn) = var./exceptionUnits{exceptionInd};
    else
        [T.(fn),~,uN{i}] = displayparser(var);
        names{i} = strcat(fn,' (',uN{i},')');
    end
end
if nargout > 1
    S = table2struct(T);
    if nargout > 2
        T = struct2table(S,'AsArray',true);
        % T.Properties.VariableUnits = uN;
        if nargout > 3
            C = [names';table2cell(T)];
        end
    end
end