function out = str2u(inStr)
% STR2U  Convert a string representing physical units to DimVar by evaluating
% the input after prepending 'u.' to valid substrings.
% 
%   If the input is a cell or string array, STR2U returns a cell array of
%   DimVars of the same size.
% 
%   Compound units are allowed with operators * and - for multiplication and /
%   for division. The characters ² and ³ are also interpreted as ^2 and ^3,
%   respectively. Other operators will be passed to the eval function.
% 
%   Grouping with parentheses for clarity is advisable. Note that
%   str2u('km/h-s') does not return the same result as str2u('km/(h-s)').
% 
%   Examples: 
%     str2u('kg-m²/s^3') returns a DimVar with units of watts (same as calling
%     u.W).
% 
%     str2u('-5km/s') or str2u('-5 km / s') is the same as calling -5*u.km/u.s.
% 
%     str2u returns a cell array for string or string array inputs. 
% 
%   See also u, eval.

%   Sky Sartorius 
%   www.mathworks.com/matlabcentral/fileexchange/authors/101715

% This first try is a shortcut as well as covers some plain number inputs.
out = str2double(inStr);
if ~isnan(out)
    return 
end

if isstring(inStr) && ~isscalar(inStr)
    out = arrayfun(@str2u,inStr,'UniformOutput',0);
    return
end
if iscell(inStr)
    out = cellfun(@str2u,inStr,'UniformOutput',0);
    return
end

if isempty(inStr)
    out = [];
    return
end

validateattributes(inStr,{'char' 'string'},{'row'},'str2u');
inStr = strtrim(inStr);

% Interpret everything prior to the first alphabetic character (incl. case of
% leading - or .) as the value.
exp = {'^[-+.0-9]+' ')('  ']['  '([A-Za-z]+\w*)' '-(?=[A-Za-z]+)'  '²'  '³' };
rep = {'$0*'        ')*(' ']*[' 'u.$0'           '*'               '^2' '^3'};
out = eval(regexprep(inStr,exp,rep));

if isa(out,'DimVar')
    unitStr = regexp(inStr,'^[-+.0-9]+','split');
    out = scd(out,strtrim(unitStr{end}));
end