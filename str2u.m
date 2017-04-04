function out = str2u(inStr)
% STR2U  Convert a string representing physical units to DimVar by evaluating
% the input after prepending 'u.' to valid substrings.
% 
%   If the input is a cell array of strings, STR2U returns a cell array of
%   DimVars of the same size.
% 
%   Compound units are allowed with operators * and - for multiplication and /
%   for division. The characters ² and ³ are also interpreted as ^2 and ^3,
%   respectively. Other operators will be passed to the eval function.
% 
%   Grouping with parentheses for clarity is advisable. Note that
%   str2u('km/h-s') does not return the same result as str2u('km/(h-s)').
% 
%   Example: str2u('kg-m²/s^3') returns a DimVar with units of watts (same as
%   calling u.W).
% 
%   See also u, eval.

%   Sky Sartorius 
%   www.mathworks.com/matlabcentral/fileexchange/authors/101715

if iscellstr(inStr)
    out = cellfun(@str2u,inStr,'UniformOutput',0);
    return
end

if isempty(inStr)
    out = 1;
    return
end

validateattributes(inStr,{'char'},{'row'},'str2u');

out = eval(regexprep(inStr,{'([A-Za-z]+\w*)' '-(?=[A-Za-z]+)'  '²'  '³'},...
                           {'u.$0'           '*'               '^2' '^3'}));
% TODO: allow '-' to have different behavior from '*' such that str2u('km/h-s')
% = str2u('km/(h-s)') or str2u('km/h-s-K') also works as str2u('km/(h-s-K)')