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
%   str2u('km/h-s') does not return the same result as str2u('km/h*s') because
%   in the former case, the hyphenated h-s is grouped in the denominator.
% 
%   The returned variable will have the unit portion of the input string as its
%   custom display unit.
% 
%   Examples: 
%     str2u('kg-m²/s^3') returns a DimVar with units of watts (u.W).
% 
%     str2u('-5km/s') or str2u('-5 km / s') is the same as calling -5*u.km/u.s.
% 
%     str2u returns a cell array for string array inputs. 
% 
%   See also u, eval.

%   Copyright Sky Sartorius 
%   www.mathworks.com/matlabcentral/fileexchange/authors/101715

% This first try is a shortcut as well as covers some plain number inputs.
out = str2double(inStr);
if ~isnan(out)
    return 
end

%% Parse inputs.
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

%% First separate out the leading number.
[number, unitStr] = regexp(inStr,'^[-+.0-9]+','match','split');
if ~isempty(number)
    number = number{1};
else
    number = '1';
end
unitStr = strtrim(unitStr{end});
if isempty(unitStr)
    out = eval(number);
    return
end

%% Build the more complex expressions.

normalExpo = '(\^-?[.0-9]+)'; % Numeric exponent.
parenExpo = '(\^\(-?[.0-9]+(/[.0-9]+)?\))'; % Exponent with parens.
validUnitStr = '([A-Za-z]+\w*)'; % Valid field names, essentially.

unitWithExponent = sprintf('(%s(%s|%s)?)',validUnitStr,normalExpo,parenExpo);
hypenated = sprintf('%s(-%s)+',unitWithExponent,unitWithExponent);

%% Regexp and eval.

exp = {
    '²'                 % 1 Squared character.
    '³'                 % 2 Cubed character.
    '(^per |^per-|^/)'  % 3 Leading 'per' special case.
    '( per |-per-)'     % 4 Replace per with /
    hypenated           % 5 Group hyphen units with parens.
    ')('                % 6 Multiply back-2-back parens.
    ']['                % 7 Multiply back-2-back brackets.
    validUnitStr        % 8 Precede alphanumeric unit w/ u.
    '-u\.'              % 9 - leading unit is *.
    };
rep = {
    '^2'                % 1
    '^3'                % 2
    '1/'                % 3
    '/'                 % 4
    '($0)'              % 5
    ')*('               % 6
    ']*['               % 7
    'u.$0'              % 8
    '*u.'               % 9
    };                

evalStr = regexprep(unitStr,exp,rep);
out = eval([number '*' evalStr]);

if isa(out,'DimVar')
    out = scd(out,strtrim(unitStr));
end