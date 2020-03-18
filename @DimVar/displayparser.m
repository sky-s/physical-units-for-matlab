function [dispVal,unitStr,numString,denString,labelStr] = displayparser(dispVar)
% Parse a DimVar into useful values and strings for display, etc. 
% [dispVal,unitStr,numString,denString,labelStr] = displayparser(v)
% 
%   See also u, DimVar.disp, DimVar.display, DimVar.plot, xlabel.

dispVal = dispVar.value;

numString = '';
denString = '';

%% Preferred units.
% First use customDisplay if there is one.
if ~isempty(dispVar.customDisplay)
    str = dispVar.customDisplay;
    if isprop(u,str)
        test = dispVar/u.(str);
    else
        test = dispVar/str2u(str);
        % avoid the overhead of str2u in most cases.
    end
    if ~isa(test, 'DimVar')
        % Units match.
        numString = str;
        dispVal = test;
        buildAppendStr();
        return
    end
end

% Determine if it matches a preferred unit. Preferred units can be list or
% 2-column cell array.
if isempty(u.dispUnits)
    % Do nothing.
elseif iscellstr(u.dispUnits)
    for i = 1:length(u.dispUnits)
        str = u.dispUnits{i};
        if isprop(u,str)
            test = dispVar/u.(str);
        else
            test = dispVar/str2u(str);
            % avoid the overhead of str2u in most cases.
        end
        if ~isa(test, 'DimVar')
            % Units match.
            numString = str;
            dispVal = test;
            buildAppendStr();
            return
        end
    end
elseif iscell(u.dispUnits)
    prefStrings = u.dispUnits(:,1);
    prefUnits = u.dispUnits(:,2);
    for i = 1:numel(prefStrings)
        test = dispVar/prefUnits{i};
        if ~isa(test, 'DimVar')
            % Units match.
            numString = prefStrings{i};
            dispVal = test;
            buildAppendStr();
            return
        end
    end
else
    error('dispUnits must be cellstr or 2-column cell array.')
end

if nargout == 1
    return
end
%% Built from base units.
names = u.baseNames;

for nd = 1:numel(names)
    currentExp = dispVar.exponents(nd);
    [n,d] = rat(currentExp);
    if currentExp > 0 % Numerator
        if d == 1
            switch currentExp
                case 1
                    numString = sprintf('%s[%s]',numString,names{nd});
                case 2
                    numString = sprintf('%s[%s²]',numString,names{nd});
                case 3
                    numString = sprintf('%s[%s³]',numString,names{nd});
                otherwise
                    numString = sprintf('%s[%s^%g]',...
                        numString,names{nd},currentExp);
            end
        else
            numString = sprintf('%s[%s^(%g/%g)]',...
                numString,names{nd},n,d);
        end
    elseif currentExp < 0 %Denominator
        if d == 1 
            switch currentExp
                case -1
                    denString = sprintf('%s[%s]',denString,names{nd});
                case -2
                    denString = sprintf('%s[%s²]',denString,names{nd});
                case -3
                    denString = sprintf('%s[%s³]',denString,names{nd});
                otherwise
                    denString = sprintf('%s[%s^%g]',...
                        denString,names{nd},-currentExp);
            end
        else
            denString = sprintf('%s[%s^(%g/%g)]',...
                denString,names{nd},-n,d);
        end
    end
end

% Trim brakets for lonely terms.
if 1 == nnz(sign(dispVar.exponents) == 1)
    numString = numString(2:end-1);
end
if 1 == nnz(sign(dispVar.exponents) == -1)
    denString = denString(2:end-1);
end
if isempty(numString)
    numString = '1';
end

buildAppendStr();

%%
    function buildAppendStr()
        if isempty(denString)
            unitStr = numString;
        else
            unitStr = sprintf('%s/%s', numString, denString);
        end
        labelStr = regexprep(unitStr,{'(' ')'},{'{' '}'});
    end
end
