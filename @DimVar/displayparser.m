function [dispVal,dispVar,unitStr,numString,denString] = displayparser(dispVar)
% Parse a DimVar into useful values and strings for display, etc. 
% [dispVal,dispVar,unitStr,numString,denString] = displayparser(v)
% 
%   To make a TeX-interpreted label: regexprep(unitStr,{'(' ')'},{'{' '}'});
% 
%   See also DimVar.disp, DimVar.display, DimVar.num2str, DimVar.plot, xlabel.

numString = '';
denString = '';

dispVal = dispVar.value;

%% Preferred units.
% Determine if it matches a preferred unit. Preferred units can be list or
% 2-column cell array.
if isempty(dispVar.dispUnits)
    % Do nothing.
elseif iscellstr(dispVar.dispUnits)
    for i = 1:length(dispVar.dispUnits)
        str = dispVar.dispUnits{i};
        test = dispVar/u.(str);
        if ~isa(test, 'DimVar')
            % Units match.
            numString = str;
            denString = '';
            dispVar.value = test;
            dispVal = test;
            buildAppendStr();
            return
        end
    end
elseif iscell(dispVar.dispUnits)
    prefStrings = dispVar.dispUnits(:,1);
    prefUnits = dispVar.dispUnits(:,2);
    for i = 1:numel(prefStrings)
        test = dispVar/prefUnits{i};
        if ~isa(test, 'DimVar')
            % Units match.
            numString = prefStrings{i};
            denString = '';
            dispVar.value = test;
            dispVal = test;
            buildAppendStr();
            return
        end
    end
else
    error('dispUnits must be cellstr or 2-column cell array.')
end

if nargout <= 2
    return
end
%% Built from base units.
names = dispVar.names;

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
    end
end
