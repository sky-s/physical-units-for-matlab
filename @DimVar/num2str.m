function [s, appendString, v, val] = num2str(v, varargin)

[numeratorString, denominatorString, v] = display(v);
val = v.value;
s = num2str(val, varargin{:});
[m, ~] = size(s);


if ~isempty(denominatorString)
    denominatorString = ['/' denominatorString];
end
appendString = sprintf('%s%s', numeratorString, denominatorString);

if m > 1
    % Multiline.
    s = char(s,appendString);
else
    % Single line.
    if isempty(s)
        s = '[]';
    end
    
    [~,n] = size(val);
    if n > 1
        % Indicate row vector using brackets.
        s = ['[' s ']'];
    end
    s = [s ' ' appendString];
end