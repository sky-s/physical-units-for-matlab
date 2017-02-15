function [s, appendString, v] = num2str(v, varargin)

[numeratorString, denominatorString, v] = display(v);
s = num2str(v.value, varargin{:});
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
    
    [~,n] = size(v.value);
    if n > 1
        % Indicate row vector using brackets.
        s = ['[' s ']'];
    end
    s = [s ' ' appendString];
end

% Created 2014-08-25
% More outputs 2014-11-05