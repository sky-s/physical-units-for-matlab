function [s, appendString, v, val] = num2str(v, varargin)
[numeratorString, denominatorString, v] = display(v);
val = v.value;

% Compensation for behavior of @tabular\disp.m trying to format with %d:
ST = dbstack(1,'-completenames');
if isequal(varargin,{'%d    '}) && ...
        contains(ST(1).file,['@tabular' filesep 'disp.m'])
    varargin = {getFloatFormats()}; % Assume double value.
end

s = num2str(val, varargin{:});

if ~isempty(denominatorString)
    denominatorString = ['/' denominatorString];
end
appendString = sprintf('%s%s', numeratorString, denominatorString);

if isempty(s)
    s = '[]';
end

s = strcat(s,[' ' appendString]);
end


function dblFmt = getFloatFormats()

switch lower(matlab.internal.display.format)
    case {'short' 'shortg' 'shorteng'}
        dblFmt  = '%.5g    ';
    case {'long' 'longg' 'longeng'}
        dblFmt  = '%.15g    ';
    case 'shorte'
        dblFmt  = '%.4e    ';
    case 'longe'
        dblFmt  = '%.14e    ';
    case 'bank'
        dblFmt  = '%.2f    ';
    otherwise % rat, hex, + fall back to shortg
        dblFmt  = '%.5g    ';
end
end