function [s] = num2str(v, varargin)

ST = dbstack(1,'-completenames');
if isequal(varargin,{'%d    '}) && ...
        contains(ST(1).file,['@tabular' filesep 'disp.m'])
    % Compensation for behavior of @tabular\disp.m trying to format with %d:
    varargin = {getFloatFormats()}; % Assume double value.
end


[dispVal,unitStr] = displayparser(v);
s = num2str(dispVal, varargin{:});

if isempty(s)
    s = '[]';
end

s = strcat(s,[' ' unitStr]);

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