function [h, xLabelString, yLabelString, zLabelString] = patch(varargin)
% DimVar.patch  Patch method for DimVar inputs.
% 
%   See also patch, DimVar.plot, DimVar.u2num.

% 2 special cases: ax as first argument; no x or y provided

args = varargin;

if ishandle(varargin{1})
    % First argument is axes.
    args = args(2:end);
    axisInput = true;
else
    axisInput = false;
end

nArgs = numel(args);

labelStrings = {'' '' ''};
for i = 1:nArgs
    if isa(args{i},'DimVar')
        [~, labelStrings{i}, jnk] = num2str(args{i});
        args{i} = jnk.value;
    end
end

s = regexprep(labelStrings,{'(' ')'},{'{' '}'});
% if nArgs < 3
%     % No X or Y provided.
%     xLabelString = s{2};
%     yLabelString = s{3};
%     zLabelString = s{1};
% else
xLabelString = s{1};
yLabelString = s{2};
zLabelString = s{3};
% end

try
    if axisInput
        h_ = patch(varargin{1},args{:});
    else
        h_ = patch(args{:});
    end
    
    a = gca;
    if ~isempty(xLabelString)
        a.XAxis.TickLabelFormat = ['%g ' xLabelString]; % R2015b
        %     xlabel(xLabelUnitString) % Prior versions.
    end
    if ~isempty(yLabelString)
        a.YAxis.TickLabelFormat = ['%g ' yLabelString]; % R2015b
        %     ylabel(yLabelUnitString) % Prior versions.
    end
    if ~isempty(zLabelString)
        a.ZAxis.TickLabelFormat = ['%g ' zLabelString]; % R2015b
        %     zlabel(zLabelUnitString) % Prior versions.
    end
    
catch  ME
    warning('%s\nSee <a href="%s">%s</a>.',...
        'Not all functionality is supported for inputs of type DimVar.',...
        'matlab:help DimVar/u2num','u2num');
    rethrow(ME)
end

if nargout
    h = h_;
end

end