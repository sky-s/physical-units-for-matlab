function [h, xLabelString, yLabelString, zLabelString] = surf(varargin)
% DimVar.surf  3D surface plotting method that automatically labels axes with
% units of DimVar inputs.
% 
%   See also DimVar.plot.

%   Copyright 2016 Sky Sartorius
%   www.mathworks.com/matlabcentral/fileexchange/authors/101715 

% 2 special cases: ax as first argument; no x or y provided

args = varargin;

if isscalar(varargin{1}) && isgraphics(varargin{1},'axes')
    % First argument is axes.
    args = args(2:end);
    axisInput = true;
else
    axisInput = false;
end

labelStrings = {'' '' ''};
numericInds = cellfun(@isnumeric,args);
numericArgs = args(numericInds);
nNumericArgs = numel(numericArgs);
for i = 1:nNumericArgs
    if isa(numericArgs{i},'DimVar')
        [numericArgs{i},~,labelStrings{i}] = displayparser(numericArgs{i});
    end
end
args(numericInds) = numericArgs;

labelStrings = regexprep(labelStrings,{'(' ')'},{'{' '}'});
if nNumericArgs < 2
    % Only Z provided.
    xLabelString = labelStrings{2};
    yLabelString = labelStrings{3};
    zLabelString = labelStrings{1};
else
    xLabelString = labelStrings{1};
    yLabelString = labelStrings{2};
    zLabelString = labelStrings{3};
end

try
    if axisInput
        h_ = surf(varargin{1},args{:});
    else
        h_ = surf(args{:});
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
    
catch ME
    warning('%s\nSee <a href="%s">%s</a>.',...
        'Not all functionality is supported for inputs of type DimVar.',...
        'matlab:help DimVar/u2num','u2num');
    rethrow(ME)
end

if nargout
    h = h_;
end

end