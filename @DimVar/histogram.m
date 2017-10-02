function [h, xLabelUnitString] = histogram(varargin)
% DimVar.histogram  Histogram method that automatically labels axes with units
% of DimVar inputs.
% 
%   [h, xLabelUnitString] = histogram(...) in addition to the handle, returns
%   the string used by the function for labeling the axes.
% 
%   See also histogram.

%   Copyright 2017 Sky Sartorius
%   www.mathworks.com/matlabcentral/fileexchange/authors/101715 

x = varargin{1};

[~, xs] = num2str(x);
xLabelUnitString = regexprep(xs,{'(' ')'},{'{' '}'});

% Go through arguments to check compatibility with arguments such as bin width,
% etc.
for i = 1:numel(varargin)
    if isa(varargin{i},'DimVar')
        compatible(varargin{i},x); % Check for compatibility.
        
        [~, ~, varargin{i}] = display(varargin{i});  
        % Adjust value for display units.
        
        varargin{i} = varargin{i}.value;
    end
end

try
    h_ = histogram(varargin{:});
    a = gca;
    a.XAxis.TickLabelFormat = ['%g ' xLabelUnitString]; % R2015b
    % xlabel(xLabelUnitString) % Prior versions.
catch ME
    warning('%s\nSee <a href="%s">%s</a>.',...
        'Not all functionality is supported for inputs of type DimVar.',...
        'matlab:help DimVar/u2num','u2num');
    rethrow(ME)
end

if nargout
    h = h_;
end