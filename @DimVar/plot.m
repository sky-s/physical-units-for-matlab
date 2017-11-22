function [h, xLabelUnitString, yLabelUnitString] = plot(varargin)
% DimVar.plot plotting method that automatically labels axes with units of
% DimVar inputs.
%   Used with most of the same syntax as plot, but accepts DimVar inputs.
%   Multiple inputs of X and Y must be compatible units.
% 
%   [h, xLabelUnitString, yLabelUnitString] = plot(...) in addition to the
%   handle, returns the strings used by the function for labeling the
%   axes. The string is empty if an axis did not have a DimVar input.
% 
%   Example:
%     edgeL = (0:5)*u.in;
%     squareA = edgeL.^2;
%     triangleA = edgeL.^2*sqrt(3)/4;
%     [~,xStr,yStr] = plot(edgeL,squareA,'s',edgeL,triangleA,'^');
% 
%   See also plot, xlabel, ylabel, text, u, displayingvalue, DimVar.plot3.

%   Copyright 2014 Sky Sartorius
%   www.mathworks.com/matlabcentral/fileexchange/authors/101715 

[varargin,props] = parseplotparams(varargin);

specInd = strncmp('',varargin,0);
args = varargin(~specInd);

nArgs = length(args);
nPairs = nArgs/2;
    
xs = '';
ys = '';

if length(args) == 1
    % Special case of only one argument.
    [args{1},~,ys] = displayparser(args{1});
else
    firstX = args{1};
    if isa(firstX,'DimVar')
        [args{1},~,xs] = displayparser(firstX);
    end
    firstY = args{2};
    if isa(firstY,'DimVar')
        [args{2},~,ys] = displayparser(firstY);
    end
    
    
    % Check to make sure nPairs is an integer.
    if nPairs ~= round(nPairs)
        error('Data must be a single matrix Y or a list of pairs X,Y.')
    end
end

for i = 2:nPairs
    x = args{2*i-1};
    y = args{2*i};
    
    % Check compatibility.
    if isa(x,'DimVar') || isa(firstX,'DimVar')
        compatible(firstX, x);
        args{2*i-1} = displayparser(x);
    end
    if isa(y,'DimVar') || isa(firstY,'DimVar')
        compatible(firstY, y);
        args{2*i} = displayparser(y);
    end
end
   
varargin(~specInd) = args;

s = regexprep({xs ys},{'(' ')'},{'{' '}'});
xLabelUnitString = s{1};
yLabelUnitString = s{2};


try
    h_ = plot(varargin{:},props{:});
    a = gca;
    if ~isempty(xLabelUnitString)
        a.XAxis.TickLabelFormat = ['%g ' xLabelUnitString]; % R2015b+
        %     xlabel(xLabelUnitString) % Prior versions.
    end
    if ~isempty(yLabelUnitString)
        a.YAxis.TickLabelFormat = ['%g ' yLabelUnitString]; % R2015b+
        %     ylabel(yLabelUnitString) % Prior versions.
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

function [args,props] = parseplotparams(args)
% Ignore all arguments from the last char preceded by multiple numerics. See
% also parseparams.
props = {};
for i = numel(args):-1:3
    if ischar(args{i}) && isnumeric(args{i-1}) && isnumeric(args{i-2})
        props = args(i:end);
        args = args(1:i-1);
        break
    end
end
end