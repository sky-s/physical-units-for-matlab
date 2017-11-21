function [h, xLabelString, yLabelString, zLabelString] = plot3(varargin)
% DimVar.plot3  3D plotting method that automatically labels axes with units of
% DimVar inputs.
%   This is the 3D analog of DimVar.plot.
% 
%   See also DimVar.plot.

%   Copyright 2014 Sky Sartorius
%   www.mathworks.com/matlabcentral/fileexchange/authors/101715 

[varargin,props] = parseplotparams(varargin);

specInd = strncmp('',varargin,0);
args = varargin(~specInd);

nArgs = length(args);
nTrips = nArgs/3;
    
xs = '';
ys = '';
zs = '';

firstX = args{1};
if isa(firstX,'DimVar')
    [args{1},~,xs] = displayparser(firstX);
end
firstY = args{2};
if isa(firstY,'DimVar')
    [args{2},~,ys] = displayparser(firstY);
end
firstZ = args{3};
if isa(firstZ,'DimVar')
    [args{3},~,zs] = displayparser(firstZ);
end


% Check to make sure nTrips is an integer.
if nTrips ~= round(nTrips)
    error('Data must be provided in X Y Z sets.')
end


for i = 2:nTrips
    x = args{3*i-2};
    y = args{3*i-1};
    z = args{3*i};
    
    % Check compatibility.
    if isa(x,'DimVar') || isa(firstX,'DimVar')
        compatible(firstX, x);
        args{3*i-2} = displayparser(x);
    end
    if isa(y,'DimVar') || isa(firstY,'DimVar')
        compatible(firstY, y);
        args{3*i-1} = displayparser(y);
    end
    if isa(z,'DimVar') || isa(firstZ,'DimVar')
        compatible(firstZ, z);
        args{3*i} = displayparser(z);
    end
end
   
varargin(~specInd) = args;

s = regexprep({xs ys zs},{'(' ')'},{'{' '}'});
xLabelString = s{1};
yLabelString = s{2};
zLabelString = s{3};

try
    h_ = plot3(varargin{:},props{:});
    a = gca;
    if ~isempty(xLabelString)
        a.XAxis.TickLabelFormat = ['%g ' xLabelString]; % R2015b+
        %     xlabel(xLabelUnitString) % Prior versions.
    end
    if ~isempty(yLabelString)
        a.YAxis.TickLabelFormat = ['%g ' yLabelString]; % R2015b+
        %     ylabel(yLabelUnitString) % Prior versions.
    end
    if ~isempty(zLabelString)
        a.ZAxis.TickLabelFormat = ['%g ' zLabelString]; % R2015b+
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