function [h, xLabelString, yLabelString, zLabelString] = plot3(varargin)
% DimVar.plot3  3D plotting method that automatically labels axes with units of
% DimVar inputs.
%   This is the 3D analog of DimVar.plot.
% 
%   See also DimVar.plot.

%   Copyright 2014 Sky Sartorius
%   www.mathworks.com/matlabcentral/fileexchange/authors/101715 

specInd = strncmp('',varargin,0);
args = varargin(~specInd);

nArgs = length(args);
nTrips = nArgs/3;
    
xs = '';
ys = '';
zs = '';

firstX = args{1};
if isa(firstX,'DimVar')
    [~, xs, firstX] = num2str(firstX);
    args{1} = firstX.value;
end
firstY = args{2};
if isa(firstY,'DimVar')
    [~, ys, firstY] = num2str(firstY);
    args{2} = firstY.value;
end
firstZ = args{3};
if isa(firstZ,'DimVar')
    [~, zs, firstZ] = num2str(firstZ);
    args{3} = firstZ.value;
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
    if isa(x,'DimVar')
        compatible(firstX, x);
        [~, ~, x] = display(x);
        args{3*i-2} = x.value;
    end
    if isa(y,'DimVar')
        compatible(firstY, y);
        [~, ~, y] = display(y);
        args{3*i-1} = y.value;
    end
    if isa(z,'DimVar')
        compatible(firstZ, z);
        [~, ~, z] = display(z);
        args{3*i} = z.value;
    end
end
   
varargin(~specInd) = args;

s = regexprep({xs ys zs},{'(' ')'},{'{' '}'});
xLabelString = s{1};
yLabelString = s{2};
zLabelString = s{3};

try
    h_ = plot3(varargin{:});
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
    
catch
    warning(['DimVar.plot in development and designed for R2015b. '...
        'Convert DimVar to double then use normal plot function.'])
end

if nargout
    h = h_;
end

end