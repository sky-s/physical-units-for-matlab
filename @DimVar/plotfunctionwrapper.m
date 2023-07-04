function varargout = plotfunctionwrapper(plotFunction,varargin)
% plotfunctionwrapper(plotFunction,varargin)  Converts all inputs in varargin
% using plottingvalue and passes to plotFunction using feval. If plotFunction is
% a plotting function (i.e., not something like or contourc),
% plotfunctionwrapper will also add appropriate unit labels to the axes returned
% by gca.
% 
% Units used for plotting are determined first by per-variable custom display
% units. 
% 
%   See also plottingvalue, feval, displayUnits,
%     AddSecondAxis - http://www.mathworks.com/matlabcentral/fileexchange/38852,
%     addaxis_unit  - http://www.mathworks.com/matlabcentral/fileexchange/26928.

% TODO: check on all functions using harmonize to make sure it's not
% over-capturing other numeric inputs such as LineWidth or color. Write test
% cases.

%% Determine if first input is axes.
args = varargin;
if numel(args) && ...
        ((isscalar(args{1}) && isgraphics(args{1},'axes')) ...
        || isa(args{1},'matlab.graphics.axis.AbstractAxes') ...
        || isa(args{1},'matlab.ui.control.UIAxes'))
    
    args = varargin(2:end);
    providedAxes = true;
    providedAxesHandle = varargin{1};
    
else
    args = varargin;
    providedAxes = false;
end

%% Get a representative X, Y, and Z.
[X,Y,Z] = deal({});

%% Easy scheme with all name-value input, indicated by char first arg

if ischar(args{1}) || isstring(args{1}) || isstruct(args{1})
    structInput = true;
    % Also works with single input struct if this ever got called for that, but
    % it won't be.
    S = struct(args{:});
    
    if isfield(S,'Vertices') 
        % Order is important. Patch will use XData, etc. instead of Vertices if
        % both are present.
        [X,Y,Z] = deal(S.Vertices);
        if size(S.Vertices,2) < 3
            Z = {};
        end
    end
    
    if isfield(S,'XData')
        X = S.XData;
    end
    if isfield(S,'YData')
        Y = S.YData;
    end
    if isfield(S,'ZData')
        Z = S.ZData;
    end
    
    nvArgs = {};
else
    structInput = false;
    [args,nvArgs] = parseplotparams(args);
end

%% Work with just the arguments preceding param/value pairs.

%% Get just the plottable arguments.
plottableArgInd = cellfun(@isplottable,args);
plottableArgs = args(plottableArgInd);
nPlottableArgs = nnz(plottableArgInd);

%% Parse out the intent of the plotting; harmonize if needed.
nonPlottingFunc = false;
if ~structInput
    switch char(plotFunction)
        case {'area'}
            if nPlottableArgs > 1 && isscalar(plottableArgs{end}) && ~((nPlottableArgs == 2) ...
                    && ~isscalar(plottableArgs{1}) && ~isscalar(plottableArgs{2}))
                % f(X,Y,base); f(Y,base)
                plottableArgs(end-1:end) = harmonize(plottableArgs(end-1:end));
                Y = plottableArgs{end-1};
                if nPlottableArgs > 2
                    X = plottableArgs{1};
                end
            elseif nPlottableArgs == 1
                % f(Y)
                Y = plottableArgs{1};
            else
                % f(X,Y)
                [X,Y] = deal(plottableArgs{1:2});
            end

        case {'bar','barh','bar3','bar3h'}
            if isscalar(plottableArgs{end})
                % Last argument is width, so ignore it.
                nPlottableArgs = nPlottableArgs - 1;
            end
            if nPlottableArgs == 2
                % First argument is bar locations.
                switch char(plotFunction)
                    case 'bar'
                        [X,Y] = deal(plottableArgs{1:2});
                    case 'barh'
                        X = plottableArgs{2};
                        Y = plottableArgs{1};
                    case 'bar3'
                        [Y,Z] = deal(plottableArgs{1:2});
                    case 'bar3h'
                        Y = plottableArgs{2};
                        Z = plottableArgs{1};
                end
            else
                switch char(plotFunction)
                    case {'bar','bar3h'}
                        Y = plottableArgs{1};
                    case 'barh'
                        X = plottableArgs{1};
                    case 'bar3'
                        Z = plottableArgs{1};
                end
            end
            
        case {'hist','histogram'}
            % Arguments that will break this: Color or Alpha of edges/faces,
            % LineWidth.
            if numel(args) > 1 && isscalar(args(2)) && isnumeric(args(2)) && ~isa(args(2),'DimVar')
                % Argument is nbins, which is the only numeric that is allowed
                % to be non-compatible.
                plottableArgs([1,3:end]) = harmonize(plottableArgs([1,3:end]));
            else
                plottableArgs = harmonize(plottableArgs);
            end
            X = plottableArgs{1};
            
        case {'histogram2'}
            [X,Y] = deal(plottableArgs{1:2});
            
        case {'contour','contourf'}
            if nPlottableArgs >= 3
                [X,Y] = deal(plottableArgs{1:2});
            end
            
        case {'surf','surfc','surfl','surface','contour3','mesh','meshc','meshz'}
            if nPlottableArgs <= 2
                % surf(z,c,...); surf(z)
                Z = plottableArgs{1};
                
            else
                % surf(x,y,z); surf(x,y,z,c)
                [X,Y,Z] = deal(plottableArgs{1:3});
                
            end
        
        case {'scatter','bubblechart','swarmchart'}
            % f(x,y); f(x,y,sz); f(x,y,sz,c)
            [X,Y] = deal(plottableArgs{1:2});

        case {'scatter3','bubblechart3','swarmchart3'}
            % f(x,y,z); f(x,y,z,sz); f(x,y,z,sz,c)
            [X,Y,Z] = deal(plottableArgs{1:3});

        case {'slice'}
            % f(x,y,z,V,Sx,Sy,Sz); f(V,Sx,Sy,Sz)
            if nPlottableArgs == 7
                [X,Y,Z] = deal(plottableArgs{1:3});
                for i = 1:3
                    plottableArgs([i,i+4]) = harmonize(plottableArgs([i,i+4]));
                end
            elseif nPlottableArgs == 4
                [X,Y,Z] = deal(plottableArgs{2:4});
            else
                error('Case not covered for DimVar inputs.')
            end

        case {'isosurface'}
            if nPlottableArgs > 3
                % f(x,y,z,v,isoval); f(x,y,z,v); f(...,colors)
                [X,Y,Z] = deal(plottableArgs{1:3});
                % Because isosurface uses a single double array to define output
                % vertices, all cardinal coordinates must be harmonized.
                plottableArgs(1:3) = harmonize(plottableArgs(1:3));
            % else
                % f(v,isoval); f(v); f(...,colors)
            end
            if nargout
                nonPlottingFunc = true;
            end

        case {'patch'}
            if nPlottableArgs <= 3
                % patch(x,y,c)
                [X,Y] = deal(plottableArgs{1:2});
                
            else
                % patch(x,y,z,c)
                [X,Y,Z] = deal(plottableArgs{1:3});
                
            end
            
        case {'line','text'}
            if nPlottableArgs <= 2
                [X,Y] = deal(plottableArgs{1:2});
            else
                [X,Y,Z] = deal(plottableArgs{1:3});
            end
            
        case {'plot','fill'}
            if nPlottableArgs == 1
                Y = plottableArgs{1};
            else
                X = plottableArgs{1};
                Y = plottableArgs{2};
                plottableArgs(1:2:end) = harmonize(plottableArgs(1:2:end));
                plottableArgs(2:2:end) = harmonize(plottableArgs(2:2:end));
            end
            
        case {'plot3','fill3'}
            plottableArgs(1:3:end) = harmonize(plottableArgs(1:3:end));
            plottableArgs(2:3:end) = harmonize(plottableArgs(2:3:end));
            plottableArgs(3:3:end) = harmonize(plottableArgs(3:3:end));
            
            [X,Y,Z] = deal(plottableArgs{1:3});
            
        case {'ribbon'}
            if nPlottableArgs == 1
                Y = plottableArgs{1};
            else
                [X,Y] = deal(plottableArgs{1:2});
            end

        case {'xline','xlim'}
            X = plottableArgs{1};
            
        case {'yline','ylim'}
            Y = plottableArgs{1};
            
        case {'zlim'}
            Z = plottableArgs{1};
            
        otherwise
            nonPlottingFunc = true;
            % Does this want to clean the units into base units, i.e. doubles?
            % Use cases are just contourc, histc
    end
    
end
%% Rebuild (consistent) arguments and run plotFunction.
args(plottableArgInd) = cellfun(@plottingvalue,plottableArgs,'UniformOutput',false);

if providedAxes
    newArgList = [{providedAxesHandle}, args, nvArgs];
else
    newArgList = [args, nvArgs];
end
[varargout{1:nargout}] = feval(plotFunction,newArgList{:});

if nonPlottingFunc
    % TODO: consider adding a step here to add units to the output vertices of
    % isosurface and perhaps other applicable functions.
    return
end

%% If a function that plots, make or check matching axis units; scenarios:

% The plotFunction doesn't actually plot or do anything with axes.

% Existing axes, but are blown away by combination of plotFunction type, hold
% state, etc., whether or not the target axes are explicitly provided.

% Existing axes have no data and are retained; lack of existing unit labels
% indicator of a blank slate.

% Prior axes have data that is retained, including units that will be used to
% enforce consistency.

% Prior axes have data and are retained, but axes are 'unitless', so very
% difficult to detect if this should throw an error for trying to plot DimVars
% onto it. In this case, old data is retained and the user is responsible for
% have previously plotted in the intended units. Similarly, putting e.g. a
% non-DimVar line onto a united plot, the user is expected to be working in the
% correct units.

% No existing axes, un-unitted axes created by plotFunction.

% Summary: use the existence of axis units AFTER executing plotFunction as
% indicator that inconsistency with prior axes should throw a warning.

if providedAxes
    ax = providedAxesHandle;
else
    % It's important to have the plot function go before any potential call to
    % gca.
    ax = gca; 
end

if ax.Type == "polaraxes"
    rulers = [ax.ThetaAxis,ax.RAxis];
    % X and Y are Theta and R.
else
    if ax.YAxisLocation == "left"
        yRuler = ax.YAxis(1);
    else
        yRuler = ax.YAxis(end);
    end
    rulers = [ax.XAxis,yRuler,ax.ZAxis];
end
    
%% Throw warnings, add labels
% if plotting onto the axis inconsistently (only if the plotFunction uses this
% axis, allowing for e.g. 2D text on existing 3D axes.)

% TODO: if ever revamping this, strongly consider the use of setappdata and
% getappdata to enforce consistency on axes. This is much more elegant than
% reading TickLabelFormat (though may want to couple both together). Regardless,
% it should somehow be documented, as there must have been a good reason the
% scheme of using TickLabelFormat was used instead of adding properties to the
% axes or using UserData.

labelUnits = regexp(get(rulers,'TickLabelFormat')','%g (?<unit>.+)','names');

if all(cellfun('isempty',labelUnits))
    % No prior axes with units, so don't check against them.
    enforcing = false;
else
    % Prior axes had units, so check for consistency and/or set units.
    enforcing = true;
end
    

data = {X,Y,Z};
for i = 1:numel(rulers)
    if (iscell(data{i}) && isempty(data{i})) ...
            || ~isa(rulers(i),'matlab.graphics.axis.decorator.NumericRuler')
        % Axis unused by plottingFunction or is not numeric/DimVar, so skip.
    else
        thisData = data{i};
        dataHasUnits = isa(thisData,'DimVar');
        if dataHasUnits
            [~,~,~,~,s] = displayparser(thisData);
        end
        
        %% Throw warnings for incompatibility
        if enforcing
            axisHasUnits = ~isempty(labelUnits{i});
            xyz = 'XYZ';
            wStr = 'Plotting a %s variable onto an established %s axis with %s.';
            if axisHasUnits && dataHasUnits
                % Check for matching plotting units.
                if strcmp(labelUnits{i}.unit, s)
                    % Good, quick, low-overhead identification of match.
                else
                    u1 = str2u(labelUnits{i}.unit);
                    u2 = str2u(s);
                    if ~iscompatible(u1,u2)
                        % Check OffsetDimVar special case.
                        if ~(isa(u1,'OffsetDimVar') && isa(u2,'OffsetDimVar') ...
                                && u1.offset == u2.offset)
                            warning('DimVar:incompatiblePlottingUnits',wStr,...
                                'dimensioned',xyz(i),'incompatible units');
                        end
                    elseif u1 ~= u2
                        warning('DimVar:incompatiblePlottingUnits',wStr,...
                            'dimensioned',xyz(i),'non-matching display units')
                        
                    %else
                        % u1 == u2% Also good, e.g. Mg and tonne.
                    end
                end
                
            elseif dataHasUnits && ~axisHasUnits
                if i == 2 && numel(ax.YAxis) == 2
                    % Assume new yyaxis, so don't warn.
                else
                    warning('DimVar:incompatiblePlottingUnits',wStr,'dimensioned',xyz(i),'no units')
                end
            elseif axisHasUnits && ~dataHasUnits
                warning('DimVar:incompatiblePlottingUnits',wStr,'non-dimensioned',xyz(i),'units')
            else %if ~axisHasUnits && ~dataHasUnits
                % All good, no units.
            end
        end
        
        %% Set TickLabelFormat
        if dataHasUnits
            rulers(i).TickLabelFormat = ['%g ' s];
        else
            rulers(i).TickLabelFormat = '%g';
        end
    end
end

end

function args = harmonize(args)
% Consider only enforcing this on DimVars, potentially allowing more errors through but enables more
% name/value arguments to functions.
compatible(args{:});
if isa(args{1},'DimVar')
    dominantUnit = args{1}.customDisplay;
    args = cellfun(@(x) scd(x,dominantUnit),args,'UniformOutput',0);
end
end

function [args,props] = parseplotparams(args)
% Ignore all arguments from the last char preceded by multiple numerics. See
% also parseparams.
props = {};
for i = numel(args):-1:3
    if (ischar(args{i}) || isstring(args{i})) && isnumeric(args{i-1}) && isnumeric(args{i-2})
        props = args(i:end);
        args = args(1:i-1);
        break
    end
end
end

function out = isplottable(x)
out = isnumeric(x) || isa(x,'datetime') || isa(x,'duration') || isa(x,'categorical');
end
