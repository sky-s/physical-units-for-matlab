function varargout = plotfunctionwrapper(plotFunction,varargin)
% plotfunctionwrapper(plotFunction,varargin)  Converts all inputs in varargin
% using plottingvalue and passes to plotFunction using feval. If plotFunction is
% a plotting function (i.e., not something like histcounts or contourc),
% plotfunctionwrapper will also add appropriate unit labels to the axes returned
% by gca.
% 
% Units used for plotting are determined by displayUnits and baseUnitSystem, not
% per-variable custom display units. This ensures that multiple elements on the
% same axes are correct.
% 
%   See also plottingvalue, feval, displayUnits,
%     AddSecondAxis - http://www.mathworks.com/matlabcentral/fileexchange/38852,
%     addaxis_unit  - http://www.mathworks.com/matlabcentral/fileexchange/26928.

%% Execute function.
% Convert all DimVar arguments to regular variables.
if numel(varargin) <= 2 && isstruct(varargin{end})
    % Special case of struct input for e.g. patch plotting. Note: because no
    % inputs are DimVar in this case, the overloaded method will not be called,
    % so this block is only here for future/alternative use.
    
    cleanedArgs = varargin;
    S = varargin{end};
    C = cellfun(@plottingvalue,struct2cell(S),'UniformOutput',false);
    cleanedArgs{end} = cell2struct(C,fieldnames(S));
    
else
    cleanedArgs = cellfun(@plottingvalue,varargin,'UniformOutput',false);
    
end
[varargout{1:nargout}] = feval(plotFunction,cleanedArgs{:});

%%
args = varargin;

%% Determine if first input is axes.
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

%% Easy scheme with struct-style input.
[X,Y,Z] = deal([]);
if ischar(args{1}) || isstruct(args{1})
    S = struct(args{:}); % Also works with single input struct.
    
    if isfield(S,'Vertices') 
        % Order is important. Patch will use XData, etc. instead of Vertices if
        % both are present.
        [X,Y,Z] = deal(S.Vertices);
        if size(S.Vertices,2) < 3
            Z = [];
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
    
    labelaxes(X,Y,Z)
    return
end

%% Find just the arguments preceding param/value pairs.
args = parseplotparams(args);

%% Get just the plottable arguments.
plottableArgInd = cellfun(@isplottable,args);
plottableArgs = args(plottableArgInd);
nPlottableArgs = nnz(plottableArgInd);

%% Parse out the intent of the plotting; check compatibility if it's easy.
warnFlag = false;

switch char(plotFunction)
    case {'bar','barh','bar3','bar3h'}
        if isscalar(plottableArgs{end})
            % Last argument is width, so ignore it.
            nPlottableArgs = nPlottableArgs - 1;
        end
        if nPlottableArgs == 2
            % First argument is bar locations.
            switch char(plotFunction)
                case 'bar'
                    labelaxes(plottableArgs{1:2},[])
                case 'barh'
                    labelaxes(plottableArgs{2},plottableArgs{1},[])
                case 'bar3'
                    labelaxes([],plottableArgs{1:2})
                case 'bar3h'
                    labelaxes([],plottableArgs{2},plottableArgs{1})
            end
        else
            switch char(plotFunction)
                case {'bar','bar3h'}
                    labelaxes([],plottableArgs{1},[])
                case 'barh'
                    labelaxes(plottableArgs{1},[],[])
                case 'bar3'
                    labelaxes([],[],plottableArgs{1})
            end
        end
        
    case {'hist','histogram'}
        dimVarArgs = varargin(cellfun('isclass',args,'DimVar'));
        if ~iscompatible(dimVarArgs{:})
            % All DimVar inputs should be compatible.
            warnFlag = true;
        end
        labelaxes(plottableArgs{1},[],[]);
        
    case {'histogram2'}
        labelaxes(plottableArgs{1:2},[])
        
    case {'contour','contourf'}
        if nPlottableArgs >= 3
            labelaxes(plottableArgs{1:2},[])
        end
    
    case {'surf','surface','contour3'}
        if nPlottableArgs <= 2
            % surf(z,c,...); surf(z)
            labelaxes([],[],plottableArgs{1})
            
        else
            % surf(x,y,z); surf(x,y,z,c)
            labelaxes(plottableArgs{1:3})
            
        end
        
    case {'patch'}
        if nPlottableArgs <= 3
            % patch(x,y,c)
            labelaxes(plottableArgs{1:2},[])
            
        else
            % patch(x,y,z,c)
            labelaxes(plottableArgs{1:3})
            
        end
        
    case {'line','text'}
        if nPlottableArgs <= 2
            labelaxes(plottableArgs{1:2},[])
        else
            labelaxes(plottableArgs{1:3})
        end
        
    case {'plot','fill'}
        % Check compatibility.
        if      ~iscompatible(plottableArgs{1:2:end}) || ...
                ~iscompatible(plottableArgs{2:2:end})
            warnFlag = true;
        end
        
        if nPlottableArgs == 1
            labelaxes([],plottableArgs{1},[])
        else
            labelaxes(plottableArgs{1:2},[])
        end
        
    case {'plot3','fill3'}
        % Check compatibility.
        if      ~iscompatible(plottableArgs{1:3:end}) || ...
                ~iscompatible(plottableArgs{2:3:end}) || ...
                ~iscompatible(plottableArgs{3:3:end})
            warnFlag = true;
        end
        
        labelaxes(plottableArgs{1:3})
        
    case {'xline'}
        labelaxes(plottableArgs{1},[],[])
        
    case {'yline'}
        labelaxes([],plottableArgs{1},[])
        
end

%% Send out warning if units might not match.
if warnFlag
    warning('DimVar:plotunitscompatibility',...
        ['Potentially incompatible units in inputs for ' plotFunction '.'])    
end


    function labelaxes(X,Y,Z)
        if providedAxes
            ax = providedAxesHandle;
        else
            ax = gca;
        end
        if isa(X,'DimVar')
            [~,~,~,~,xs] = displayparser(scd(X));
            ax.XAxis.TickLabelFormat = ['%g ' xs]; % R2015b+
        end
        if isa(Y,'DimVar')
            [~,~,~,~,ys] = displayparser(scd(Y));
            % Prepare for the possibility of a yyaxis plot with a 2-element
            % YAxis NumericRuler.
            if strcmp(ax.YAxisLocation,'right') 
                ax.YAxis(end).TickLabelFormat = ['%g ' ys]; % R2015b+
            else
                ax.YAxis(1).TickLabelFormat = ['%g ' ys]; % R2015b+
            end
        end
        if isa(Z,'DimVar')
            [~,~,~,~,zs] = displayparser(scd(Z));
            ax.ZAxis.TickLabelFormat = ['%g ' zs]; % R2015b+
        end
        
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

function out = isplottable(x)
out = isnumeric(x) || isa(x,'datetime') || ...
    isa(x,'duration') || isa(x,'categorical');
end