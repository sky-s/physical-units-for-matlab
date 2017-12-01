function varargout = plotfunctionwrapper(plotFunction,varargin)
% plotFunction is string or handle for feval.

%% Possible formats
%{
surf(z,c,...); surf(z)

plot, plot3, fill, and fill3 can't take name/value pairs (XData, YData, etc.),
so these simply get data labels based on the first plottable arguments.


things that might want dimvar in follow-on param/value pairs:
histogram, histogram2

nothing else that I know of:
patch
plot
surf
fill
%}

%% Execute function.
% Convert all DimVar arguments to regular variables.
cleanedArgs = cellfun(@displayingvalue,varargin,'UniformOutput',false);
[varargout{1:nargout}] = feval(plotFunction,cleanedArgs{:});
ax = gca;

%%
args = varargin;

%% Determine if first input is axes.
if numel(args) && ...
        ((isscalar(args{1}) && isgraphics(args{1},'axes')) ...
        || isa(args{1},'matlab.graphics.axis.AbstractAxes') ...
        || isa(args{1},'matlab.ui.control.UIAxes'))
    
    args = varargin(2:end);
    
else
    args = varargin;
    
end

%% Easy scheme.
[X,Y,Z] = deal([]);
if isstruct(args{1}) || ischar(args{1})
    S = struct(args{:}); % Works with single input struct.
    
    % struct input scheme.
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
    
    labelaxes(ax,X,Y,Z)
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
    case {'hist','histogram'}
        dimVarArgs = varargin(cellfun('isclass',args,'DimVar'));
        if ~iscompatible(dimVarArgs{:})
            % All DimVar inputs should be compatible.
            warnFlag = true;
        end
        labelaxes(ax,plottableArgs{1},[],[]);
        
    case {'histogram2'}
        labelaxes(ax,plottableArgs{1:2},[])
        
    case {'surf','surface'}
        if nPlottableArgs <= 2
            % surf(z,c,...); surf(z)
            labelaxes(ax,[],[],plottableArgs{1})
            
        else
            % surf(x,y,z); surf(x,y,z,c)
            labelaxes(ax,plottableArgs{1:3})
            
        end
        
    case {'patch'}
        if nPlottableArgs <= 3
            % patch(x,y,c)
            labelaxes(ax,plottableArgs{1:2},[])
            
        else
            % patch(x,y,z,c)
            labelaxes(ax,plottableArgs{1:3})
            
        end
        
    case {'line'}
        if nPlottableArgs <= 2
            labelaxes(ax,plottableArgs{1:2},[])
        else
            labelaxes(ax,plottableArgs{1:3})
        end
        
    case {'plot','fill'}
        % Check compatibility.
        if      ~iscompatible(plottableArgs{1:2:end}) || ...
                ~iscompatible(plottableArgs{2:2:end})
            warnFlag = true;
        end
        
        labelaxes(ax,plottableArgs{1:2},[])
        
    case {'plot3','fill3'}
        % Check compatibility.
        if      ~iscompatible(plottableArgs{1:3:end}) || ...
                ~iscompatible(plottableArgs{2:3:end}) || ...
                ~iscompatible(plottableArgs{3:3:end})
            warnFlag = true;
        end
        
        labelaxes(ax,plottableArgs{1:3})
        
end

%% Send out warning if units might not match.
if warnFlag
    warning('DimVar:plotunitscompatibility',...
        ['Potentially incompatible units in inputs for ' plotFunction '.'])    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% Parse out the intent of the plotting; check compatibility if it's easy.
% warnFlag = false;
% 
% switch char(plotFunction)
%     case {'hist','histogram','histcounts'}
%         % All DimVar inputs should be compatible.
%         if ~iscompatible(dimVarArgs{:})
%             warnFlag = true;
%         end
%         labelaxes(ax,plottableArgs{1},[],[]);
%         
%     case {'histogram2','histcounts2'}
%         if ~iscompatible(plottableArgs{1:2:end}) ||...
%                 ~iscompatible(plottableArgs{2:2:end})
%             warnFlag = true;
%         end
%         labelaxes(ax,plottableArgs{1:2},[])
%         
%     case {'surf','surface'}
%         if nPlottableArgs <= 2
%             % surf(z,c,...); surf(z)
%             labelaxes(ax,[],[],plottableArgs{1})
%             
%         else
%             % surf(x,y,z); surf(x,y,z,c)
%             labelaxes(ax,plottableArgs{1:3})
%             
%         end
%         
%     case {'patch'}
%         % (x,y,c),(x,y,z,c), 
%         if nPlottableArgs == 0
%             % Assume struct input.
%         elseif nPlo
%         end
% end       
%% Label axes

% if ~isempty(ys)
%     ax.YAxis.TickLabelFormat = ['%g ' ys]; % R2015b+
%     %     ylabel(ys) % Prior versions.
% end
% if ~isempty(zs)
%     ax.ZAxis.TickLabelFormat = ['%g ' zs]; % R2015b+
%     %     zlabel(zs) % Prior versions.
% end


% switch nDimVarArgs
% for i = 1:nDimVars
%     ind = dimVarArgLocations(i);
%     
%     
% end


% should work with plot(Y), plot(X,Y), surf(Z), surf(X,Y,Z), ...
% patch, surf, histogram, plot, plot3, histogram2, 


% axescheck
% parseparams
% (InferiorClasses = {?matlab.graphics.axis.Axes}) DimVAr
% 
% function [newArgList,labels] = parseplotunits(varargin)
% 
end

function labelaxes(ax,X,Y,Z)
if isa(X,'DimVar')
    [~,~,~,~,~,xs] = displayparser(X);
    ax.XAxis.TickLabelFormat = ['%g ' xs]; % R2015b+
end
if isa(Y,'DimVar')
    [~,~,~,~,~,ys] = displayparser(Y);
    ax.YAxis.TickLabelFormat = ['%g ' ys]; % R2015b+
end
if isa(Z,'DimVar')
    [~,~,~,~,~,zs] = displayparser(Z);
    ax.ZAxis.TickLabelFormat = ['%g ' zs]; % R2015b+
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