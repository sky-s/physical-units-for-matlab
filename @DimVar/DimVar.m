classdef (InferiorClasses = {?matlab.graphics.axis.Axes}) DimVar
% See also u.

% Copyright (c) 2012-2017, Sky Sartorius.
properties (Access = protected)
    exponents
    value
end
properties
    customDisplay
end

methods
    function v = scd(v,val)
        % scd  Set Custom Display units on a per-variable basis. 
        %   v = scd(v,str) uses str as the preferred custom display unit for v.
        %   str must be a valid field of u or be evaluable by str2u.
        %
        %   v = scd(v) with only one input returns v with custom display units
        %   cleared. Custom display units are also cleared by most operations
        %   that change the units.
        %   
        %   See also str2u, displayUnits.
        if nargin == 1
            v.customDisplay = '';
        else
            v.customDisplay = val;
        end
    end
    %% Core methods (not overloads).
    % Constructor:
    function v = DimVar(expos,val)
        % See also u.
        v.exponents = expos;
        v.value = val;
    end
    
    function v = clearcanceledunits(v)
        % If all DimVar unit exponents are zero, return normal (double)
        % variable. Exponent tolerance is to fifth decimal.
        
        if ~any(round(1e5*v.exponents)) % Seems to be faster than round(x,5).
            v = v.value;
        else
            v.customDisplay = '';
            % The customDisplay property is invalid after e.g. a multiply or
            % divide operation, so clean it up from the new variable to prevent
            % any undesirable side effects later.
        end
    end
    
    function compatible(v,varargin)
        % See also compatible, iscompatible.
        
        % Note: This must be both a function and a DimVar method so that the
        % method can access the exponents property.
        
        if ~isa(v,'DimVar')
            
            ME = MException('DimVar:incompatibleUnits',...
                ['Incompatible units. Cannot perform operation on '...
                'variables with different units.']);
            throwAsCaller(ME);
            
        end
        
        vExpos = v.exponents;
        
        for i = 1:numel(varargin)
            
            if ~isa(varargin{i},'DimVar') ...
                    || ~isequal(vExpos,varargin{i}.exponents)
                
                ME = MException('DimVar:incompatibleUnits',...
                    ['Incompatible units. Cannot perform operation on '...
                    'variables with different units.']);
                throwAsCaller(ME);
                
            end
        end
    end
    
    %% Concatenation.
    function v = cat(dim,v,varargin)
        if ~isa(v,'DimVar')
            ME = MException('DimVar:incompatibleUnits',...
                'Incompatible units. All inputs must be DimVar.');
            throwAsCaller(ME);
        end
        
        vExpos = v.exponents;
        
        for i = 1:numel(varargin)
            vi = varargin{i};
            
            if ~isa(vi,'DimVar') || ~isequal(vExpos,vi.exponents)
                ME = MException('DimVar:incompatibleUnits',...
                    ['Incompatible units. Cannot perform operation on '...
                    'variables with different units.']);
                throwAsCaller(ME);
            end
            
            v.value = cat(dim,v.value,vi.value);
        end
        
        % Functionality of compatible method is integrated for sake of speed.
    end
    function vOut = horzcat(varargin)
        vOut = cat(2,varargin{:});
    end
    function vOut = vertcat(varargin)
        vOut = cat(1,varargin{:});
    end
    
    %% Validation functions (mustBe__).
    function mustBeGreaterThan(v1,v2)
        compatible(v1,v2);
        mustBeGreaterThan(v1.value,v2.value);
    end
    function mustBeGreaterThanOrEqual(v1,v2)
        compatible(v1,v2);
        mustBeGreaterThanOrEqual(v1.value,v2.value);
    end
    function mustBeLessThan(v1,v2)
        compatible(v1,v2);
        mustBeLessThan(v1.value,v2.value);
    end
    function mustBeLessThanOrEqual(v1,v2)
        compatible(v1,v2);
        mustBeLessThanOrEqual(v1.value,v2.value);
    end
    function mustBeNegative(v);     mustBeNegative(v.value);    end
    function mustBeNonnegative(v);  mustBeNonnegative(v.value); end
    function mustBeNonpositive(v);  mustBeNonpositive(v.value); end
    function mustBeNonzero(v);      mustBeNonzero(v.value);     end
    function mustBePositive(v);     mustBePositive(v.value);    end
    
    %% is__ functions.
    function result = isempty(v);   result = isempty(v.value);      end
    function result = isfinite(v);  result = isfinite(v.value);     end
    function result = isinf(v);     result = isinf(v.value);        end
    function result = isnan(v);     result = isnan(v.value);        end
    function result = isnumeric(v); result = isnumeric(v.value);    end
    function result = isreal(v);    result = isreal(v.value);       end
    
    %% Logical operators (>, <, ==, ~, etc.).
    function result = eq(v1,v2)
        compatible(v1,v2);
        result = v1.value == v2.value;
    end
    function result = ge(v1,v2)
        compatible(v1,v2);
        result = v1.value >= v2.value;
    end
    function result = gt(v1,v2)
        compatible(v1,v2);
        result = v1.value > v2.value;
    end
    function result = le(v1,v2)
        compatible(v1,v2);
        result = v1.value <= v2.value;
    end
    function result = lt(v1,v2)
        compatible(v1,v2);
        result = v1.value < v2.value;
    end
    function result = ne(v1,v2)
        compatible(v1,v2);
        result = v1.value ~= v2.value;
    end
    function result = not(v)
        result = ~v.value;
    end
    
    %% Class conversions.
    function result = double(v)
        % DimVar.double(V)  Returns value property (not diplayingvalue) of V.
        %
        % See also u2num, displayingvalue, displayparser.
        result = v.value;
    end
    function result = logical(v)
        result = logical(v.value);
    end
    function s = string(v)
        s = string(cellfun(@num2str,num2cell(v),'UniformOutput',false));
    end
    
    %% Simple functions that return DimVar.
    function v = abs(v);            v.value = abs(v.value);                 end
    function v = circshift(v,varargin)
        v.value = circshift(v.value,varargin{:});
    end
    function v = conj(v);           v.value = conj(v.value);                end
    function v = ctranspose(v);     v.value = v.value';                     end
    function v = cumsum(v,varargin);v.value = cumsum(v.value,varargin{:});  end
    function v = diag(v,varargin);  v.value = diag(v.value,varargin{:});    end
    function v = diff(v,varargin);  v.value = diff(v.value,varargin{:});    end
    function v = full(v);           v.value = full(v.value);                end
    function v = imag(v);           v.value = imag(v.value);                end
    function v = mean(v,varargin);  v.value = mean(v.value,varargin{:});    end
    function v = median(v,varargin);v.value = median(v.value,varargin{:});  end
    function v = norm(v,varargin);  v.value = norm(v.value,varargin{:});    end
    function v = permute(v,varargin);v.value = permute(v.value,varargin{:});end
    function v = real(v);           v.value = real(v.value);                end
    function v = reshape(v,varargin);v.value = reshape(v.value,varargin{:});end
    function [v,I] = sort(v,varargin)
        [v.value, I] = sort(v.value,varargin{:});
    end
    function v = std(v,varargin);   v.value = std(v.value,varargin{:});     end
    function v = subsref(v,varargin);v.value = subsref(v.value,varargin{:});end
    function v = sum(v,varargin);   v.value = sum(v.value,varargin{:});     end
    function v = trace(v);          v.value = trace(v.value);               end
    function v = transpose(v);      v.value = v.value.';                    end
    function v = tril(v,varargin);  v.value = tril(v.value,varargin{:});    end
    function v = triu(v,varargin);  v.value = triu(v.value,varargin{:});    end
    function v = uminus(v);         v.value = -v.value;                     end
    function v = uplus(v);                                                  end
    
    %% Functions that require compatibility check.
    function out = atan2(v1,v2)
        compatible(v1,v2);
        out = atan2(v1.value, v2.value);
    end
    function v1 = hypot(v1,v2)
        compatible(v1,v2);
        v1.value = hypot(v1.value,v2.value);
    end
    function v1 = minus(v1,v2)
        compatible(v1,v2);
        v1.value = v1.value - v2.value;
    end
    function v1 = plus(v1,v2)
        compatible(v1,v2);
        v1.value = v1.value + v2.value;
    end
    
    %% Other simple functions.
    function n = length(v);     n = length(v.value);    end
    function n = ndims(v);      n = ndims(v.value);     end
    function n = nnz(v);        n = nnz(v.value);       end
    function n = numel(v);      n = numel(v.value);     end
    function out = sign(v);     out = sign(v.value);    end
    function varargout = size(x,varargin)
        [varargout{1:nargout}] = size(x.value,varargin{:});
    end
    function out = issorted(v,varargin)
        out = issorted(v.value,varargin{:});
    end
    
    %% Plot-like functions.
    function varargout = bar(varargin)
        [varargout{1:nargout}] = plotfunctionwrapper('bar',varargin{:});
    end
    function varargout = bar3(varargin)
        [varargout{1:nargout}] = plotfunctionwrapper('bar3',varargin{:});
    end
    function varargout = barh(varargin)
        [varargout{1:nargout}] = plotfunctionwrapper('barh',varargin{:});
    end
    function varargout = bar3h(varargin)
        [varargout{1:nargout}] = plotfunctionwrapper('bar3h',varargin{:});
    end
    function varargout = contour(varargin)
        [varargout{1:nargout}] = plotfunctionwrapper('contour',varargin{:});
    end
    function varargout = contour3(varargin)
        [varargout{1:nargout}] = plotfunctionwrapper('contour3',varargin{:});
    end
    function varargout = contourf(varargin)
        [varargout{1:nargout}] = plotfunctionwrapper('contourf',varargin{:});
    end
    function varargout = contourc(varargin)
        [varargout{1:nargout}] = plotfunctionwrapper('contourc',varargin{:});
    end
    function varargout = fill(varargin)
        [varargout{1:nargout}] = plotfunctionwrapper('fill',varargin{:});
    end
    function varargout = fill3(varargin)
        [varargout{1:nargout}] = plotfunctionwrapper('fill3',varargin{:});
    end
    function varargout = hist(varargin)
        [varargout{1:nargout}] = plotfunctionwrapper('hist',varargin{:});
    end
    function varargout = histcounts(varargin)
        [varargout{1:nargout}] = plotfunctionwrapper('histcounts',varargin{:});
    end
     function varargout = histcounts2(varargin)
        [varargout{1:nargout}] = plotfunctionwrapper('histcounts2',varargin{:});
    end
    function varargout = histogram(varargin)
        [varargout{1:nargout}] = plotfunctionwrapper('histogram',varargin{:});
    end
    function varargout = histogram2(varargin)
        [varargout{1:nargout}] = plotfunctionwrapper('histogram2',varargin{:});
    end
    function varargout = line(varargin)
        [varargout{1:nargout}] = plotfunctionwrapper('line',varargin{:});
    end
    % TODO: add mesh, meshc, meshz.
    function varargout = patch(varargin)
        [varargout{1:nargout}] = plotfunctionwrapper('patch',varargin{:});
    end
    function varargout = plot(varargin)
        [varargout{1:nargout}] = plotfunctionwrapper('plot',varargin{:});
    end
    function varargout = plot3(varargin)
        [varargout{1:nargout}] = plotfunctionwrapper('plot3',varargin{:});
    end
    % TODO: add ribbon, slice.
    function varargout = surf(varargin)
        [varargout{1:nargout}] = plotfunctionwrapper('surf',varargin{:});
    end
    % TODO: add surfc, surfl.
    function varargout = surface(varargin)
        [varargout{1:nargout}] = plotfunctionwrapper('surface',varargin{:});
    end
    function varargout = text(varargin)
        [varargout{1:nargout}] = plotfunctionwrapper('text',varargin{:});
    end
    function varargout = xline(varargin)
        [varargout{1:nargout}] = plotfunctionwrapper('xline',varargin{:});
    end
    function varargout = yline(varargin)
        [varargout{1:nargout}] = plotfunctionwrapper('yline',varargin{:});
    end
    
    %% Axis limits _lim functions.
    function xlim(varargin)
        varargin = cellfun(@plottingvalue,varargin,'UniformOutput',false);
        xlim(varargin{:});
    end
    function ylim(varargin)
        varargin = cellfun(@plottingvalue,varargin,'UniformOutput',false);
        ylim(varargin{:});
    end
    function zlim(varargin)
        varargin = cellfun(@plottingvalue,varargin,'UniformOutput',false);
        zlim(varargin{:});
    end

% Certain saveobj and loadobj functionality may be desirable to turn on in some
% circumstances.

%     function s = saveobj(v)
%         s = struct('exponents',v.exponents,...
%             'value',v.value,'customDisplay',v.customDisplay,...
%             'base',u.baseUnitSystem);
%     end
end

methods (Static)
    function v = loadobj(v)
        % Loading DimVar objects that may have been saved with an older version
        % of the Physical Units Toolbox that had a shorter vector of base units.
        % This will have undesirable behavior if the saved and loaded variables
        % do not share the same base unit system.
        
        base = u.baseUnitSystem;
        
        ex = v.exponents;
        nEx = numel(ex);
        nBase = length(base);
        if nEx < nBase
            v.exponents = zeros(1,nBase);
            v.exponents(1:nEx) = ex;
        end
        
%         if isa(v,'DimVar')
%             return
%         end
%         
%         v_ = DimVar(v.exponents,v.value);
%         
%         if isfield(v,'customDisplay')
%             v_ = scd(v_,v.customDisplay);
%         end
%         if isfield(v,'base')
%             % Check that base unit system (at least the multipliers) is the
%             % same.
%             if ~isequal(base(1:nEx,2),v.base(1:nEx,2))
%                 warning(['Base unit system of loaded object does not match '...
%                 'current workspace base units system.'])
%             end
%         end
    end
end
end