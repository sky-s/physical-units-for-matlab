classdef (InferiorClasses = {?DimVar}) OffsetDimVar
% OffsetDimVar  A special case of DimVar for physical units that have an offset
% zero reference (i.e., Fahrenheit and Celsius temperatures).
%
%   OffsetDimVars should be used only for setting and converting units through
%   multiplication and division. 
% 
%     Use cases:
%       Set a value with multiplication (before or after scalar): 
%           roomTemp = 20*u.degC        -> 293.15 K
%           hotDay = u.degC*40          -> 313.15 K
%           str2u('20 degC')            -> 293.15 K
%       Convert units with division:
%           reallyCold = 200*u.K/u.degC -> -73.15
%           convertToK = 20*u.degC/u.K  -> 293      treated as (20*u.degC)/u.K
%           special case: u.degC/u.degF -> 1.8
% 
%     Pretty much all other use cases should be avoided and mostly throw an
%     error, but not always (usually in cases where it will be interpreted at 1
%     degC = 274.15 K, e.g.), so be careful.
%       u.kg*u.degC
%       u.degC*u.degF
%       u.degC/u.K
%       str2u('20 degC/s')
%       u.degC + 5*u.K 
%       unitconversionfactor(u.K,u.degC)
%       anything else with unitconversionfactor (since it won't be a "factor"),
%       e.g. unitconversionfactor('degC','K'), even though this doesn't throw an
%       error.
%       20*u.degC + 20*u.degF also doesn't error but should be avoided. Use
%       u.deltaDegF instead, for example.
% 
%   See also DimVar, u, str2u, unitconversionfactor, u.deltaDegC, u.deltaDegF.

    properties (Access = protected)
        dv
        customDisplay = ''
    end
    properties
        offset = 0
    end
    

    methods
        function v = OffsetDimVar(dv,offset)
            v.dv = dv;
            v.offset = offset;
        end
        function v = scd(v,val)
            if nargin == 1
                v.customDisplay = '';
            else
                v.customDisplay = val;
            end
        end
        function disp(varargin)
            dispdisp(varargin{:})
        end
        
        function v = times(v1,v2)
            if ~isa(v2,'OffsetDimVar') % v1 is only OffsetDimVar
                v = scd(v2 .* v1.dv + v1.offset,v1.customDisplay);
            elseif ~isa(v1,'OffsetDimVar') % v2 is only OffsetDimVar
                v = scd(v1 .* v2.dv + v2.offset,v2.customDisplay);
            else % both
                error('OffsetDimVar:incompatibleUnits',...
                    'Multiplication should only be used for setting values.')
            end
        end
        function v = mtimes(v1,v2)
            v = times(v1,v2);
        end
        
        function v = rdivide(v1,v2)
            if ~isa(v2,'OffsetDimVar') % v1 is only OffsetDimVar
                error('OffsetDimVar:undefined',...
                    'Division of an offset physical unit is undefined.')
            elseif ~isa(v1,'OffsetDimVar') % v2 is only OffsetDimVar
                v = (v1 - v2.offset)./v2.dv;
            else % both, very special case.
                v = v1.dv./v2.dv;
            end
        end
        function v = mrdivide(v1,v2)
            v = rdivide(v1,v2);
        end
        function compatible(varargin)
            ME = MException('DimVar:incompatibleUnits',...
                ['Incompatible units. Cannot perform operation on '...
                'variables with different or offset units.']);
            throwAsCaller(ME);
        end
    end
end
