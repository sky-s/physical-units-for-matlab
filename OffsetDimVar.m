classdef (InferiorClasses = {?DimVar}) OffsetDimVar
    properties (Access = protected)
        dv
    end
    properties
        offset = 0
    end
    
% Use cases:
%   Set a value: 
%       roomTemp = 20*u.degC -> 293 K
%       hotDay = u.degC*40
%       str2u('20 degC')
%   Convert units:
%       reallyCold = 200*u.K/u.degC -> -73
%       convertToSci = 20*u.degC/u.K -> 293
%       special case: u.degC/u.degF
%   
% Pretty much all other use cases should throw an error:
%   u.kg*u.degC
%   u.degC*u.degF
%   u.degC/u.K
%   anything with unitconversionfactor (since it won't be a "factor")
%   str2u('20 degC/s')
    methods
        function v = OffsetDimVar(dv,offset)
            v.dv = dv;
            v.offset = offset;
        end
        function disp(varargin)
            dispdisp(varargin{:})
        end
        
        function v = times(v1,v2)
            if ~isa(v2,'OffsetDimVar') % v1 is only OffsetDimVar
                v = v2 .* v1.dv + v1.offset;
            elseif ~isa(v1,'OffsetDimVar') % v2 is only OffsetDimVar
                v = v1 .* v2.dv + v2.offset;
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
                    'Dividing by an offset physical unit is undefined.')
            elseif ~isa(v1,'OffsetDimVar') % v2 is only OffsetDimVar
                v = (v1 - v2.offset)./v2.dv;
            else % both, very special case.
                v = v1.dv./v2.dv;
            end
        end
        function v = mrdivide(v1,v2)
            v = rdivide(v1,v2);
        end
    end
end
