classdef (InferiorClasses = {?DimVar}) LogDimVar
% LogDimVar  A special case of DimVar for logarithmic physical units
% (e.g., decibels).
%
%   LogDimVars handle units that use logarithmic scales, such as decibels (dB),
%   nepers, and similar units. These units represent ratios on a logarithmic
%   scale. LogDimVars should be used only for setting (through multiplication),
%   converting (through division), and displaying units, so DO NOT USE WITHOUT A
%   MULTIPLICATION OR DIVISION OPERATOR. The math in the background is always
%   operating on the underlying linear value.
% 
%     Use cases:
%       Set a value with multiplication (before or after scalar): 
%           soundLevel = 20*u.dBSPL -> underlying units of pressure
%           gain = u.dB*20          -> yields unitless (double) = 100 since
%                                      reference value is not a DimVar.
%           str2u('20 dBV')         -> a 10 V quantity, displays as dBV
%       Convert units with division:
%           20*u.dBSPL/u.Pa -> 0.0002
%           voltageRatio = 20*u.dBV/u.V -> 10 (voltage ratio)
%       Display as linear instead of dB:
%           12*u.dB (already a double and displays as such)
%           scd(5*u.dBW)
% 
%     Math operations of LogDimVars: These are all DimVars once they are set, so
%     sath is done on the underlying linear units, but in many cases the custom
%     log unit display will carry through the operation. This means that LOG
%     MATH DOESN'T HAPPEN. So, 20*u.dB + 10*u.dB does NOT yield 30 dB, but
%     rather yeld 100 + 10 = 110 (~20.414 dB).
% 
%     Many other use cases should be avoided and mostly throw an error:
%       u.kg*u.dB
%       u.dB*u.dB
%       str2u('20 dB/s')
% 
%   See also DimVar, u, str2u, OffsetDimVar, scd.

    properties (Access = protected)
        customDisplay = ''
    end
    properties
        reference = 1   % reference value with units (DimVar or scalar)
        divisor = 10    % divisor for logarithm (10 for dB power, 20 for dB voltage)
    end
    

    methods
        function v = LogDimVar(reference,divisor)
            % LogDimVar  Construct a logarithmic unit
            %   v = LogDimVar(reference) creates a log unit with default divisor=10
            %     reference should be a DimVar with both value and units.
            %   v = LogDimVar(reference,divisor) sets both reference and divisor.
            if nargin >= 1
                v.reference = reference;
            end
            if nargin >= 2
                v.divisor = divisor;
            end
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
            % Multiplication: scalar * logUnit -> convert from dB to linear scale
            if ~isa(v2,'LogDimVar') % v1 is only LogDimVar
                % scalar * u.dB: interpret scalar as dB value
                linearValue = v1.reference .* 10.^(v2./v1.divisor);
                v = scd(linearValue, v1.customDisplay);
            elseif ~isa(v1,'LogDimVar') % v2 is only LogDimVar
                % u.dB * scalar: interpret scalar as dB value
                linearValue = v2.reference .* 10.^(v1./v2.divisor);
                v = scd(linearValue, v2.customDisplay);
            else % both LogDimVar
                error('LogDimVar:incompatibleUnits',...
                    'Multiplication of two logarithmic units is undefined.')
            end
        end
        
        function v = mtimes(v1,v2)
            v = times(v1,v2);
        end
        
        function v = rdivide(v1,v2)
            % Division: handle conversion from linear to dB scale
            if ~isa(v2,'LogDimVar') % v1 is only LogDimVar
                error('LogDimVar:undefined',...
                    'Division of a logarithmic unit by a scalar is undefined.')
            elseif ~isa(v1,'LogDimVar') % v2 is only LogDimVar
                % linearValue / u.dB: convert linear value to dB
                dbValue = v2.divisor .* log10(v1 ./ v2.reference);
                v = dbValue;
            else % both LogDimVar
                error('LogDimVar:undefined',...
                    'Conversion of a LogDimVar is undefined.')
            end
        end
        
        function v = mrdivide(v1,v2)
            v = rdivide(v1,v2);
        end
        
        function compatible(varargin)
            ME = MException('DimVar:incompatibleUnits',...
                ['Incompatible units. Cannot perform operation on '...
                'variables with different or logarithmic units.']);
            throwAsCaller(ME);
        end
    end
end
