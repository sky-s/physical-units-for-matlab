classdef (InferiorClasses = {?DimVar}) LogDimVar
% LogDimVar  A special case of DimVar for logarithmic physical units
% (e.g., decibels).
%
%   LogDimVars handle units that use logarithmic scales, such as decibels (dB),
%   nepers, and similar units. These units represent ratios on a logarithmic
%   scale. LogDimVars should be used only for setting (through multiplication),
%   converting (through division), and displaying units, so DO NOT USE WITHOUT A
%   MULTIPLICATION OR DIVISION OPERATOR. The math in the background is always
%   operating on the underlying linearized value after the unit is set.
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
%     math is done on the underlying linear units, but in many cases the custom
%     log unit display will carry through the operation. This means that LOG
%     MATH DOESN'T HAPPEN. So, 20*u.dB + 10*u.dB does NOT yield 30 dB, but
%     rather yields 100 + 10 = 110 (~20.414 dB).
% 
%     IMPORTANT - Adding field quantities (root-power/amplitude units): Field
%     quantities like voltage, current, and pressure / SPL use a divisor of 20
%     (e.g., dBV, dBuV, dBSPL). When adding these quantities in code, the result
%     is LINEAR addition of the underlying physical quantities, NOT logarithmic
%     addition of the dB values:
%       - INCORRECT thinking: 20*u.dBV + 20*u.dBV = 40 dBV
%       - ACTUAL result: 10V + 10V = 20V (displays as ~26.02 dBV)
%     This behavior is correct for coherent (in-phase) signals. For addition of
%     incoherent/uncorrelated signals, first convert to power quantity or use
%     specialized functions like squaring the linear quantities. The same
%     applies to power quantities (divisor=10, e.g., dBW, dBm): adding them
%     performs linear addition of the underlying power, not the dB values.
% 
%     Examples: 
%       % Power quantities (divisor=10):
%       10*u.dBW + 10*u.dBW = 100 W + 100 W = 200 W (displays ~13 dBW)
%       
%       % Field quantities (divisor=20) - coherent (in-phase) addition:
%       10*u.dBV + 10*u.dBV = 10 V + 10 V = 20 V (displays ~16 dBV)
%       
%       % Field quantities - incoherent/uncorrelated addition:
%       sqrt((10*u.dBV)^2 + (10*u.dBV)^2) = ~14.14 V (displays ~13 dBV)
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
        reference = 1   % Reference value with units (DimVar or scalar).

        % Divisor for logarithm (10 for dB power, 20 root-power / field quantities).
        % See wikipedia.org/wiki/Power,_root-power,_and_field_quantities.
        divisor = 10    
    end
    

    methods
        function v = LogDimVar(reference,divisor)
            % LogDimVar  Construct a logarithmic unit
            %   v = LogDimVar(reference) creates a log unit with default divisor=10.
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
