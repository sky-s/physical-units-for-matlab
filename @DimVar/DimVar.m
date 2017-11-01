classdef DimVar
    % See also u, myUnits, units.
    properties (Access = protected)
        names
        exponents
        value = 1
        dispUnits
    end
    
    methods
        function v = DimVar(dimensionNames,dimensionToCreate,dispUnits)
            % See also u.
            v.names = dimensionNames;
            v.exponents = zeros(size(v.names));
            dimensionIndex = strcmp(v.names,dimensionToCreate);
            v.exponents(dimensionIndex) = 1;
            if nargin == 3
                v.dispUnits = dispUnits;
            end
        end
        
        function validateattributes(a,classes,varargin)
            % Check to make sure that DimVar was specified as an okay input.
            if any(strcmp('DimVar',classes))
                % DimVar okay. Value should be double, so allow it.
                classes = strrep(classes,'DimVar','double');
                validateattributes(a.value,classes,varargin{:})
            else
                builtin('validateattributes',a,classes,varargin{:})
            end
        end
    end
    methods (Static)
        [listString,list] = unitslist(varargin)
        [cTo,cInverse] = strconv(sFrom,sTo,varargin)
        
        function displayunitsstruct(u)
            % Display the struct created by calling units.
            % 
            %   See also units.
            if nargin < 1
                u = units;
            end
            uNames = fieldnames(u);
            nFields = length(uNames);
            lastUnits = 'wR!2V!17r@OenDg6wH&t';
            for nf = 1:nFields
                name = uNames{nf};
                [~,currentUnits] = unitsOf(u.(name));
                if(~strcmp(currentUnits,lastUnits))
                    fprintf('\n');
                    lastUnits = currentUnits;
                end
                
                fprintf('%16s: %-s\n',name,num2str(u.(name)));
            end
        end
    end
end