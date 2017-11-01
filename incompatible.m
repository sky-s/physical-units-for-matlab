function tf = incompatible(varargin)
% Returns true if all inputs are not DimVars with the same units (if compatible
% throws an error) and otherwise returns false.
% 
%   See also u, compatible.

try
    tf = ~compatible(varargin{:}); 
    % Return false if compatible.
catch ME
    if strcmp(ME.identifier,'DimVar:incompatibleUnits')
        tf = true;
    else
        rethrow(ME);
    end
end