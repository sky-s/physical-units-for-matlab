function [tf,ME] = iscompatible(varargin)
% Returns true if all inputs are DimVars with the same units (compatible
% executes successfully) and otherwise returns false.
% 
%   See also u, compatible.

try 
    compatible(varargin{:});
    tf = true;
catch ME
    % Capture all-double case.
    if all(cellfun('isclass',varargin,'double'))
        tf = true;
        return
    end

    if strcmp(ME.identifier,'DimVar:incompatibleUnits')
        tf = false;
    else
        rethrow(ME)
    end

end