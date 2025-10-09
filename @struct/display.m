function display(s,varargin)
% Overload of display method to use dispdisp to show DimVars.
% 
% See also dispdisp.

try 
    dispdisp(s);
catch
    builtin('display',s,varargin{:});
end
