function varargout = unitsOf(v)
% uv = UNITSOF(x) returns a value = 1 dimensioned variable with same units
% as x.
% 
%   [uv,us] = UNITSOF(x) returns the the same uv as with s single output
%   argument but in addition returns a string, us, describing the units of
%   x.
% 
%   Example: Use unitsOf to use MATLAB functions that are undefined for
%   dimensioned variables. (Note: MESHGRID is actually defined for
%   dimensioned variables - this is just an example.)
%       V = u.kph*(0:5:80); P = u.kW*(15:10:65);
%       tempV = u2num(V); tempP = u2num(P);
%       [tempV,tempP] = meshgrid(tempV,tempP);
%       V = tempV*unitsOf(V); P = tempP*unitsOf(P); 
% 
%   See also UNITS, U2NUM.

NDimensions = length(v.names);
nameString = '';
for nd = 1:NDimensions
    if(v.exponents(nd)~=0)
        [n,d] = rat(v.exponents(nd));
        if(d==1)
            nameString = sprintf...
                ('%s[%s^%g]',nameString,v.names{nd},v.exponents(nd));
        else
            nameString = ...
                sprintf('%s[%s^(%g/%g)]',nameString,v.names{nd},n,d);
        end
    end
end

unitsVar = v;
unitsVar.value = 1;

switch nargout % changed from if list to switch 2014-02-16
    case 0
        disp(nameString);
    case 1
        varargout = {unitsVar};
    case 2
        varargout = {unitsVar,nameString};
end