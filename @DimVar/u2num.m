function num = u2num(inVariable,inUnit)
% U2NUM converts a dimensioned variable to an undimensioned variable.
%
%   num = U2NUM(inVariable) converts a dimensioned variable to a normal
%   variable. This is essentially a typecast to a matlab numeric variable
%   and is very useful for forcing a dimensioned variable to work with
%   arbitrary matlab functions.
%
%   num = U2NUM(inVariable,inUnit) first converts inVariable to be in units
%   of inUnit. This is the same as num = inVariable/inUnit except that the
%   result is checked to make sure that the units of inVariable and inUnit
%   cancel.
%
%   Example:
%     t = (0:.01:2)*u.s; w = (2*pi)*2*u.Hz; y = (3*u.in)*sin(w.*t);
%     % plot(t,y); % BAD: plot function not defined for dimensioned input.
%     plot(u2num(t),u2num(y)); % Good
%     plot(t/u.s,y/u.cm); % Better: Explicitly state desired units to plot.
%     plot(u2num(t,u.s),u2num(y,u.cm)); % Also better.
%
%   See also U, UNITSOF, DISPLAYINGVALUE.

if nargin == 1
%     if isa(inVariable, 'DimVar')
%     % This "if isa" switch is for when this function is standalone and
%     % not a DimVar method.
        num = inVariable.value;
%     else
%         num = inVariable;
%     end
else
    num = inVariable/inUnit;
    if isa(num,'DimVar')
        error('Incompatible input dimensions')
    end
end
    
end