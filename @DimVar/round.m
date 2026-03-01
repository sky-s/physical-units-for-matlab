function v = round(v,varargin)
% ROUND Round DimVar based on its diplaying value.
%   Syntax is the same as built-in ROUND, with the exception of the special case
%   of specifying a DimVar instead of number of digits.
%
%   V = ROUND(V,INCREMENT) rounds V to the nearest multiple of INCREMENT
%   when INCREMENT is a DimVar. This is a special-case overload for DimVar
%   increments (for example, rounding 1.23 m to the nearest 0.05 m).
%
%   Note:
%   Using ROUND with DimVars may yield unexpected results depending on
%   display scaling.
%
%   Examples:
%       v = [1.23456, 12.3456, 123.456]*u.ft;
%       round(v)                          % Round to nearest foot.
%       round(scd(v,'m'), 2)              % Round to nearest cm, display in m.
%       round(v,u.in)                     % Round to nearest inch.
%       round(v,1,"significant")          % Round in feet, 1 significant digit.
%       round(v, 0.5*u.ft)                % Round to nearest half foot.
%       round(v, [1/12, 1, 5]*u.ft)       % Round to varying increments.
%
%   See also round.

% Handle special case where second argument is a DimVar (round to nearest increment)
if nargin > 1 && iscompatible(v,varargin{1})
    increment = varargin{1};
    v.value = round(v.value./increment.value).*increment.value;
    return
end
% warning('DimVar:round','Using round with DimVars may yield unexpected results.')

dispVal = displayparser(v);
k = round(dispVal,varargin{:})./dispVal;

v.value = v.value.*k;
