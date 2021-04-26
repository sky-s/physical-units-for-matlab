% A quick test of baseUnitSystem = 'none'.
%   <a href="matlab:runtests('testScript_noBaseUnits')">run tests</a>

% Set up by clearing class
clear all

baseUnitSystem = 'none';

%% Normal
assert(1==u.m)
assert(0.3048==u.ft)

%% Offset
assert(isa(u.degC,'OffsetDimVar'))
assert(278.15==5*u.degC)
assert((5*u.degC)/u.degC==5)

%% Cleanup
clear all