% A quick test of baseUnitSystem = 'none'.
%   <a href="matlab:runtests('testScript_noBaseUnits')">run tests</a>

% Set up by clearing class
clear all

baseUnitSystem = 'none';

%% Normal
assert(1==u.m)  %FIXME: fails with `Incompatible units. Cannot perform operation on variables with different units.` what is an expected behaviour of comparing a number with unit? magnitude comparison?
assert(0.3048==u.ft)%FIXME: fails with `Incompatible units. Cannot perform operation on variables with different units.` what is an expected behaviour of comparing a number with unit? magnitude comparison?

%% Offset
assert(isa(u.degC,'OffsetDimVar'))
assert(278.15==5*u.degC)%FIXME: fails with `Incompatible units. Cannot perform operation on variables with different units.` what is an expected behaviour of comparing a number with unit? magnitude comparison?
assert((5*u.degC)/u.degC==5)

%% Cleanup
clear all