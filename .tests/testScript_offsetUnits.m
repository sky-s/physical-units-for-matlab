% A test script for OffsetDimVar.
%   <a href="matlab:runtests('testScript_offsetUnits')">run tests</a>

% Set up by clearing class
clear all

%%   Set a value with multiplication (before or after scalar): 
assert(20*u.degC==293.15*u.K)
assert(u.degC*40== 313.15*u.K)

%% Set with str2u
assert(str2u('20 degC')==293.15*u.K)
assert(str2u('20 °C')==293.15*u.K)

%% Special str2u case
assert(isa(str2u('degC'),'OffsetDimVar'))

%%       Convert units with division:
assert(abs(200*u.K/u.degC - (-73.15)) < sqrt(eps))
assert(isequal(20*u.degC/u.K,293.15))
assert(abs(u.degC/u.degF - 1.8) < sqrt(eps))

%% Should errors
%     Pretty much all other use cases should be avoided and mostly throw an
%     error, but not always (usually in cases where it will be interpreted at 1
%     degC = 274.15 K, e.g.), so be careful.
shoulderror('DimVar:incompatibleUnits','u.kg*u.degC');
shoulderror('DimVar:incompatibleUnits','u.degC*u.kg');
shoulderror('OffsetDimVar:incompatibleUnits','u.degC*u.degF');
shoulderror('OffsetDimVar:undefined','u.degC/u.K');
shoulderror('OffsetDimVar:undefined',"str2u('20 degC/s')");
shoulderror('u.degC + 5*u.K');
shoulderror('DimVar:incompatibleUnits','unitconversionfactor(u.K,u.degC)');

%       anything else with unitconversionfactor (since it won't be a "factor"),
%       e.g. unitconversionfactor('degC','K').

shoulderror('DimVar:incompatibleUnits',"unitconversionfactor('degC','K')");

%% this syntax is bad practice, but it works
20*u.degC + 20*u.degF;
