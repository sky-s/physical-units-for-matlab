% A test script for OffsetDimVar and LogDimVar.
%   <a href="matlab:runtests('testScript_offsetUnits')">run tests</a>

% Set up by clearing class
clear all

%% ========== OffsetDimVar Tests ==========

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

%% this syntax is bad practice, but it unfortunately works
20*u.degC + 20*u.degF;

%% ========== LogDimVar Tests ==========

%% Generic dB: Set a value with multiplication
% 20 dB power ratio = 10^(20/10) = 100
assert(abs(20*u.dB - 100) < sqrt(eps))
assert(abs(u.dB*30 - 1000) < sqrt(eps))

%% dBW: Power in decibels relative to 1 W
% 30 dBW = 10^(30/10) * 1W = 1000 W
assert(abs(30*u.dBW - 1000*u.W) < 1e-10*u.W)
assert(abs(0*u.dBW - 1*u.W) < 1e-10*u.W)

%% dBm: Power in decibels relative to 1 mW
% 30 dBm = 10^(30/10) * 1mW = 1000 mW = 1 W
assert(abs(30*u.dBm - 1*u.W) < 1e-10*u.W)
assert(abs(0*u.dBm - 1*u.mW) < 1e-13*u.W)

%% dBV: Voltage in decibels relative to 1 V (20*log10)
% 20 dBV = 10^(20/20) * 1V = 10 V
assert(abs(20*u.dBV - 10*u.V) < 1e-10*u.V)
assert(abs(0*u.dBV - 1*u.V) < 1e-10*u.V)

%% dBSPL: Sound pressure level (20 µPa reference)
% 94 dBSPL = 10^(94/20) * 20µPa ≈ 1 Pa
assert(abs(94*u.dBSPL - 1*u.Pa) < 0.01*u.Pa)

%% Convert from linear to dB with division
% 100 W / dBW should give 20 dB
powerInW = 100*u.W;
powerIndB = powerInW/u.dBW;
assert(abs(powerIndB - 20) < 1e-10)

% 1 W / dBm should give 30 dBm
assert(abs(1*u.W/u.dBm - 30) < 1e-10)

% 10 V / dBV should give 20 dBV
voltageInV = 10*u.V;
voltageIndB = voltageInV/u.dBV;
assert(abs(voltageIndB - 20) < 1e-10)

%% Round trip: dB to linear to dB
dbValue = 25;
linear = dbValue*u.dB;
backToDb = linear/u.dB;
assert(abs(backToDb - dbValue) < 1e-10)

%% Addition of dB values is addition of underlying linear quantities (same units)
% 10 dBW + 10 dBW = 100 W + 100 W = 200 W
% Note: This is not 20 dBW, which would be the case if adding ratios
result = 10*u.dBW + 10*u.dBW;
assert(abs(result - 20*u.W) < 1e-10*u.W) 

%% Subtraction of dB values 
assert(abs(30*u.dBm - 30*u.dBm) < eps*u.W)

%% LogDimVar error cases
% Cannot multiply two LogDimVars
shoulderror('LogDimVar:incompatibleUnits','u.dB*u.dB');
shoulderror('LogDimVar:incompatibleUnits','u.dBm*u.dBW');

% Cannot divide a LogDimVar by a scalar
shoulderror('LogDimVar:undefined','u.dBW/2');

% Cannot add different LogDimVar types
shoulderror('u.dBW + u.dBV');
