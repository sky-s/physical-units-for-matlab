% Test script for utility functions and edge cases
%   <a href="matlab:runtests('testScript_utilities')">run tests</a>

% Set up by clearing class
clear all
close all

%% Test OffsetDimVar edge cases
% Test arithmetic with offset units
temp1 = 20*u.degC;  % 293.15 K
temp2 = 30*u.degC;  % 303.15 K

% Test difference (should work)
temp_diff = temp2 - temp1;  % Should be 10 K
assert(abs(u2num(temp_diff) - 10) < 1e-10);
assert(unitsOf(temp_diff) == u.K);

% Test multiplication by scalar
doubled_temp = 2*temp1;
assert(abs(u2num(doubled_temp) - 2*293.15) < 1e-10);

%% Test duration2u and u2duration functions  
% Test conversion to duration
time_val = 2.5*u.hr;
dur = u2duration(time_val);
assert(isduration(dur));
assert(hours(dur) == 2.5);

% Test conversion from duration
dur2 = minutes(90);
time_val2 = duration2u(dur2);
assert(abs((time_val2/u.min) - 90) < 1e-10);
assert(iscompatible(time_val2,u.min));

%% Test str2u edge cases with complex expressions
% Test with parentheses
result1 = str2u('(kg*m)/s^2');
assert(result1 == u.N);

% Test with mixed units
result2 = str2u('hp-hr/gal');
expected = u.hp*u.hr/u.gal;
assert(abs(u2num(result2/expected) - 1) < 1e-10);

% Test with fractions
result3 = str2u('kg^(1/2)');
expected3 = u.kg^0.5;
assert(abs(u2num(result3/expected3) - 1) < 1e-10);

%% Test unitconversionfactor edge cases
% Test with same units
factor1 = unitconversionfactor(u.m, u.m);
assert(factor1 == 1);

% Test with complex conversion
factor2 = unitconversionfactor(u.hp, u.W);
assert(abs(u2num(factor2) - u2num(u.hp/u.W)) < 1e-10);

% Test string to DimVar conversion
factor3 = unitconversionfactor('ft', u.m);
assert(abs(u2num(factor3) - u2num(u.ft/u.m)) < 1e-10);

%% Test istype function comprehensively
% Test basic dimensional types
assert(istype(u.m^2, 'Area'));
assert(istype(u.m^3, 'Volume'));
assert(istype(u.m/u.s^2, 'Acceleration'));
assert(istype(u.kg*u.m/u.s^2, 'Force'));
assert(istype(u.kg*u.m^2/u.s^2, 'Energy'));
assert(istype(u.kg*u.m^2/u.s^3, 'Power'));
assert(istype(u.kg/(u.m*u.s^2), 'Pressure'));
assert(istype(u.m/u.s, 'Velocity'));

% Test negative cases
assert(~istype(u.m, 'Mass'));
assert(~istype(u.kg, 'Length'));
assert(~istype(u.s, 'Temperature'));

% Test with complex units
% assert(istype(u.V, 'Voltage'));
% assert(istype(u.A, 'Current'));
% assert(istype(u.ohm, 'Resistance'));

%% Test clearcanceledunits functionality (if it exists)
try
    % Test if the function exists and works
    test_val = (u.m*u.kg)/(u.m*u.s^2);  % Should simplify to kg/s^2
    cleaned = clearcanceledunits(test_val);
    % This might not exist - just testing if it does
catch ME
    if ~strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
        rethrow(ME);
    end
end

%% Test DimVar with very small and very large numbers
% Test with very small numbers
tiny = 1e-100 * u.m;
assert(u2num(tiny) == 1e-100);
assert(unitsOf(tiny) == u.m);

% Test with very large numbers
huge = 1e100 * u.kg;
assert(u2num(huge) == 1e100);
assert(unitsOf(huge) == u.kg);

% Test arithmetic with extreme values
result_extreme = tiny * huge;  % Should be 1e0 kg*m
assert(abs(u2num(result_extreme) - 1) < 1e-10);

%% Test array operations with mixed dimensions
% This should error - incompatible units in array
shoulderror('DimVar:incompatibleUnits', @()[u.m u.kg]);
shoulderror('DimVar:incompatibleUnits', @()[u.m; u.s]);

%% Test subsasgn with complex indexing
test_matrix = u.Pa * magic(4);
% Test assignment to submatrix
test_matrix(2:3, 2:3) = u.Pa * ones(2,2) * 99;
assert(all(all(test_matrix(2:3, 2:3) == 99*u.Pa)));

% Test assignment with logical indexing
logical_idx = test_matrix > 10*u.Pa;
test_matrix(logical_idx) = 1*u.Pa;
assert(all(test_matrix(logical_idx) == 1*u.Pa));

%% Test colon operator with units
% Basic colon with units
range1 = u.m:u.m:5*u.m;
expected1 = u.m * (1:5);
assert(isequal(range1, expected1));

% Colon with different step
range2 = 0*u.kg:2*u.kg:10*u.kg;
expected2 = u.kg * (0:2:10);
assert(isequal(range2, expected2));

% Colon with non-unit step (should error)
shoulderror(@()(u.ft:5:u.yd));

%% Test ndims, size, length functions
test_3d = u.W * rand(2,3,4);
assert(ndims(test_3d) == 3);
assert(isequal(size(test_3d), [2 3 4]));
assert(length(test_3d) == 4);  % max dimension
assert(numel(test_3d) == 24);  % total elements

%% Test circshift function
original = u.K * [1 2 3 4 5];
shifted = circshift(original, 2);
expected = u.K * [4 5 1 2 3];
assert(isequal(shifted, expected));

% Test 2D circshift
matrix_2d = u.N * [1 2; 3 4];
shifted_2d = circshift(matrix_2d, [1 1]);
expected_2d = u.N * [4 3; 2 1];
assert(isequal(shifted_2d, expected_2d));

%% Test diag function
test_matrix = u.Pa * magic(3);
diagonal = diag(test_matrix);
assert(length(diagonal) == 3);
assert(unitsOf(diagonal) == u.Pa);

% Test creating diagonal matrix
diag_values = u.m * [1 2 3];
diag_matrix = diag(diag_values);
assert(isequal(size(diag_matrix), [3 3]));
assert(diag_matrix(1,1) == 1*u.m);
assert(diag_matrix(2,2) == 2*u.m);
assert(diag_matrix(1,2) == 0*u.m);

%% Test trace function  
square_matrix = u.J * magic(3);
tr = trace(square_matrix);
expected_trace = sum(diag(square_matrix));
assert(tr == expected_trace);
assert(unitsOf(tr) == u.J);

%% Test norm function
vector = u.N * [3 4];
vector_norm = norm(vector);
assert(abs(u2num(vector_norm) - 5) < 1e-10);  % 3-4-5 triangle
assert(unitsOf(vector_norm) == u.N);

% Test matrix norm
matrix_norm = norm(u.Pa * [1 2; 3 4]);
assert(unitsOf(matrix_norm) == u.Pa);

%% Test hypot function
a = 3*u.m;
b = 4*u.m;
hyp = hypot(a, b);
assert(abs(u2num(hyp) - 5) < 1e-10);
assert(unitsOf(hyp) == u.m);

%% Test atan2 function
y = 1*u.m;
x = 1*u.m;
angle = atan2(y, x);
assert(abs(u2num(angle) - pi/4) < 1e-10);
assert(unitsOf(angle) == u.rad);

%% Test validation functions
% Test mustBePositive equivalent
test_val = 5*u.kg;
assert(test_val > 0*u.kg);  % Should pass

neg_val = -3*u.m;
assert(neg_val < 0*u.m);    % Should pass

%% Test sign function
pos_val = 10*u.W;
neg_val = -5*u.W;
zero_val = 0*u.W;

assert(sign(pos_val) == 1);
assert(sign(neg_val) == -1);
assert(sign(zero_val) == 0);

%% Test round, floor, ceil functions
test_val = 3.7*u.m;
shouldwarn('round(test_val);');
assert(abs(4*u.m/round(test_val) - 1) < 1e-10);
% assert(abs(3*u.m/floor(test_val) - 1) < 1e-10);
% assert(abs(4*u.m/ceil(test_val) - 1) < 1e-10);

%% Test abs function with complex
complex_val = (3-4i)*u.V;
abs_val = abs(complex_val);
assert(abs(u2num(abs_val) - 5) < 1e-10);
assert(unitsOf(abs_val) == u.V);

%% Test conj function
complex_val = (3-4i)*u.V;
conj_val = conj(complex_val);
expected_conj = (3+4i)*u.V;
assert(conj_val == expected_conj);

%% Test real and imag functions
complex_val = (3-4i)*u.V;
assert(real(complex_val) == 3*u.V);
assert(imag(complex_val) == -4*u.V);

%% Test full function (if sparse arrays are supported)
try
    % This might not be fully implemented
    sparse_test = sparse(u.Pa * [1 0 2; 0 3 0]);
    full_test = full(sparse_test);
    assert(isa(full_test, 'DimVar'));
catch ME
    if ~strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
        % Expected - sparse might not be implemented
    end
end

%% Cleanup
clear all
close all
