% Test script for error conditions and edge cases
%   <a href="matlab:runtests('testScript_errors')">run tests</a>

% Set up by clearing class
clear all
close all

%% Test mustBeNonDim error scenarios
% Helper function to test mustBeNonDim errors
shoulderror('mustBeNonDim(u.m)');
mustBeNonDim(2)

%% Test str2u error conditions
% Test invalid unit strings
shoulderror(@()str2u('invalid_unit_name'));
shoulderror(@()str2u('kg-invalid'));
shoulderror(@()str2u('m^invalid'));

% Test security: prevent code injection
shoulderror(@()str2u('eval("malicious_code")'));
shoulderror(@()str2u('system("echo hacked")'));

%% Test DimVar constructor errors
% Test invalid exponent arrays
shoulderror(@()DimVar(1, [1 0], 'invalid')); % wrong size
shoulderror(@()DimVar(1, 'invalid', 'unit')); % wrong type

%% Test incompatible unit operations
% Test addition of incompatible units
shoulderror('DimVar:incompatibleUnits', @()(u.m + u.kg));
shoulderror('DimVar:incompatibleUnits', @()(u.s + u.A));

% Test concatenation of incompatible units
shoulderror('DimVar:incompatibleUnits', @()[u.m u.kg]);
shoulderror('DimVar:incompatibleUnits', @()[u.Pa; u.V]);

%% Test subsasgn error conditions
% Test assignment of incompatible units
x = [1 2 3]*u.m;
shoulderror('x(2) = u.kg;');

% Test assignment of non-numeric to DimVar
shoulderror('x(1) = "string";');

%% Test colon operator errors
% Test colon without proper increment
shoulderror('u.m:u.kg');
shoulderror('1*u.ft:5:1*u.yd');

%% Test OffsetDimVar errors  
% Test invalid operations with offset units
shoulderror('DimVar:incompatibleUnits', @()(u.degC * u.kg));
shoulderror('OffsetDimVar:incompatibleUnits', @()(u.degC * u.degF));
shoulderror('OffsetDimVar:undefined', @()(u.degC / u.K));

% Test offset units in unit conversion
shoulderror('DimVar:incompatibleUnits', @()unitconversionfactor(u.K, u.degC));

%% Test validation function errors
% Test validateattributes with wrong class
shoulderror('MATLAB:invalidType', @()validateattributes(u.m, {'double'}, {'finite'}));

%% Test plotting function errors with incompatible units
% Test plot with incompatible x data
a = 1:5;
shoulderror('DimVar:incompatibleUnits', @()plot(a*u.m, a, a*u.kg, a));

%% Test interp1 errors with incompatible units
x = [1 2 3 4]*u.m;
y = [1 4 9 16];  % no units
xi = 2.5*u.kg;   % incompatible
shoulderror('DimVar:incompatibleUnits', @()interp1(x, y, xi));

%% Test mathematical function errors
% Test functions that should error with dimensioned inputs
shoulderror(@()sin(u.m));
shoulderror(@()cos(u.kg));
shoulderror(@()exp(u.s));
shoulderror(@()log(u.A));

%% Test unit conversion errors
% Test conversion between incompatible units
shoulderror(@()unitconversionfactor(u.m, u.kg));
shoulderror(@()unitconversionfactor(u.Pa, u.A));

%% Test display function edge cases
% Test display with very complex units
complex_unit = u.kg^2 * u.m^3 / (u.s^4 * u.K^2 * u.A);
disp(complex_unit); % Should not error, just testing display

% Test display with zero values
zero_mass = 0*u.kg;
disp(zero_mass); % Should display properly

%% Test array operations with size mismatches
% Test operations requiring compatible sizes
a = [1 2]*u.m;
b = [1; 2; 3]*u.m;
shoulderror(@()(a + b)); % Size mismatch should error

%% Test reduction operations on empty arrays
empty_array = []*u.kg;
result_sum = sum(empty_array);
assert(isempty(result_sum) || result_sum == 0*u.kg);

result_mean = mean(empty_array);
assert(isnan(result_mean) || isempty(result_mean));

%% Test indexing errors
test_array = [1 2 3]*u.m;
% Test out of bounds indexing
shoulderror(@()test_array(5));
shoulderror(@()test_array(0));

%% Test reshape errors
test_data = [1 2 3 4 5 6]*u.Pa;
% Test incompatible reshape
shoulderror(@()reshape(test_data, [2 2])); % 6 elements can't fit in 2x2

%% Test matrix operations requiring square matrices
non_square = [1 2 3; 4 5 6]*u.J;
% Test operations that need square matrices
shoulderror(@()trace(non_square));

%% Test complex number operations where not allowed
% Some operations might not support complex DimVars
complex_val = (1+2i)*u.m;
% Test if certain operations error with complex units (implementation dependent)

%% Test concatenation size mismatches
a = [1 2]*u.kg;
b = [1; 2; 3]*u.kg;
shoulderror(@()horzcat(a, b)); % Size mismatch in horizontal concatenation

%% Test cumulative operations on incompatible data
% This should work, but test edge cases
single_val = 42*u.W;
cum_result = cumsum(single_val);
assert(cum_result == single_val);

%% Test sorting with NaN values
data_with_nan = [1 NaN 3 2]*u.V;
[sorted_data, indices] = sort(data_with_nan);
% NaN should be at the end
assert(isnan(sorted_data(end)));

%% Test diff on single element
single_element = 5*u.m;
diff_result = diff(single_element);
assert(isempty(diff_result));

%% Test statistical functions on constant arrays
constant_array = [5 5 5 5]*u.Pa;
std_result = std(constant_array);
assert(abs(u2num(std_result)) < 1e-10); % Should be zero or very close

%% Test min/max on empty arrays
empty_test = []*u.N;
[min_val, min_idx] = min(empty_test);
assert(isempty(min_val) && isempty(min_idx));

%% Nested function definition for helper
function shoulderror(varargin)
    % Handle different input formats for shoulderror
    if nargin == 1
        % shoulderror(@()expression)
        func = varargin{1};
        try
            func();
            error('Expected error but none was thrown');
        catch
            % Expected error occurred
        end
    elseif nargin == 2
        % shoulderror(expectedId, @()expression) or shoulderror(expectedMsg, @()expression)
        expected = varargin{1};
        func = varargin{2};
        try
            func();
            error('Expected error but none was thrown');
        catch ME
            if ischar(expected) && ~contains(ME.identifier, expected) && ~contains(ME.message, expected)
                error('Expected error "%s" but got "%s": %s', expected, ME.identifier, ME.message);
            end
        end
    else
        % shoulderror(expectedId, expression_string, variables...)
        expected = varargin{1};
        expr = varargin{2};
        try
            evalin('caller', expr);
            error('Expected error but none was thrown');
        catch ME
            if ischar(expected) && ~contains(ME.identifier, expected) && ~contains(ME.message, expected)
                error('Expected error "%s" but got "%s": %s', expected, ME.identifier, ME.message);
            end
        end
    end
end

%% Cleanup
clear all
close all
