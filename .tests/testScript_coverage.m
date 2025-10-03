% Additional test cases to improve code coverage
%   <a href="matlab:runtests('testScript_coverage')">run tests</a>

% Set up by clearing class
clear all
close all

%% Test scale function
% Test with DimVar
x = 5*u.m;
scaled_x = scale(x, 1000);
assert(u2num(scaled_x) == 5000);
assert(unitsOf(scaled_x) == u.m);

% Test with numeric
y = [1 2 3];
scaled_y = scale(y, 10);
assert(isequal(scaled_y, y));

%% Test displayingvalue function
% Test basic functionality
dv = displayingvalue(100*u.cm);
assert(isequal(dv, 100)); 

% Test with custom display units
set_length = scd(150*u.cm, 'ft');
dv2 = displayingvalue(set_length);
% Should be in feet
assert(abs(u2num(dv2) - 150*u.cm/u.ft) < 1e-10);

%% Test plottingvalue function  
% Test basic conversion
pv = plottingvalue(1000*u.mm);
assert(pv == 1000);

pv2 = plottingvalue(scd(3*u.ft));
expected = 3*u.ft/u.m;
assert(abs((pv2) - (expected)) < 1e-10);

%% Test compatible function with multiple inputs
% Should pass without error
compatible(u.m, u.ft, u.km, u.in);
compatible(u.kg, u.lb, u.g);

% Single input should pass
compatible(u.Pa);

%% Test iscompatible function edge cases
% Test with empty inputs
assert(iscompatible([]));

% Test with mixed numeric and DimVar
assert(~iscompatible(5, u.m));
assert(~iscompatible(u.m, 'string'));

% Test with OffsetDimVar
assert(iscompatible(5*u.degC, u.K));
assert(~iscompatible(u.degC, u.m));

%% Test unitsOf function edge cases
% Test with scalar
assert(unitsOf(42) == 1);

% Test with empty arrays
assert(unitsOf([]) == 1);

% Test with complex numbers
assert(unitsOf(3+4i) == 1);
assert(unitsOf((3+4i)*u.m) == u.m);


%% Test baseUnitSystem function
% Test setting and getting
old_base = baseUnitSystem();
baseUnitSystem('SI');

% Test 'none' setting
baseUnitSystem('FPS');

% Restore original
baseUnitSystem('SI');

% %% Test displayUnits function
% % Test setting and getting
% old_display = displayUnits();
% displayUnits({'ft', 'lb', 'deg'});
% current = displayUnits();
% assert(any(strcmp(current, 'ft')));
% assert(any(strcmp(current, 'lb')));
% 
% % Test clearing
% displayUnits({});
% assert(isempty(displayUnits()));
% 
% % Restore original
% displayUnits(old_display);

%% Test DimVar constructor edge cases
% Test with zero
zero_mass = DimVar([1 0 0 0 0 0 0 0 0 0],0);
assert(u2num(zero_mass) == 0);

% Test with negative values  
neg_temp = DimVar([0 0 0 0 1 0 0 0 0 0], -10);
assert(u2num(neg_temp) == -10);


%% Test subsref edge cases with DimVar
test_array = u.m*[1 2; 3 4];

% Test linear indexing
assert(test_array(3) == 2*u.m);

% Test end keyword
assert(test_array(end) == 4*u.m);
assert(test_array(1,end) == 2*u.m);

% Test colon operator
assert(isequal(test_array(:), u.m*[1; 3; 2; 4]));
assert(isequal(test_array(1,:), u.m*[1 2]));

%% Test arithmetic edge cases
% Test with inf and nan
inf_length = inf*u.m;
assert(isinf(inf_length));
assert(unitsOf(inf_length) == u.m);

nan_mass = nan*u.kg;
assert(isnan(nan_mass));
assert(unitsOf(nan_mass) == u.kg);

% Test mixed operations
result = inf_length + 5*u.m;
assert(isinf(result));

%% Test string representation functions
% Test mat2str
test_val = 3.14159*u.m;
str_rep = mat2str(test_val);
assert(ischar(str_rep));
assert(contains(str_rep, '3.14159'));

% Test num2str  
num_str = num2str(test_val);
assert(ischar(num_str));

%% Test logical operations edge cases
% Test with zero
assert(~(0*u.m));
assert(logical(1*u.kg));

% Test comparisons with mixed signs
assert(-5*u.m < 3*u.m);
assert(abs(-5*u.m) > abs(3*u.m));

%% Test reduction operations
test_matrix = u.Pa*magic(3);

% Test sum
col_sum = sum(test_matrix);
assert(isequal(size(col_sum), [1 3]));
assert(unitsOf(col_sum) == u.Pa);

row_sum = sum(test_matrix, 2);
assert(isequal(size(row_sum), [3 1]));

% Test mean
mean_val = mean(test_matrix(:));
assert(isscalar(mean_val));
assert(unitsOf(mean_val) == u.Pa);

% Test std
std_val = std(test_matrix(:));
assert(isscalar(std_val));
assert(unitsOf(std_val) == u.Pa);

%% Test sorting operations
unsorted = u.kg*[3 1 4 1 5];
[sorted_vals, indices] = sort(unsorted);
assert(isequal(u2num(sorted_vals), [1 1 3 4 5]));
assert(unitsOf(sorted_vals) == u.kg);

%% Test min/max operations
test_data = u.W*[10 5 20 15];
[min_val, min_idx] = min(test_data);
assert(min_val == 5*u.W);
assert(min_idx == 2);

[max_val, max_idx] = max(test_data);
assert(max_val == 20*u.W);
assert(max_idx == 3);

%% Test reshape operations
original = u.m*[1 2 3 4 5 6];
reshaped = reshape(original, [2 3]);
assert(isequal(size(reshaped), [2 3]));
assert(unitsOf(reshaped) == u.m);
assert(reshaped(1,1) == 1*u.m);
assert(reshaped(2,3) == 6*u.m);

%% Test permute operations
test_3d = u.K*rand(2,3,4);
permuted = permute(test_3d, [3 1 2]);
assert(isequal(size(permuted), [4 2 3]));
assert(unitsOf(permuted) == u.K);

%% Test transpose operations  
test_matrix = u.N*[1 2; 3 4];
transposed = test_matrix';
assert(isequal(size(transposed), [2 2]));
assert(transposed(1,2) == 3*u.N);
assert(transposed(2,1) == 2*u.N);

% Test ctranspose with complex
complex_matrix = (1+1i)*u.V*[1 2; 3 4];
ctrans = complex_matrix';
assert(ctrans(1,1) == (1-1i)*u.V);

%% Test concatenation edge cases
% Test with empty arrays
a = [1 2]*u.m;
b = []*u.m;
c = [3 4]*u.m;
result = [a b c];
assert(isequal(result, [1 2 3 4]*u.m));

% Test vertical concatenation
vert_result = [a; c];
assert(isequal(size(vert_result), [2 2]));

%% Test fillmissing function
data_with_missing = u.Pa*[1 NaN 3 NaN 5];
filled = fillmissing(data_with_missing, 'linear');
assert(~any(isnan(filled)));
assert(unitsOf(filled) == u.Pa);

%% Test diff function  
sequence = u.m*[1 4 9 16 25];  % squares
differences = diff(sequence);
expected_diff = u.m*[3 5 7 9];  % differences of squares
assert(isequal(differences, expected_diff));

%% Test cumulative operations
data = u.J*[1 2 3 4];
cumsum_result = cumsum(data);
assert(isequal(cumsum_result, u.J*[1 3 6 10]));
assert(unitsOf(cumsum_result) == u.J);

%% Cleanup
clear all
close all
