% Test script to verify SI prefix functionality using new prefix methods
% This script tests the new SI prefix static method implementation

% Clear any cached class information
clear classes

fprintf('Testing SI prefix functionality with new static methods...\n');

% Test 1: Try accessing an existing unit
fprintf('Test 1: Accessing existing unit u.meter\n');
try
    result1 = u.meter;
    fprintf('  SUCCESS: u.meter = %s\n', string(result1));
catch ME
    fprintf('  ERROR: %s\n', ME.message);
end

% Test 2: Try accessing a prefix value (no arguments)
fprintf('Test 2: Getting prefix value u.micro()\n');
try
    result2 = u.micro();
    fprintf('  SUCCESS: u.micro() = %g\n', result2);
catch ME
    fprintf('  ERROR: %s\n', ME.message);
end

% Test 3: Try creating a prefixed unit with one argument
fprintf('Test 3: Creating prefixed unit u.micro(''inch'')\n');
try
    result3 = u.micro('inch');
    fprintf('  SUCCESS: u.micro(''inch'') = %s\n', string(result3));
catch ME
    fprintf('  ERROR: %s\n', ME.message);
end

% Test 4: Try another prefix with different unit
fprintf('Test 4: Creating u.kilo(''joule'')\n');
try
    result4 = u.kilo('joule');
    fprintf('  SUCCESS: u.kilo(''joule'') = %s\n', string(result4));
catch ME
    fprintf('  ERROR: %s\n', ME.message);
end

% Test 5: Try with full name option
fprintf('Test 5: Creating u.kilo(''joule'', true) for full name\n');
try
    result5 = u.kilo('joule', true);
    fprintf('  SUCCESS: u.kilo(''joule'', true) = %s\n', string(result5));
catch ME
    fprintf('  ERROR: %s\n', ME.message);
end

% Test 6: Try non-SI unit with prefix
fprintf('Test 6: Creating u.kilo(''acre'') for imperial unit\n');
try
    result6 = u.kilo('acre');
    fprintf('  SUCCESS: u.kilo(''acre'') = %s\n', string(result6));
catch ME
    fprintf('  ERROR: %s\n', ME.message);
end

% Test 7: Test multiple prefix methods
fprintf('Test 7: Testing various prefixes\n');
prefixes = {'nano', 'mega', 'giga', 'tera'};
for i = 1:length(prefixes)
    try
        prefix_val = u.(prefixes{i})();
        unit_val = u.(prefixes{i})('meter');
        fprintf('  SUCCESS: u.%s() = %g, u.%s(''meter'') = %s\n', ...
            prefixes{i}, prefix_val, prefixes{i}, string(unit_val));
    catch ME
        fprintf('  ERROR with %s: %s\n', prefixes{i}, ME.message);
    end
end

fprintf('Test completed.\n');