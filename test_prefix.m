% Test script to verify SI prefix functionality
% This script tests if the subsref method is being called

% Clear any cached class information
clear classes

fprintf('Testing SI prefix functionality...\n');

% Test 1: Try accessing an existing unit
fprintf('Test 1: Accessing existing unit u.meter\n');
try
    result1 = u.meter;
    fprintf('  SUCCESS: u.meter = %s\n', string(result1));
catch ME
    fprintf('  ERROR: %s\n', ME.message);
end

% Test 2: Try accessing a non-existent unit that should trigger subsref
fprintf('Test 2: Accessing non-existent unit u.microinch\n');
try
    result2 = u.microinch;
    fprintf('  SUCCESS: u.microinch = %s\n', string(result2));
catch ME
    fprintf('  ERROR: %s\n', ME.message);
end

% Test 3: Try the static method approach
fprintf('Test 3: Using static method u.get(''microinch'')\n');
try
    result3 = u.get('microinch');
    fprintf('  SUCCESS: u.get(''microinch'') = %s\n', string(result3));
catch ME
    fprintf('  ERROR: %s\n', ME.message);
end

% Test 4: Try with a different prefix
fprintf('Test 4: Accessing u.kilojoule\n');
try
    result4 = u.kilojoule;
    fprintf('  SUCCESS: u.kilojoule = %s\n', string(result4));
catch ME
    fprintf('  ERROR: %s\n', ME.message);
end

fprintf('Test completed.\n');