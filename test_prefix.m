% Test script to verify SI prefix functionality
% This script tests the working SI prefix implementation

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

% Test 2: Try accessing a prefixed unit using u.get()
fprintf('Test 2: Accessing u.get(''microinch'')\n');
try
    result2 = u.get('microinch');
    fprintf('  SUCCESS: u.get(''microinch'') = %s\n', string(result2));
catch ME
    fprintf('  ERROR: %s\n', ME.message);
end

% Test 3: Try the static method approach
fprintf('Test 3: Using static method u.microinch()\n');
try
    result3 = u.microinch();
    fprintf('  SUCCESS: u.microinch() = %s\n', string(result3));
catch ME
    fprintf('  ERROR: %s\n', ME.message);
end

% Test 4: Try with a different prefix
fprintf('Test 4: Accessing u.get(''kilojoule'')\n');
try
    result4 = u.get('kilojoule');
    fprintf('  SUCCESS: u.get(''kilojoule'') = %s\n', string(result4));
catch ME
    fprintf('  ERROR: %s\n', ME.message);
end

% Test 5: Demonstrate the MATLAB limitation
fprintf('Test 5: Attempting u.microinch (should fail due to MATLAB limitation)\n');
try
    result5 = u.microinch;  % This will fail
    fprintf('  UNEXPECTED SUCCESS: u.microinch = %s\n', string(result5));
catch ME
    fprintf('  EXPECTED ERROR: %s\n', ME.message);
end

fprintf('Test completed.\n');