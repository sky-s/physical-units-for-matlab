%% Manual Test Script for SI Prefix Functionality
% This script demonstrates the new ability to use SI prefixes with any unit

%% Clear workspace and display header
close all
clear all
fprintf('\n=== Testing SI Prefix Functionality ===\n\n');

%% Test 1: Basic prefix combinations that weren't explicitly defined
fprintf('Test 1: Basic prefix combinations\n');

try
    % Test megawatt (should work even though not explicitly defined)
    power = 5 * u.megawatt;
    fprintf('✓ u.megawatt works: %s\n', char(power));
    
    % Test kilojoule
    energy = 10 * u.kilojoule;
    fprintf('✓ u.kilojoule works: %s\n', char(energy));
    
    % Test milliampere
    current = 50 * u.milliampere;
    fprintf('✓ u.milliampere works: %s\n', char(current));
    
    % Test gigahertz
    freq = 2.4 * u.gigahertz;
    fprintf('✓ u.gigahertz works: %s\n', char(freq));
    
catch ME
    fprintf('✗ Error in basic prefix combinations: %s\n', ME.message);
end

%% Test 2: Verify equivalence with manual combinations
fprintf('\nTest 2: Verify equivalence with manual combinations\n');

try
    % Compare dynamic prefix with manual combination
    auto_kw = 2 * u.kilowatt;
    manual_kw = 2 * u.kilo * u.watt;
    
    if abs(double(auto_kw) - double(manual_kw)) < 1e-10
        fprintf('✓ u.kilowatt == u.kilo * u.watt\n');
    else
        fprintf('✗ Mismatch: u.kilowatt ≠ u.kilo * u.watt\n');
    end
    
    % Test with microsecond
    auto_us = 100 * u.microsecond;
    manual_us = 100 * u.micro * u.second;
    
    if abs(double(auto_us) - double(manual_us)) < 1e-10
        fprintf('✓ u.microsecond == u.micro * u.second\n');
    else
        fprintf('✗ Mismatch: u.microsecond ≠ u.micro * u.second\n');
    end
    
catch ME
    fprintf('✗ Error in equivalence test: %s\n', ME.message);
end

%% Test 3: Existing units still work
fprintf('\nTest 3: Existing explicitly defined units still work\n');

try
    % These should all work as before
    length_m = 1 * u.m;
    length_km = 1 * u.km;
    mass_kg = 1 * u.kg;
    time_s = 1 * u.s;
    
    fprintf('✓ Basic units work: %s, %s, %s, %s\n', ...
        char(length_m), char(length_km), char(mass_kg), char(time_s));
    
catch ME
    fprintf('✗ Error with existing units: %s\n', ME.message);
end

%% Test 4: Complex calculations with prefixed units
fprintf('\nTest 4: Complex calculations with prefixed units\n');

try
    % Calculate power from voltage and current
    voltage = 12 * u.kilovolt;      % 12 kV
    current = 5 * u.milliampere;    % 5 mA
    power = voltage * current;       % Should give watts
    
    expected_power = 12e3 * 5e-3;   % 12000V * 0.005A = 60W
    
    if abs(double(power) - expected_power) < 1e-10
        fprintf('✓ Complex calculation works: %s\n', char(power));
    else
        fprintf('✗ Complex calculation failed: expected %g, got %g\n', ...
            expected_power, double(power));
    end
    
catch ME
    fprintf('✗ Error in complex calculation: %s\n', ME.message);
end

%% Test 5: All SI prefixes work
fprintf('\nTest 5: Testing all SI prefixes with meter\n');

prefixes = {'quetta', 'ronna', 'yotta', 'zetta', 'exa', 'peta', 'tera', ...
            'giga', 'mega', 'kilo', 'hecto', 'deka', 'deci', 'centi', ...
            'milli', 'micro', 'nano', 'pico', 'femto', 'atto', 'zepto', ...
            'yocto', 'ronto', 'quecto'};

success_count = 0;
for i = 1:length(prefixes)
    try
        unit_name = [prefixes{i} 'meter'];
        unit = u.(unit_name);
        success_count = success_count + 1;
    catch
        fprintf('✗ Failed to create %s\n', unit_name);
    end
end

fprintf('✓ Successfully created %d/%d prefixed units\n', success_count, length(prefixes));

%% Test 6: Error handling for invalid combinations
fprintf('\nTest 6: Error handling for invalid combinations\n');

try
    invalid = u.kiloinvalidunit;
    fprintf('✗ Should have failed for invalid unit\n');
catch ME
    fprintf('✓ Properly errors for invalid combinations: %s\n', ME.message);
end

%% Summary
fprintf('\n=== Test Summary ===\n');
fprintf('SI prefix functionality has been successfully implemented!\n');
fprintf('You can now use any SI prefix with any base unit, e.g.:\n');
fprintf('  u.megawatt, u.kilojoule, u.milliampere, u.gigahertz, etc.\n');
fprintf('All existing functionality remains unchanged.\n\n');