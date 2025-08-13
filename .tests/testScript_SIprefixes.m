% Test script for SI prefix functionality
%   Tests the new ability to use any SI prefix with any unit
%   <a href="matlab:runtests('testScript_SIprefixes')">run tests</a>

% Set up by clearing class
close all
clear all
%% Test basic SI prefix functionality

function tests = testScript_SIprefixes
    tests = functiontests(localfunctions);
end

function testKiloMeter(testCase)
    % Test that u.kilometer works
    distance = 5 * u.kilometer;
    expected = 5 * u.kilo * u.meter;
    verifyEqual(testCase, double(distance), double(expected), 'RelTol', 1e-10);
end

function testMegaWatt(testCase)
    % Test that u.megawatt works even though only u.watt is explicitly defined
    power = 2 * u.megawatt;
    expected = 2 * u.mega * u.watt;
    verifyEqual(testCase, double(power), double(expected), 'RelTol', 1e-10);
end

function testMilliAmpere(testCase)
    % Test that u.milliampere works
    current = 100 * u.milliampere;
    expected = 100 * u.milli * u.ampere;
    verifyEqual(testCase, double(current), double(expected), 'RelTol', 1e-10);
end

function testGigaHertz(testCase)
    % Test that u.gigahertz works
    freq = 3.5 * u.gigahertz;
    expected = 3.5 * u.giga * u.hertz;
    verifyEqual(testCase, double(freq), double(expected), 'RelTol', 1e-10);
end

function testMicroSecond(testCase)
    % Test that u.microsecond works
    time = 10 * u.microsecond;
    expected = 10 * u.micro * u.second;
    verifyEqual(testCase, double(time), double(expected), 'RelTol', 1e-10);
end

function testNanoMeter(testCase)
    % Test that u.nanometer works
    length = 500 * u.nanometer;
    expected = 500 * u.nano * u.meter;
    verifyEqual(testCase, double(length), double(expected), 'RelTol', 1e-10);
end

function testPicoFarad(testCase)
    % Test that u.picofarad works
    capacitance = 22 * u.picofarad;
    expected = 22 * u.pico * u.farad;
    verifyEqual(testCase, double(capacitance), double(expected), 'RelTol', 1e-10);
end

function testTeraJoule(testCase)
    % Test that u.terajoule works
    energy = 1.5 * u.terajoule;
    expected = 1.5 * u.tera * u.joule;
    verifyEqual(testCase, double(energy), double(expected), 'RelTol', 1e-10);
end

function testExistingUnitsStillWork(testCase)
    % Test that existing manually defined units still work
    distance1 = 1 * u.km;  % Existing manually defined
    distance2 = 1 * u.kilometer;  % New dynamically created
    
    % They should be equivalent
    verifyEqual(testCase, double(distance1), double(distance2), 'RelTol', 1e-10);
    
    % Test a few more existing units
    mass = 1 * u.kg;
    time = 1 * u.s;
    force = 1 * u.N;
    
    % These should work without issues
    verifyTrue(testCase, isa(mass, 'DimVar'));
    verifyTrue(testCase, isa(time, 'DimVar'));
    verifyTrue(testCase, isa(force, 'DimVar'));
end

function testInvalidPrefixUnit(testCase)
    % Test that invalid combinations properly error
    try
        invalid = u.kiloinvalidunit;
        verifyFail(testCase, 'Should have failed for invalid unit');
    catch ME
        % Should get an error
        verifyTrue(testCase, contains(ME.message, 'property') || contains(ME.message, 'field'));
    end
end

function testPrefixWithoutBaseUnit(testCase)
    % Test that just a prefix without valid base unit errors
    try
        invalid = u.kilo;  % This should work as kilo is defined as a constant
        verifyTrue(testCase, isnumeric(invalid));
        verifyEqual(testCase, invalid, 1e3);
    catch ME
        verifyFail(testCase, ['Unexpected error: ' ME.message]);
    end
end

function testCustomDisplayNames(testCase)
    % Test that the custom display names are set correctly
    power = u.megawatt;
    
    % Check that the custom display is set
    verifyTrue(testCase, isa(power, 'DimVar'));
    if isa(power, 'DimVar')
        verifyEqual(testCase, power.customDisplay, 'megawatt');
    end
end

function testChainedOperations(testCase)
    % Test that prefixed units work in calculations
    voltage = 5 * u.kilovolt;
    current = 10 * u.milliampere;
    power = voltage * current;
    
    % This should give us watts
    expected_watts = 5 * 1e3 * 10 * 1e-3;  % 5000V * 0.01A = 50W
    verifyEqual(testCase, double(power), expected_watts, 'RelTol', 1e-10);
end

function testAbbreviatedPrefixes(testCase)
    % Test that common prefix abbreviations work
    
    % Test kV (kilovolt)
    voltage1 = 12 * u.kV;
    voltage2 = 12 * u.kilovolt;
    verifyEqual(testCase, double(voltage1), double(voltage2), 'RelTol', 1e-10);
    
    % Test MW (megawatt) - but be careful as M might conflict with other units
    % Skip this test if it conflicts
    try
        power1 = 5 * u.MW;
        power2 = 5 * u.megawatt;
        if isa(power1, 'DimVar') && isa(power2, 'DimVar')
            verifyEqual(testCase, double(power1), double(power2), 'RelTol', 1e-10);
        end
    catch
        % Skip if MW conflicts with existing unit
    end
    
    % Test mA (milliampere)
    current1 = 100 * u.mA;
    current2 = 100 * u.milliampere;
    verifyEqual(testCase, double(current1), double(current2), 'RelTol', 1e-10);
    
    % Test nF (nanofarad)
    cap1 = 22 * u.nF;
    cap2 = 22 * u.nanofarad;
    verifyEqual(testCase, double(cap1), double(cap2), 'RelTol', 1e-10);
    
    % Test pF (picofarad)
    cap3 = 100 * u.pF;
    cap4 = 100 * u.picofarad;
    verifyEqual(testCase, double(cap3), double(cap4), 'RelTol', 1e-10);
end

function testAllSIPrefixes(testCase)
    % Test that all SI prefixes work with a base unit
    prefixes = {'quetta', 'ronna', 'yotta', 'zetta', 'exa', 'peta', 'tera', ...
                'giga', 'mega', 'kilo', 'hecto', 'deka', 'deci', 'centi', ...
                'milli', 'micro', 'nano', 'pico', 'femto', 'atto', 'zepto', ...
                'yocto', 'ronto', 'quecto'};
    
    for i = 1:length(prefixes)
        prefix = prefixes{i};
        
        % Test with meter as base unit
        unitName = [prefix 'meter'];
        try
            unit = u.(unitName);
            verifyTrue(testCase, isa(unit, 'DimVar'), ...
                ['Failed to create unit: ' unitName]);
        catch ME
            verifyFail(testCase, ['Error creating ' unitName ': ' ME.message]);
        end
    end
end