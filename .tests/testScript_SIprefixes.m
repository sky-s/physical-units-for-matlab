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
    % Test that u.get('kilometer') works
    distance = 5 * u.get('kilometer');
    expected = 5 * u.kilo * u.meter;
    verifyEqual(testCase, double(distance), double(expected), 'RelTol', 1e-10);
end

function testMegaWatt(testCase)
    % Test that u.get('megawatt') works even though only u.watt is explicitly defined
    power = 2 * u.get('megawatt');
    expected = 2 * u.mega * u.watt;
    verifyEqual(testCase, double(power), double(expected), 'RelTol', 1e-10);
end

function testMilliAmpere(testCase)
    % Test that u.get('milliampere') works
    current = 100 * u.get('milliampere');
    expected = 100 * u.milli * u.ampere;
    verifyEqual(testCase, double(current), double(expected), 'RelTol', 1e-10);
end

function testGigaHertz(testCase)
    % Test that u.get('gigahertz') works
    freq = 3.5 * u.get('gigahertz');
    expected = 3.5 * u.giga * u.hertz;
    verifyEqual(testCase, double(freq), double(expected), 'RelTol', 1e-10);
end

function testMicroSecond(testCase)
    % Test that u.microsecond works using static method
    time = 10 * u.get('microsecond');
    expected = 10 * u.micro * u.second;
    verifyEqual(testCase, double(time), double(expected), 'RelTol', 1e-10);
end

function testNanoMeter(testCase)
    % Test that u.get('nanometer') works
    length = 500 * u.get('nanometer');
    expected = 500 * u.nano * u.meter;
    verifyEqual(testCase, double(length), double(expected), 'RelTol', 1e-10);
end

function testPicoFarad(testCase)
    % Test that u.get('picofarad') works
    capacitance = 22 * u.get('picofarad');
    expected = 22 * u.pico * u.farad;
    verifyEqual(testCase, double(capacitance), double(expected), 'RelTol', 1e-10);
end

function testTeraJoule(testCase)
    % Test that u.get('terajoule') works
    energy = 1.5 * u.get('terajoule');
    expected = 1.5 * u.tera * u.joule;
    verifyEqual(testCase, double(energy), double(expected), 'RelTol', 1e-10);
end

function testExistingUnitsStillWork(testCase)
    % Test that existing manually defined units still work
    distance1 = 1 * u.km;  % Existing manually defined
    distance2 = 1 * u.get('kilometer');  % New dynamically created
    
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
        invalid = u.get('kiloinvalidunit');
        verifyFail(testCase, 'Should have failed for invalid unit');
    catch ME
        % Should get an error
        verifyTrue(testCase, contains(ME.message, 'Invalid unit') || contains(ME.message, 'field'));
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
    power = u.get('megawatt');
    
    % Check that the custom display is set
    verifyTrue(testCase, isa(power, 'DimVar'));
    if isa(power, 'DimVar')
        verifyEqual(testCase, power.customDisplay, 'megawatt');
    end
end

function testChainedOperations(testCase)
    % Test that prefixed units work in calculations
    voltage = 5 * u.get('kilovolt');
    current = 10 * u.get('milliampere');
    power = voltage * current;
    
    % This should give us watts
    expected_watts = 5 * 1e3 * 10 * 1e-3;  % 5000V * 0.01A = 50W
    verifyEqual(testCase, double(power), expected_watts, 'RelTol', 1e-10);
end

function testAbbreviatedPrefixes(testCase)
    % Test that common prefix abbreviations work (avoiding conflicts)
    
    % Test kV (kilovolt) - using existing constant
    voltage1 = 12 * u.kV;
    voltage2 = 12 * u.get('kilovolt');
    verifyEqual(testCase, double(voltage1), double(voltage2), 'RelTol', 1e-10);
    
    % Test GHz (gigahertz) - using existing constant
    freq1 = 2.4 * u.GHz;
    freq2 = 2.4 * u.get('gigahertz');
    verifyEqual(testCase, double(freq1), double(freq2), 'RelTol', 1e-10);
    
    % Test mA (milliampere) - using existing constant
    current1 = 100 * u.mA;
    current2 = 100 * u.get('milliampere');
    verifyEqual(testCase, double(current1), double(current2), 'RelTol', 1e-10);
    
    % Test nF (nanofarad) - using existing constant
    cap1 = 22 * u.nF;
    cap2 = 22 * u.get('nanofarad');
    verifyEqual(testCase, double(cap1), double(cap2), 'RelTol', 1e-10);
    
    % Test pF (picofarad) - using existing constant
    cap3 = 100 * u.pF;
    cap4 = 100 * u.get('picofarad');
    verifyEqual(testCase, double(cap3), double(cap4), 'RelTol', 1e-10);
    
    % Test uF (microfarad) - using existing constant
    cap5 = 10 * u.uF;
    cap6 = 10 * u.get('microfarad');
    verifyEqual(testCase, double(cap5), double(cap6), 'RelTol', 1e-10);
end

function testConflictAvoidance(testCase)
    % Test that existing units take precedence over prefix combinations
    
    % Test that u.M returns the existing Molar unit, not mega + something
    molar = u.M;
    verifyTrue(testCase, isa(molar, 'DimVar'));
    % M should be molar concentration (mol/L)
    expected = u.mol / u.L;
    verifyEqual(testCase, double(molar), double(expected), 'RelTol', 1e-10);
    
    % Test that u.T returns tesla, not tera + something
    tesla = u.T;
    verifyTrue(testCase, isa(tesla, 'DimVar'));
    % T should be tesla (N/(A*m))
    expected = u.N / (u.A * u.m);
    verifyEqual(testCase, double(tesla), double(expected), 'RelTol', 1e-10);
    
    % Test that existing km still works (manually defined)
    km_existing = u.km;
    km_dynamic = u.get('kilometer');
    verifyEqual(testCase, double(km_existing), double(km_dynamic), 'RelTol', 1e-10);
end

function testAllSIPrefixes(testCase)
    % Test that all SI prefixes work with a base unit using u.get()
    prefixes = {'quetta', 'ronna', 'yotta', 'zetta', 'exa', 'peta', 'tera', ...
                'giga', 'mega', 'kilo', 'hecto', 'deka', 'deci', 'centi', ...
                'milli', 'micro', 'nano', 'pico', 'femto', 'atto', 'zepto', ...
                'yocto', 'ronto', 'quecto'};
    
    for i = 1:length(prefixes)
        prefix = prefixes{i};
        
        % Test with meter as base unit
        unitName = [prefix 'meter'];
        try
            unit = u.get(unitName);
            verifyTrue(testCase, isa(unit, 'DimVar'), ...
                ['Failed to create unit: ' unitName]);
        catch ME
            verifyFail(testCase, ['Error creating ' unitName ': ' ME.message]);
        end
    end
end

function testNonSIUnitsWithPrefixes(testCase)
    % Test that SI prefixes work with non-SI units (imperial, other units)
    % This demonstrates the key capability: prefix + ANY unit, not just SI
    
    % Test kiloacre - combining SI prefix with imperial area unit
    area = 2.5 * u.get('kiloacre');
    expected = 2.5 * u.kilo * u.acre;
    verifyEqual(testCase, double(area), double(expected), 'RelTol', 1e-10);
    
    % Test nanoinch - combining SI prefix with imperial length unit
    length = 500 * u.get('nanoinch');
    expected = 500 * u.nano * u.inch;
    verifyEqual(testCase, double(length), double(expected), 'RelTol', 1e-10);
    
    % Test megapound - combining SI prefix with imperial mass unit
    mass = 1.2 * u.get('megapound');
    expected = 1.2 * u.mega * u.pound;
    verifyEqual(testCase, double(mass), double(expected), 'RelTol', 1e-10);
    
    % Test microinch - very small imperial unit (static method)
    precision = 25 * u.microinch();
    expected = 25 * u.micro * u.inch;
    verifyEqual(testCase, double(precision), double(expected), 'RelTol', 1e-10);
    
    % Test kilogallon - liquid volume (static method)
    volume = 5 * u.kilogallon();
    expected = 5 * u.kilo * u.gallon;
    verifyEqual(testCase, double(volume), double(expected), 'RelTol', 1e-10);
    
    % Test megabyte (should work for digital units too)
    data = 100 * u.get('megabyte');
    expected = 100 * u.mega * u.byte;
    verifyEqual(testCase, double(data), double(expected), 'RelTol', 1e-10);
    
    % Test milligallon (small liquid volume)
    droplet = 2.5 * u.get('milligallon');
    expected = 2.5 * u.milli * u.gallon;
    verifyEqual(testCase, double(droplet), double(expected), 'RelTol', 1e-10);
    
    % Test display names are set correctly for non-SI units
    verifyEqual(testCase, area.customDisplay, 'kiloacre');
    verifyEqual(testCase, length.customDisplay, 'nanoinch');
    verifyEqual(testCase, mass.customDisplay, 'megapound');
end

function testPrefixWithImperialCalculations(testCase)
    % Test that prefixed non-SI units work in real calculations
    
    % Area calculation: convert kiloacre to square miles
    area_kiloacre = 1 * u.get('kiloacre');
    area_sqmi = area_kiloacre / u.sqmi;
    
    % 1 kiloacre = 1000 acres, 1 acre = 43560 sqft, 1 sqmi = 5280^2 sqft
    expected = 1000 * 43560 / (5280^2);
    verifyEqual(testCase, double(area_sqmi), expected, 'RelTol', 1e-10);
    
    % Precision calculation: nanoinch to millimeter
    precision_nanoinch = 1000 * u.get('nanoinch');
    precision_mm = precision_nanoinch / u.mm;
    
    % 1 inch = 25.4 mm, so 1000 nanoinch = 1000 * 1e-9 * 25.4 mm
    expected = 1000 * 1e-9 * 25.4;
    verifyEqual(testCase, double(precision_mm), expected, 'RelTol', 1e-10);
end

function testDisplayParser(testCase)
    % Test that displayparser works correctly with prefixed units
    
    % Create a prefixed unit
    precision = 25 * u.get('microinch');
    
    % Test that displayparser can handle it
    [dispVal, unitStr, numString, denString, labelStr] = displayparser(precision);
    
    % Should return the value and the custom display name
    verifyEqual(testCase, dispVal, 25, 'RelTol', 1e-10);
    verifyEqual(testCase, unitStr, 'microinch');
    verifyEqual(testCase, numString, 'microinch');
    verifyTrue(testCase, isempty(denString));
    
    % Test with another prefixed unit  
    power = 2.5 * u.get('megawatt');
    [dispVal2, unitStr2] = displayparser(power);
    
    verifyEqual(testCase, dispVal2, 2.5, 'RelTol', 1e-10);
    verifyEqual(testCase, unitStr2, 'megawatt');
    
    % Test display methods don't throw errors
    % Just verify they run without error
    try
        evalc('disp(precision)'); % capture output to avoid cluttering test results
        evalc('display(power)');
        success = true;
    catch
        success = false;
    end
    verifyTrue(testCase, success, 'Display methods should work with prefixed units');
end

function testDisplayParserWithStrTou(testCase)
    % Test displayparser integration when str2u can't handle prefixed units
    
    % Create a unit and manually set a prefixed custom display
    length_unit = 5 * u.inch;
    length_unit = scd(length_unit, 'microinch'); % Set custom display to a prefixed unit
    
    % This should fallback to u.get() when str2u fails
    [dispVal, unitStr] = displayparser(length_unit);
    
    % The displayparser should convert the value to microinch units
    expected_val = 5 * 1e6; % 5 inches = 5 million microinches
    verifyEqual(testCase, dispVal, expected_val, 'RelTol', 1e-10);
    verifyEqual(testCase, unitStr, 'microinch');
end
end