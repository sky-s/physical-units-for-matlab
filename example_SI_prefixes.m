% Example usage of SI prefix functionality
% Demonstrates the new capability to use any SI prefix with any unit

%% Basic Examples
fprintf('=== SI Prefix Examples ===\n\n');

% These all work automatically now:
fprintf('Length units:\n');
fprintf('  1 kilometer = %g meters\n', double(1 * u.kilometer / u.meter));
fprintf('  1 millimeter = %g meters\n', double(1 * u.millimeter / u.meter));
fprintf('  1 nanometer = %g meters\n', double(1 * u.nanometer / u.meter));

fprintf('\nElectrical units:\n');
fprintf('  1 kilovolt = %g volts\n', double(1 * u.kilovolt / u.volt));
fprintf('  1 milliampere = %g amperes\n', double(1 * u.milliampere / u.ampere));
fprintf('  1 microfarad = %g farads\n', double(1 * u.microfarad / u.farad));

fprintf('\nPower and Energy:\n');
fprintf('  1 megawatt = %g watts\n', double(1 * u.megawatt / u.watt));
fprintf('  1 kilojoule = %g joules\n', double(1 * u.kilojoule / u.joule));
fprintf('  1 gigawatt = %g watts\n', double(1 * u.gigawatt / u.watt));

fprintf('\nFrequency:\n');
fprintf('  1 kilohertz = %g hertz\n', double(1 * u.kilohertz / u.hertz));
fprintf('  1 megahertz = %g hertz\n', double(1 * u.megahertz / u.hertz));
fprintf('  1 gigahertz = %g hertz\n', double(1 * u.gigahertz / u.hertz));

%% Practical Calculations
fprintf('\n=== Practical Calculations ===\n\n');

% Power calculation
voltage = 12 * u.kilovolt;      % 12 kV
current = 5 * u.milliampere;    % 5 mA
power = voltage * current;
fprintf('Power = %s × %s = %s\n', char(voltage), char(current), char(power));

% Energy calculation
time = 2 * u.hour;
energy = power * time;
fprintf('Energy = %s × %s = %s\n', char(power), char(time), char(energy));

% Frequency and wavelength
freq = 2.4 * u.gigahertz;       % 2.4 GHz (WiFi frequency)
wavelength = u.c / freq;        % Speed of light / frequency
fprintf('WiFi wavelength = %s / %s = %s\n', char(u.c), char(freq), char(wavelength));

%% All SI Prefixes Demonstration
fprintf('\n=== All SI Prefixes with Meter ===\n');

prefixes = {'quetta', 'ronna', 'yotta', 'zetta', 'exa', 'peta', 'tera', ...
            'giga', 'mega', 'kilo', 'hecto', 'deka', 'deci', 'centi', ...
            'milli', 'micro', 'nano', 'pico', 'femto', 'atto', 'zepto', ...
            'yocto', 'ronto', 'quecto'};

values = [1e30, 1e27, 1e24, 1e21, 1e18, 1e15, 1e12, 1e9, 1e6, 1e3, ...
          1e2, 1e1, 1e-1, 1e-2, 1e-3, 1e-6, 1e-9, 1e-12, 1e-15, ...
          1e-18, 1e-21, 1e-24, 1e-27, 1e-30];

for i = 1:length(prefixes)
    unit_name = [prefixes{i} 'meter'];
    try
        unit_val = u.(unit_name);
        actual = double(unit_val / u.meter);
        expected = values(i);
        if abs(actual - expected) < 1e-15 * max(abs(actual), abs(expected))
            status = '✓';
        else
            status = '✗';
        end
        fprintf('%s %s = %g m\n', status, unit_name, actual);
    catch ME
        fprintf('✗ %s failed: %s\n', unit_name, ME.message);
    end
end

%% Abbreviated Prefixes
fprintf('\n=== Abbreviated Prefixes (Non-conflicting) ===\n');

try
    fprintf('k (kilo): 1 kV = %g V\n', double(u.kV / u.V));
    fprintf('G (giga): 1 GHz = %g Hz\n', double(u.GHz / u.Hz));
    fprintf('m (milli): 1 mA = %g A\n', double(u.mA / u.A));
    fprintf('u (micro): 1 uF = %g F\n', double(u.uF / u.F));
    fprintf('n (nano): 1 nF = %g F\n', double(u.nF / u.F));
    fprintf('p (pico): 1 pF = %g F\n', double(u.pF / u.F));
    
    % Note: M and T abbreviations are not supported due to conflicts:
    fprintf('\nNote: M and T abbreviations avoided due to existing units:\n');
    fprintf('  u.M = %s (molar concentration)\n', char(u.M));
    fprintf('  u.T = %s (tesla)\n', char(u.T));
    fprintf('  Use u.megawatt, u.terawatt instead of MW, TW\n');
catch ME
    fprintf('Error with abbreviated prefixes: %s\n', ME.message);
end

fprintf('\n=== Summary ===\n');
fprintf('The SI prefix functionality allows you to use any SI prefix\n');
fprintf('with any base unit, greatly expanding the available units\n');
fprintf('without cluttering the u class with explicit definitions.\n');