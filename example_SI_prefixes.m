% Example usage of SI prefix functionality using new static methods
% Demonstrates the new capability to use any SI prefix with any unit
% Use u.prefix('unit') methods for maximum flexibility

%% Basic Examples
fprintf('=== SI Prefix Examples with Static Methods ===\n\n');

fprintf('Length units:\n');
fprintf('  1 kilometer = %g meters\n', double(1 * u.kilo('meter') / u.meter));
fprintf('  1 millimeter = %g meters\n', double(1 * u.milli('meter') / u.meter));
fprintf('  1 nanometer = %g meters\n', double(1 * u.nano('meter') / u.meter));

fprintf('\nElectrical units:\n');
fprintf('  1 kilovolt = %g volts\n', double(1 * u.kilo('volt') / u.volt));
fprintf('  1 milliampere = %g amperes\n', double(1 * u.milli('ampere') / u.ampere));
fprintf('  1 microfarad = %g farads\n', double(1 * u.micro('farad') / u.farad));

fprintf('\nPower and Energy:\n');
fprintf('  1 megawatt = %g watts\n', double(1 * u.mega('watt') / u.watt));
fprintf('  1 kilojoule = %g joules\n', double(1 * u.kilo('joule') / u.joule));
fprintf('  1 gigawatt = %g watts\n', double(1 * u.giga('watt') / u.watt));

fprintf('\nFrequency:\n');
fprintf('  1 kilohertz = %g hertz\n', double(1 * u.kilo('hertz') / u.hertz));
fprintf('  1 megahertz = %g hertz\n', double(1 * u.mega('hertz') / u.hertz));
fprintf('  1 gigahertz = %g hertz\n', double(1 * u.giga('hertz') / u.hertz));

%% Prefix Values (No Arguments)
fprintf('\n=== Prefix Values (No Arguments) ===\n\n');
fprintf('u.kilo() = %g\n', u.kilo());
fprintf('u.mega() = %g\n', u.mega());
fprintf('u.micro() = %g\n', u.micro());
fprintf('u.nano() = %g\n', u.nano());

%% Display Name Control
fprintf('\n=== Display Name Control ===\n\n');
fprintf('Short names (default):\n');
short_kV = u.kilo('volt');
fprintf('  u.kilo(''volt'') displays as: %s\n', char(short_kV));

fprintf('\nFull names (second argument = true):\n');
full_kV = u.kilo('volt', true);
fprintf('  u.kilo(''volt'', true) displays as: %s\n', char(full_kV));

%% Practical Calculations
fprintf('\n=== Practical Calculations ===\n\n');

% Power calculation
voltage = 12 * u.kilo('volt');      % 12 kV
current = 5 * u.milli('ampere');    % 5 mA
power = voltage * current;
fprintf('Power = %s × %s = %s\n', char(voltage), char(current), char(power));

% Energy calculation
time = 2 * u.hour;
energy = power * time;
fprintf('Energy = %s × %s = %s\n', char(power), char(time), char(energy));

% Frequency and wavelength
freq = 2.4 * u.giga('hertz');       % 2.4 GHz (WiFi frequency)
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
    try
        unit_val = u.(prefixes{i})('meter');  % Use prefix method
        actual = double(unit_val / u.meter);
        expected = values(i);
        if abs(actual - expected) < 1e-15 * max(abs(actual), abs(expected))
            status = '✓';
        else
            status = '✗';
        end
        fprintf('%s %smeter = %g m\n', status, prefixes{i}, actual);
    catch ME
        fprintf('✗ %smeter failed: %s\n', prefixes{i}, ME.message);
    end
end

%% Working with Existing Units
fprintf('\n=== Compatibility with Existing Units ===\n');

% These existing units still work normally
fprintf('Existing constant units:\n');
fprintf('  u.kV = %s (existing constant)\n', char(u.kV));
fprintf('  u.mA = %s (existing constant)\n', char(u.mA)); 
fprintf('  u.GHz = %s (existing constant)\n', char(u.GHz));

% And these can also be created with prefix methods
fprintf('\nEquivalent with prefix methods:\n');
fprintf('  u.kilo(''V'') = %s\n', char(u.kilo('V')));
fprintf('  u.milli(''A'') = %s\n', char(u.milli('A')));
fprintf('  u.giga(''Hz'') = %s\n', char(u.giga('Hz')));

fprintf('\n=== Summary ===\n');
fprintf('The SI prefix static methods allow you to:\n');
fprintf('1. Get prefix values: u.kilo() returns 1000\n');
fprintf('2. Create prefixed units: u.kilo(''volt'') creates kV\n');
fprintf('3. Control display names: u.kilo(''volt'', true) creates ''kilovolt''\n');
fprintf('This greatly expands the available units while maintaining\n');
fprintf('compatibility with existing units and providing discoverable methods.\n');