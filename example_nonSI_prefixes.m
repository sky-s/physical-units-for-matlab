%% Examples of SI Prefixes with Non-SI Units using Static Methods
% This demonstrates the key capability: SI prefixes work with ANY unit,
% not just SI base units. You can combine any SI prefix with any base unit.
% Use u.prefix('unit') methods for maximum flexibility.

% Clear the class to ensure fresh state
clear u

fprintf('=== SI Prefixes with Imperial/US Customary Units ===\n\n');

%% Area measurements with imperial units
fprintf('Area Examples:\n');
large_area = 2.5 * u.kilo('acre');
fprintf('Large farm: %.1f kiloacre = %.0f acre\n', ...
    displayingvalue(large_area, 'kiloacre'), ...
    displayingvalue(large_area, 'acre'));

small_plot = 500 * u.milli('acre');
fprintf('Small plot: %.0f milliacre = %.2f acre\n', ...
    displayingvalue(small_plot, 'milliacre'), ...
    displayingvalue(small_plot, 'acre'));

fprintf('\n');

%% Length measurements with imperial units  
fprintf('Length/Precision Examples:\n');
thickness = 250 * u.nano('inch');
fprintf('Surface roughness: %.0f nanoinch = %.3f micrometer\n', ...
    displayingvalue(thickness, 'nanoinch'), ...
    displayingvalue(thickness, 'micrometer'));

pipe_length = 2.5 * u.kilo('inch');
fprintf('Pipe length: %.1f kiloinch = %.1f foot\n', ...
    displayingvalue(pipe_length, 'kiloinch'), ...
    displayingvalue(pipe_length, 'foot'));

tolerance = 5 * u.micro('inch');  % Precision manufacturing
fprintf('Machining tolerance: %.0f microinch = %.2f micrometer\n', ...
    displayingvalue(tolerance, 'microinch'), ...
    displayingvalue(tolerance, 'micrometer'));

fprintf('\n');

%% Mass measurements
fprintf('Mass Examples:\n');
cargo = 1.2 * u.mega('pound');  % Large mass
fprintf('Ship cargo: %.1f megapound = %.0f ton\n', ...
    displayingvalue(cargo, 'megapound'), ...
    displayingvalue(cargo, 'ton'));

sample = 25 * u.milli('gram');  % Small mass
fprintf('Lab sample: %.0f milligram = %.3f grain\n', ...
    displayingvalue(sample, 'milligram'), ...
    displayingvalue(sample, 'grain'));

fprintf('\n');

%% Volume measurements
fprintf('Volume Examples:\n');
fuel_tank = 50 * u.kilo('gallon');  % Large volume
fprintf('Fuel tank: %.0f kilogallon = %.0f liter\n', ...
    displayingvalue(fuel_tank, 'kilogallon'), ...
    displayingvalue(fuel_tank, 'liter'));

droplet = 0.5 * u.milli('gallon');
fprintf('Paint droplet: %.1f milligallon = %.2f milliliter\n', ...
    displayingvalue(droplet, 'milligallon'), ...
    displayingvalue(droplet, 'milliliter'));

fprintf('\n');

%% Digital units (non-imperial, but also non-SI)
fprintf('Digital Storage Examples:\n');
storage = 2 * u.tera('byte');  % Large storage
fprintf('Storage drive: %.0f terabyte = %.0f gigabyte\n', ...
    displayingvalue(storage, 'terabyte'), ...
    displayingvalue(storage, 'gigabyte'));

cache = 512 * u.kilo('byte');
fprintf('Cache size: %.0f kilobyte = %.0f byte\n', ...
    displayingvalue(cache, 'kilobyte'), ...
    displayingvalue(cache, 'byte'));

fprintf('\n');

%% Display name control
fprintf('Display Name Control Examples:\n');
voltage_short = 12 * u.kilo('volt');  % Default: short name 'kV'
voltage_full = 12 * u.kilo('volt', true);  % Full name: 'kilovolt'
fprintf('Short name: %s\n', char(voltage_short));
fprintf('Full name: %s\n', char(voltage_full));

area_short = 100 * u.micro('acre');  % Default: short name 'uacre' 
area_full = 100 * u.micro('acre', true);  % Full name: 'microacre'
fprintf('Short name: %s\n', char(area_short));
fprintf('Full name: %s\n', char(area_full));

fprintf('\n');

%% Real-world calculation examples
fprintf('=== Real-World Calculations ===\n\n');

% Convert farm area from kiloacre to hectare
farm_area = 5 * u.kilo('acre');
farm_hectare = farm_area / u.hectare;
fprintf('Farm conversion: %.0f kiloacre = %.0f hectare\n', ...
    displayingvalue(farm_area, 'kiloacre'), double(farm_hectare));

% Precision manufacturing: microinch to nanometer
surface_finish = 32 * u.micro('inch');
surface_nm = surface_finish / u.nanometer;
fprintf('Surface finish: %.0f microinch = %.0f nanometer\n', ...
    displayingvalue(surface_finish, 'microinch'), double(surface_nm));

% Fuel efficiency with mixed units
distance = 500 * u.mile;
fuel = 25 * u.kilo('gallon');
efficiency = distance / fuel;
fprintf('Fuel efficiency: %.0f mile / %.0f kilogallon = %.3f mile/gallon\n', ...
    displayingvalue(distance, 'mile'), ...
    displayingvalue(fuel, 'kilogallon'), ...
    displayingvalue(efficiency, 'mpg'));

% Demonstrate prefix value access
fprintf('\nPrefix values (no arguments):\n');
fprintf('u.kilo() = %g\n', u.kilo());
fprintf('u.micro() = %g\n', u.micro());
fprintf('u.nano() = %g\n', u.nano());

fprintf('\n=== Summary ===\n');
fprintf('SI prefix static methods work with ANY unit:\n');
fprintf('• u.prefix(''unit'') creates prefixed unit with smart naming\n');
fprintf('• u.prefix(''unit'', true) forces full prefix name\n');
fprintf('• u.prefix() returns prefix multiplier value\n');
fprintf('• Works with Imperial: u.kilo(''acre''), u.nano(''inch''), u.mega(''pound'')\n');
fprintf('• Works with US customary: u.kilo(''gallon''), u.milli(''gallon'')\n');  
fprintf('• Works with digital: u.tera(''byte''), u.kilo(''byte'')\n');
fprintf('• Works with any unit in the toolbox!\n');