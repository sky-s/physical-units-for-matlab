%% Examples of SI Prefixes with Non-SI Units
% This demonstrates the key capability: SI prefixes work with ANY unit,
% not just SI base units. You can combine any SI prefix with any base unit.
% Due to MATLAB limitations, use u.get('prefixunit') or static methods.

% Clear the class to ensure fresh state
clear u

fprintf('=== SI Prefixes with Imperial/US Customary Units ===\n\n');

%% Area measurements with imperial units
fprintf('Area Examples:\n');
large_area = 2.5 * u.get('kiloacre');
fprintf('Large farm: %.1f kiloacre = %.0f acre\n', ...
    displayingvalue(large_area, 'kiloacre'), ...
    displayingvalue(large_area, 'acre'));

small_plot = 500 * u.get('milliacre');
fprintf('Small plot: %.0f milliacre = %.2f acre\n', ...
    displayingvalue(small_plot, 'milliacre'), ...
    displayingvalue(small_plot, 'acre'));

fprintf('\n');

%% Length measurements with imperial units  
fprintf('Length/Precision Examples:\n');
thickness = 250 * u.get('nanoinch');
fprintf('Surface roughness: %.0f nanoinch = %.3f micrometer\n', ...
    displayingvalue(thickness, 'nanoinch'), ...
    displayingvalue(thickness, 'micrometer'));

pipe_length = 2.5 * u.get('kiloinch');
fprintf('Pipe length: %.1f kiloinch = %.1f foot\n', ...
    displayingvalue(pipe_length, 'kiloinch'), ...
    displayingvalue(pipe_length, 'foot'));

tolerance = 5 * u.microinch();  % Static method for common units
fprintf('Machining tolerance: %.0f microinch = %.2f micrometer\n', ...
    displayingvalue(tolerance, 'microinch'), ...
    displayingvalue(tolerance, 'micrometer'));

fprintf('\n');

%% Mass measurements
fprintf('Mass Examples:\n');
cargo = 1.2 * u.megapound();  % Static method for common units
fprintf('Ship cargo: %.1f megapound = %.0f ton\n', ...
    displayingvalue(cargo, 'megapound'), ...
    displayingvalue(cargo, 'ton'));

sample = 25 * u.milligram;  % This was already defined, but shows consistency
fprintf('Lab sample: %.0f milligram = %.3f grain\n', ...
    displayingvalue(sample, 'milligram'), ...
    displayingvalue(sample, 'grain'));

fprintf('\n');

%% Volume measurements
fprintf('Volume Examples:\n');
fuel_tank = 50 * u.kilogallon();  % Static method for common units
fprintf('Fuel tank: %.0f kilogallon = %.0f liter\n', ...
    displayingvalue(fuel_tank, 'kilogallon'), ...
    displayingvalue(fuel_tank, 'liter'));

droplet = 0.5 * u.get('milligallon');
fprintf('Paint droplet: %.1f milligallon = %.2f milliliter\n', ...
    displayingvalue(droplet, 'milligallon'), ...
    displayingvalue(droplet, 'milliliter'));

fprintf('\n');

%% Digital units (non-imperial, but also non-SI)
fprintf('Digital Storage Examples:\n');
storage = 2 * u.terabyte();  % Static method for common units
fprintf('Storage drive: %.0f terabyte = %.0f gigabyte\n', ...
    displayingvalue(storage, 'terabyte'), ...
    displayingvalue(storage, 'gigabyte'));

cache = 512 * u.get('kilobyte');
fprintf('Cache size: %.0f kilobyte = %.0f byte\n', ...
    displayingvalue(cache, 'kilobyte'), ...
    displayingvalue(cache, 'byte'));

fprintf('\n');

# Real-world calculation examples
fprintf('=== Real-World Calculations ===\n\n');

% Convert farm area from kiloacre to hectare
farm_area = 5 * u.get('kiloacre');
farm_hectare = farm_area / u.hectare;
fprintf('Farm conversion: %.0f kiloacre = %.0f hectare\n', ...
    displayingvalue(farm_area, 'kiloacre'), double(farm_hectare));

% Precision manufacturing: microinch to nanometer
surface_finish = 32 * u.microinch();  % Static method
surface_nm = surface_finish / u.nanometer;
fprintf('Surface finish: %.0f microinch = %.0f nanometer\n', ...
    displayingvalue(surface_finish, 'microinch'), double(surface_nm));

% Fuel efficiency with mixed units
distance = 500 * u.mile;
fuel = 25 * u.get('kilogallon');
efficiency = distance / fuel;
fprintf('Fuel efficiency: %.0f mile / %.0f kilogallon = %.3f mile/gallon\n', ...
    displayingvalue(distance, 'mile'), ...
    displayingvalue(fuel, 'kilogallon'), ...
    displayingvalue(efficiency, 'mpg'));

fprintf('\n=== Summary ===\n');
fprintf('SI prefixes (kilo, mega, nano, micro, etc.) work with ANY unit:\n');
fprintf('• Use u.get(''prefixunit'') for any combination\n');
fprintf('• Static methods available for common units: u.microinch(), u.kiloacre()\n');
fprintf('• Imperial units: kiloacre, nanoinch, megapound\n');
fprintf('• US customary: kilogallon, milligallon\n');  
fprintf('• Digital units: terabyte, kilobyte\n');
fprintf('• Time units: kilosecond, millisecond\n');
fprintf('• Any other unit in the toolbox!\n');