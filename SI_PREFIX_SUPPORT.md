# SI Prefix Support

## New Feature: Individual SI Prefix Static Methods

The Physical Units for MATLAB toolbox now supports using any SI prefix with **ANY unit** through individual static methods for each prefix. This works with SI units, imperial units, digital units, and any other unit in the toolbox.

## Individual Prefix Methods

Each SI prefix is now available as its own static method, providing three usage modes:

### Mode 1: Get Prefix Multiplier (No Arguments)
```matlab
multiplier = u.kilo();      % Returns 1000
multiplier = u.mega();      % Returns 1e6  
multiplier = u.micro();     % Returns 1e-6
multiplier = u.nano();      % Returns 1e-9
```

### Mode 2: Create Prefixed Unit (One Argument)
```matlab
% Use u.prefix('unit') to create prefixed units
distance = 5 * u.kilo('meter');     % 5 km
voltage = 12 * u.kilo('volt');      % 12 kV  
power = 2 * u.mega('watt');         % 2 MW
wavelength = 500 * u.nano('meter'); % 500 nm

% Imperial/US Customary units
area = 2.5 * u.kilo('acre');        % 2.5 thousand acres
thickness = 250 * u.nano('inch');   % 250 billionths of an inch
mass = 1.2 * u.mega('pound');       % 1.2 million pounds
precision = 5 * u.micro('inch');    % 5 millionths of an inch

% Volume measurements
tank = 50 * u.kilo('gallon');       % 50 thousand gallons
droplet = 0.5 * u.milli('gallon');  % half a thousandth of a gallon

% Digital units
storage = 2 * u.tera('byte');       % 2 TB
cache = 512 * u.kilo('byte');       % 512 KB
```

### Mode 3: Control Display Name (Two Arguments)
```matlab
% Second argument controls whether to use full prefix name
voltage_short = 12 * u.kilo('volt');       % Default: displays as 'kV'
voltage_full = 12 * u.kilo('volt', true);  % Displays as 'kilovolt'

area_short = 100 * u.micro('acre');        % Default: displays as 'uacre'
area_full = 100 * u.micro('acre', true);   % Displays as 'microacre'
```

### Key Capability: Works with ANY Unit

The power of this feature is that SI prefixes work with **every unit** in the toolbox:

- **SI units**: meter, gram, second, ampere, kelvin, etc.
- **Imperial units**: inch, foot, yard, mile, acre, pound, gallon, etc.
- **Digital units**: byte, bit
- **Time units**: hour, day, year
- **Specialized units**: hertz, pascal, joule, watt, ohm, etc.

### Supported SI Prefixes

All SI prefixes are available as individual static methods:

- **Large**: quetta (1e30), ronna (1e27), yotta (1e24), zetta (1e21), exa (1e18), peta (1e15), tera (1e12), giga (1e9), mega (1e6), kilo (1e3), hecto (1e2), deka (1e1)
- **Small**: deci (1e-1), centi (1e-2), milli (1e-3), micro (1e-6), nano (1e-9), pico (1e-12), femto (1e-15), atto (1e-18), zepto (1e-21), yocto (1e-24), ronto (1e-27), quecto (1e-30)

### Smart Display Names

The methods automatically choose appropriate display names to avoid conflicts:

```matlab
voltage1 = 5 * u.kilo('V');         % Displays as 'kV' (no conflict)
voltage2 = 5 * u.kilo('volt');      % Displays as 'kV' (smart abbreviation)
voltage3 = 5 * u.kilo('volt', true); % Displays as 'kilovolt' (forced full name)
```

### Compatibility with Existing Units

Existing units continue to work exactly as before:
```matlab
voltage = 5 * u.kV;                 % Existing constant still works
freq = 1 * u.GHz;                   % Existing constant still works
current = 10 * u.mA;                % Existing constant still works
```

And equivalent units can be created with prefix methods:
```matlab
voltage = 5 * u.kilo('V');          % Equivalent to u.kV
freq = 1 * u.giga('Hz');            % Equivalent to u.GHz  
current = 10 * u.milli('A');        % Equivalent to u.mA
```

### Real-World Examples

#### Imperial Unit Conversions
```matlab
% Farm area management
farm = 5 * u.kilo('acre');                  % 5000 acres
farm_hectares = farm / u.hectare;           % Convert to hectares

% Precision manufacturing  
tolerance = 32 * u.micro('inch');           % 32 millionths of an inch
tolerance_nm = tolerance / u.nanometer;     % Convert to nanometers
```

#### Mixed Unit Calculations
```matlab
% Fuel efficiency calculation
distance = 500 * u.mile;
fuel = 25 * u.kilo('gallon');               % 25,000 gallons
efficiency = distance / fuel;               % miles per gallon equivalent
```

#### Digital Storage
```matlab
drive = 4 * u.tera('byte');                 % 4 TB drive
file = 250 * u.mega('byte');                % 250 MB file
files_per_drive = drive / file;             % How many files fit
```

### Backward Compatibility

All existing units continue to work exactly as before. The new functionality only adds support for previously undefined prefix-unit combinations.

### Examples

#### Power Calculation
```matlab
voltage = 12 * u.kilo('volt');      % 12 kV
current = 5 * u.milli('ampere');    % 5 mA
power = voltage * current;          % = 60 watts
```

#### Unit Conversion
```matlab
% These are equivalent
dist1 = 1 * u.kilo('meter');
dist2 = 1 * u.kilo() * u.meter;
dist3 = 1000 * u.meter;
```

#### All prefixes work with any unit
```matlab
% With SI base units
mass = 5 * u.kilogram;              % Works (already existed)
mass = 5 * u.mega('gram');          % Works (new)
mass = 5 * u.giga('gram');          % Works (new)

% With imperial units  
area = 10 * u.acre;                 % Works (already existed)
area = 10 * u.kilo('acre');         % Works (new)
area = 10 * u.milli('acre');        % Works (new)

% With digital units
data = 1 * u.byte;                  % Works (already existed)
data = 1 * u.kilo('byte');          % Works (new)
data = 1 * u.tera('byte');          % Works (new)
```

This feature greatly expands the available units (from hundreds to potentially thousands) while providing discoverable static methods for all SI prefixes and maintaining full backward compatibility.