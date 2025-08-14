# SI Prefix Support

## New Feature: Automatic SI Prefix Support

The Physical Units for MATLAB toolbox now supports using any SI prefix with **ANY unit** automatically, without requiring explicit definitions for each combination. This works with SI units, imperial units, digital units, and any other unit in the toolbox.

## MATLAB Limitation and Solution

Due to MATLAB's class system limitations with constant properties, the direct `u.prefixunit` syntax (like `u.microinch`) cannot work. Instead, use one of these approaches:

### Method 1: u.get() Function (Recommended)

```matlab
% Use u.get('prefixunit') for any combination
distance = 5 * u.get('kilometer');     % 5 km
voltage = 12 * u.get('kilovolt');      % 12 kV  
power = 2 * u.get('megawatt');         % 2 MW
wavelength = 500 * u.get('nanometer'); % 500 nm

% Imperial/US Customary units
area = 2.5 * u.get('kiloacre');        % 2.5 thousand acres
thickness = 250 * u.get('nanoinch');   % 250 billionths of an inch
mass = 1.2 * u.get('megapound');       % 1.2 million pounds
precision = 5 * u.get('microinch');    % 5 millionths of an inch

% Volume measurements
tank = 50 * u.get('kilogallon');       % 50 thousand gallons
droplet = 0.5 * u.get('milligallon');  % half a thousandth of a gallon

% Digital units
storage = 2 * u.get('terabyte');       % 2 TB
cache = 512 * u.get('kilobyte');       % 512 KB
```

### Method 2: Static Methods for Common Units

```matlab
% Common prefixed units are available as static methods
precision = 25 * u.microinch();       % Works as static method
area = 2.5 * u.kiloacre();            % Works as static method
storage = 1 * u.terabyte();           % Works as static method
mass = 1.2 * u.megapound();           % Works as static method
volume = 50 * u.kilogallon();         % Works as static method
length = 500 * u.nanoinch();          % Works as static method
```

### Key Capability: Works with ANY Unit

The power of this feature is that SI prefixes work with **every unit** in the toolbox:

- **SI units**: meter, gram, second, ampere, kelvin, etc.
- **Imperial units**: inch, foot, yard, mile, acre, pound, gallon, etc.
- **Digital units**: byte, bit
- **Time units**: hour, day, year
- **Specialized units**: hertz, pascal, joule, watt, ohm, etc.

### Supported SI Prefixes

All SI prefixes are supported:

- **Large**: quetta (1e30), ronna (1e27), yotta (1e24), zetta (1e21), exa (1e18), peta (1e15), tera (1e12), giga (1e9), mega (1e6), kilo (1e3), hecto (1e2), deka (1e1)
- **Small**: deci (1e-1), centi (1e-2), milli (1e-3), micro (1e-6), nano (1e-9), pico (1e-12), femto (1e-15), atto (1e-18), zepto (1e-21), yocto (1e-24), ronto (1e-27), quecto (1e-30)

### Abbreviations

Common abbreviations are supported where they don't conflict with existing units:

```matlab
voltage = 5 * u.kV;    % kilovolt (k = kilo)
freq = 1 * u.GHz;      % gigahertz (G = giga)
current = 10 * u.mA;   % milliampere (m = milli)
cap = 100 * u.uF;      % microfarad (u = micro)
area = 2 * u.kacre;    % kiloacre (k + acre)
```

Note: `M` (mega) and `T` (tera) abbreviations are not supported as they conflict with existing units (Molar and Tesla).

### Real-World Examples

#### Imperial Unit Conversions
```matlab
% Farm area management
farm = 5 * u.kiloacre;                    % 5000 acres
farm_hectares = farm / u.hectare;         % Convert to hectares

% Precision manufacturing  
tolerance = 32 * u.microinch;             % 32 millionths of an inch
tolerance_nm = tolerance / u.nanometer;   % Convert to nanometers
```

#### Mixed Unit Calculations
```matlab
% Fuel efficiency calculation
distance = 500 * u.mile;
fuel = 25 * u.kilogallon;                 % 25,000 gallons
efficiency = distance / fuel;             % miles per gallon equivalent
```

#### Digital Storage
```matlab
drive = 4 * u.terabyte;                   % 4 TB drive
file = 250 * u.megabyte;                  % 250 MB file
files_per_drive = drive / file;           % How many files fit
```

### Backward Compatibility

All existing units continue to work exactly as before. The new functionality only adds support for previously undefined prefix-unit combinations.

### Examples

#### Power Calculation
```matlab
voltage = 12 * u.kilovolt;      % 12 kV
current = 5 * u.milliampere;    % 5 mA
power = voltage * current;      % = 60 watts
```

#### Unit Conversion
```matlab
% These are equivalent
dist1 = 1 * u.kilometer;
dist2 = 1 * u.kilo * u.meter;
dist3 = 1000 * u.meter;
```

#### All prefixes work with any unit
```matlab
% With SI base units
mass = 5 * u.kilogram;          % Works (already existed)
mass = 5 * u.megagram;          % Works (new)
mass = 5 * u.gigagram;          % Works (new)

% With imperial units  
area = 10 * u.acre;             % Works (already existed)
area = 10 * u.kiloacre;         % Works (new)
area = 10 * u.milliacre;        % Works (new)

% With digital units
data = 1 * u.byte;              % Works (already existed)
data = 1 * u.kilobyte;          % Works (already existed)
data = 1 * u.terabyte;          % Works (new)
```

This feature greatly expands the available units (from hundreds to potentially thousands) without cluttering the codebase with explicit definitions for every possible combination.