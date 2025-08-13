# SI Prefix Support

## New Feature: Automatic SI Prefix Support

The Physical Units for MATLAB toolbox now supports using any SI prefix with any unit automatically, without requiring explicit definitions for each combination.

### Usage

Simply combine any SI prefix name with any base unit:

```matlab
% Length units
distance = 5 * u.kilometer;     % 5 km
wavelength = 500 * u.nanometer; % 500 nm
thickness = 0.1 * u.millimeter; % 0.1 mm

% Electrical units
voltage = 12 * u.kilovolt;      % 12 kV
current = 50 * u.milliampere;   % 50 mA
capacitance = 22 * u.picofarad; % 22 pF

% Power and energy
power = 2 * u.megawatt;         % 2 MW
energy = 10 * u.kilojoule;      % 10 kJ

% Frequency
freq = 2.4 * u.gigahertz;       % 2.4 GHz
```

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
```

Note: `M` (mega) and `T` (tera) abbreviations are not supported as they conflict with existing units (Molar and Tesla).

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
mass = 5 * u.kilogram;          % Works (already existed)
mass = 5 * u.megagram;          % Works (new)
mass = 5 * u.gigagram;          % Works (new)
mass = 5 * u.milligram;         % Works (already existed)
mass = 5 * u.microgram;         % Works (already existed)
mass = 5 * u.nanogram;          % Works (new)
```

This feature greatly expands the available units without cluttering the codebase with explicit definitions for every possible combination.