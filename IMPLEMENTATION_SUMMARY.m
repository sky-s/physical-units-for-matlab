%% SUMMARY: SI Prefix Support Implementation
%
% This file documents the implementation of automatic SI prefix support
% for the Physical Units for MATLAB toolbox.
%
% PROBLEM SOLVED:
% The issue requested the ability to use any SI prefix with any unit,
% such as u.megawatt, u.kilojoule, etc., without explicitly defining
% every possible combination.
%
% SOLUTION IMPLEMENTED:
% 1. Added subsref method to u class (lines ~1191-1254 in u.m)
% 2. Created parseSIPrefix function (lines ~1258-1334 in u.m) 
% 3. Added comprehensive test suite (.tests/testScript_SIprefixes.m)
% 4. Updated documentation and examples
%
% HOW IT WORKS:
% 1. When u.someunit is accessed, subsref method is called
% 2. First checks if 'someunit' exists as a direct property
% 3. If not, calls parseSIPrefix to check for prefix + base unit
% 4. parseSIPrefix tries to match SI prefixes at start of name
% 5. If valid prefix + base unit found, creates scaled DimVar
% 6. Returns the new unit with custom display name set
%
% FEATURES:
% - All 24 SI prefixes supported (quetta to quecto)
% - Common abbreviations: k, G, m, u, n, p (avoiding conflicts)
% - Handles special unit types (DimVar, OffsetDimVar, constants)
% - Proper error handling for invalid combinations
% - Maintains complete backward compatibility
% - Custom display names automatically set
%
% EXAMPLES:
% u.kilowatt    -> 1000 * u.watt
% u.milliampere -> 0.001 * u.ampere  
% u.gigahertz   -> 1e9 * u.hertz
% u.nanometer   -> 1e-9 * u.meter
% u.kV          -> 1000 * u.volt (abbreviation)
%
% TESTING:
% Run the test suite with: runtests('testScript_SIprefixes')
% See example usage in: example_SI_prefixes.m
% Manual demo with: demo_SI_prefixes.m
%
% IMPLEMENTATION NOTES:
% - Existing units take precedence (u.M = molar, not megaanything)
% - Only non-conflicting abbreviations enabled
% - Longest prefix names matched first to avoid conflicts
% - Single-letter abbreviations checked more carefully
% - Works with chained operations and unit conversions
%
% This implementation satisfies the original issue requirements while
% maintaining the robust design principles of the existing toolbox.

disp('SI Prefix Support Implementation Summary');
disp('==========================================');
disp(' ');
disp('✓ Dynamic SI prefix support implemented');
disp('✓ All 24 SI prefixes supported');  
disp('✓ Common abbreviations available');
disp('✓ Conflict resolution with existing units');
disp('✓ Comprehensive test suite added');
disp('✓ Documentation and examples provided');
disp('✓ Backward compatibility maintained');
disp(' ');
disp('Examples:');
disp('  u.megawatt, u.kilojoule, u.milliampere');  
disp('  u.gigahertz, u.nanometer, u.picofarad');
disp('  u.kV, u.GHz, u.mA, u.uF, u.nF, u.pF');
disp(' ');
disp('See SI_PREFIX_SUPPORT.md for full documentation');
disp(' ');