function baseUnitSystem = baseUnitSystem()
% baseUnitSystem  Sets the base set of fundamental units for variables with
% physical units. Copy this file to your project directory and modify to
% customize physical units for a particular project.
% 
%   Most of the time it will be unnecessary to use a custom baseUnitSystem, as
%   most users preferences can often be captured by defining only preferred
%   display units in displayUnits.
%   
%   baseUnitSystem should return a 9 by 2 cell array with characters in the
%   first column and doubles in the second. The format of the list must
%   correspond to the default base SI system (m, kg, s, A, K, mol, cd, bit, ¤)
%   and have the same order. The value in the second column is the factor by
%   which the user-provided base unit must be multiplied to have the value of
%   the base SI unit. For example, to use cm as the default length, the value in
%   the second column would be 100 (100 cm in one meter). See in-code examples.
% 
%   Alternatively, make baseUnitSystem return 'none' (or set a variable
%   baseUnitSystem = 'none' in the base workspace) to have the u class contain
%   non-DimVars. Use this option for when you need a bit more speed on code that
%   has already been developed and debugged by getting rid of the DimVar method
%   overhead, including e.g. unit compatibility checking.
% 
%   Remember to clear classes whenever changing baseUnitSystem.
% 
%   See also u, displayUnits.

%% Basic SI.
baseUnitSystem = {
    'm'        1
    'kg'       1
    's'        1
    'A'        1
    'K'        1
    'mol'      1
    'cd'       1
    'bit'      1
    '¤'        1 
    };

% return
%% Menu
unitSystem = 'SI'; % Change this value to select system.

switch upper(unitSystem)
    case {'SI' 'MKS' 'METRIC' 'INTERNATIONAL'}
        unitSystem = baseUnitSystem;
    case 'IPS'
        unitSystem = {
            'in'    1/0.0254
            'lbm'   2.2046
            's'     1
            'A'     1
            'R'     1.8};
    case 'IPSK'
        unitSystem = {
            'in'    1/0.0254
            'lbm'   2.2046};
    case {'FPS' 'EE' 'AE' 'IMPERIAL' 'AMERICAN'}
        % EE = English Engineering; AE = Absolute English
        unitSystem = {
            'ft'    1/0.3048
            'lbm'   2.2046
            's'     1
            'A'     1
            'R'     1.8};
    case 'FPSK'
        unitSystem = {
            'ft'    1/0.3048
            'lbm'   2.2046};
    case {'FSS' 'BG' 'GRAVITATIONAL FPS' 'TECHNICAL FPS'}
        % BG = British Gravitational
        unitSystem = {
            'ft'    1/0.3048
            'slug'  0.3048/4.4482216152605 %exact
            's'     1
            'A'     1
            'R'     1.8};
    case {'ISS' 'GRAVITATIONAL IPS' 'TECHNICAL IPS'}
        unitSystem = {
            'in'    1/0.0254
            'slug'  0.3048/4.4482216152605 %exact
            's'     1
            'A'     1
            'R'     1.8};
    case 'SLINCH'
        unitSystem = {
            'in'    1/0.0254
            'slinch'  0.3048/4.4482216152605/12 %exact
            's'     1
            'A'     1
            'R'     1.8};
    case {'VERBOSE' 'VERBOSE MKS' 'VERBOSE SI' 'SI VERBOSE' ...
            'MKS VERBOSE'}
        unitSystem = {
            'meter'     1
            'kilogram'  1
            'second'    1
            'ampere'    1
            'kelvin'    1
            'mole'      1
            'candela'   1
            'bit'       1
            'dollar'    1};
    case {'VERBOSE FPS' 'FPS VERBOSE' 'VERBOSE EE' 'EE VERBOSE'}
        unitSystem = {
            'foot'      1/0.3048
            'poundMass' 2.2046
            'second'    1
            'ampere'    1
            'Rankine'   1.8
            'mole'      1
            'candela'   1
            'bit'       1
            'dollar'    1};
    case 'MTS' % Meter–tonne–second
        unitSystem = {
            'm'     1
            't'     1/1000};
    case 'CGS' % Centimeter–gram–second
        unitSystem = {
            'cm'   100
            'g'    1000};
    case 'MMGS' % Milimeter–gram–second
        unitSystem = {
            'mm'   1000
            'g'    1000};
    case {'GM' 'GRAVITATIONAL METRIC'}
        unitSystem = {
            'm'     1
            'hyl'   1/9.80665};
    case {'MKH'}
        unitSystem = {
            'm'   1
            'kg'  1
            'hr'  1/3600};
    otherwise
        error('Unknown unit system.')
end

baseUnitSystem(1:size(unitSystem,1),:) = unitSystem;

%% %% CUSTOM EXAMPLES BELOW %%%%

%% No units.
% baseUnitSystem = 'none';

%% GCS with cents for currency, µA for current, ms for time.
% baseUnitSystem = {
%     'cm'    100     % 100   cm  / m
%     'g'     1000    % 1000  g   / kg
%     'ms'    1000
%     'µA'    1e6
%     'K'     1
%     'mol'   1
%     'cd'    1
%     'bit'   1
%     '¢'     100     % 100 cents per dollar
%     };