function u = units(unitSystem,varargin)
% u = UNITS  Returns a struct for use in myUnits.m. Call this function from a
% project's myUnits.m to define the behavior of the physical units class u.
% 
%   u = UNITS(baseUnitStr) uses the base units as indicated by baseUnitStr. Some
%   common possibilities are:
%     SI,MKS    - Standard SI/metric (default).
%     CGS       - cm-g-s.
%     MMGS      - mm-g-s.
%     FPS,EE,AE - ft-lb-s / English Engineering / Absolute English.
%     FSS,BG    - ft-slug-s / Gravitational FPS / Technical FPS / 
%                 British Gravitational.
%     Verbose   - Spelled-out base units (meter, kilogram, second, etc.).
% 
%   Uncommon: u = UNITS(baseUnitSystem) creates a units struct based on the base
%   units provided in the n-by-2 cell array baseUnitSystem. The first element in
%   each row of the cell array is the string for the base unit, and the second
%   element is the number that the unit must be multiplied by to get to the
%   corresponding SI base unit. The rows must be in order of length, mass, time,
%   current, absolute temperature, .... Not all rows must be provided if only,
%   e.g., non-SI length and mass are desired.
% 
%   u = UNITS(..., dispUnit1, dispUnit2, dispUnit2, ....) allows for displaying
%   units using preferred display strings corresponding to a supported unit.
%   u = UNITS(..., '-SI') is a special shorthand for using standard derived SI
%   units (N, Pa, J, etc.). 
% 
%   u = UNITS(..., dispUnits) allows to displaying units using strings,
%   characters, and symbols that are not valid MATLAB identifiers and/or don't
%   correspond to a pre-defined unit. In this case it is necessary to call UNITS
%   twice, with the first call used to generate building blocks necessary for
%   defining derived display units. dispUnits is an n-by-2 cell array with
%   strings in the first column and corresponding unit in the second.
% 
%   See default myUnits.m for many examples.
% 
%   See also myUnits.


%% ----- Set up the fundamental dimensions over which to assign units -----
% See also www.mathworks.com/help/physmod/simscape/ug/unit-definitions.html.

defaultUnitSystem = 'SI';

% Check input
if nargin < 1
    unitSystem = defaultUnitSystem;
end

% Everything is based off of SI:
siUnitSystem =  {'m'    1
                 'kg'   1
                 's'    1
                 'A'    1 
                 'K'    1
                 'mol'  1
                 'cd'   1
                 'USD'  1};

siDimensionNames = siUnitSystem(:,1);
siMultipliers = [siUnitSystem{:,2}];
nSIunits = length(siDimensionNames);

%% Pre-defined shorthand unit systems
if ischar(unitSystem)
    % Shorthand used for a predefined system - set unitSystem to
    % corresponding cell array. 
    switch upper(unitSystem)
        case {'SI' 'MKS' 'METRIC' 'INTERNATIONAL'}
            unitSystem = siUnitSystem;
        case 'IPS'
            unitSystem = {'in'  1/0.0254
                'lbm'   2.2046
                's'    1
                'A'    1
                'R'    1.8};
        case 'IPSK'
            unitSystem = {'in'  1/0.0254
                'lbm'   2.2046};
        case {'FPS' 'EE' 'AE' 'IMPERIAL' 'AMERICAN'} 
            % EE = English Engineering; AE = Absolute English
            unitSystem = {'ft'  1/0.3048
                'lbm'   2.2046 
                's'    1
                'A'    1
                'R'    1.8};
        case 'FPSK'
            unitSystem = {'ft'  1/0.3048
                'lbm'   2.2046};
        case {'FSS' 'BG' 'GRAVITATIONAL FPS' 'TECHNICAL FPS'} 
            % BG = British Gravitational
            unitSystem = {'ft'  1/0.3048
                'slug' 0.3048/4.4482216152605 %exact
                's'    1
                'A'    1
                'R'    1.8};
        case {'VERBOSE' 'VERBOSE MKS' 'VERBOSE SI' 'SI VERBOSE' ...
                'MKS VERBOSE'}
            unitSystem = {'meter'  1
                'kilogram'   1
                'second'    1
                'ampere'    1
                'kelvin'    1
                'mole'      1
                'candela'   1
                'dollar'    1};
        case {'VERBOSE FPS' 'FPS VERBOSE' 'VERBOSE EE' 'EE VERBOSE'}
            unitSystem = {'foot'  1/0.3048
                'poundMass'   2.2046
                'second'    1
                'ampere'    1
                'Rankine'    1.8
                'mole'      1
                'candela'   1
                'dollar'    1};
        case 'MTS' % Meter–tonne–second
            unitSystem = {'m'   1
                 't'    1/1000};
        case 'CGS' % Centimeter–gram–second
            unitSystem = {'cm'   100
                 'g'    1000};
        case 'MMGS' % Milimeter–gram–second
            unitSystem = {'mm'   1000
                 'g'    1000};
        case {'GM' 'GRAVITATIONAL METRIC'}
            unitSystem = {'m'   1
                 'hyl'  1/9.80665};
        case {'MKH'}
            unitSystem = {'m'   1
                 'kg'  1
                 'hr'  1/3600};
        otherwise
            error('Unknown unit system.')
    end
end

if iscell(unitSystem)
    providedNames = unitSystem(:,1);
    providedMultipliers = [unitSystem{:,2}];
    
    nProvidedUnits = length(providedNames);
else
    error('Base unit system must be short-hand string or cell array.')
end

%% Preferred display units
% Allow for special shorthand cases.
if nargin == 2 && ischar(varargin{1}) && strncmpi('-',varargin{1},1)
    switch upper(varargin{1})
        case '-SI'
            % Standard set of common derived SI units.
            varargin = {'N' 'Pa' 'J' 'W' 'C' 'V' 'F'...
                'Ohm' 'S' 'Wb' 'T' 'H'};
            
        case '-FPS'
            varargin = {'lbf' 'psf' 'hp' 'V' 'Ohm'};
            
        % Add shorthand cases for your particular application here.
    end
end

% Process the input.
if nargin == 2 && iscell(varargin{1})
    % Assume dispUnits is already in correct form and ready to go.
    dispUnits = varargin{1};
elseif nargin >=2
    if iscellstr(varargin)
        % Just a list of valid unit strings.
        u = units(unitSystem);
        nPrefs = length(varargin);
        dispUnits = varargin';
        for i = 1:nPrefs
            dispUnits{i,2} = u.(dispUnits{i});
        end
    else
        % Assume varargin is a list of str, DimVar pairs.
        nPrefs = length(varargin)/2;
        dispUnits = reshape(varargin,2,nPrefs)';
    end
else
    dispUnits = [];
end

% Check valid dispUnits.
if ~isempty(dispUnits)
    % Check that cell is m x 2.
    [~,m2] = size(dispUnits);
    if m2 ~= 2
        error('Provided dispUnits in a cell array must be n x 2.')
    end
    
    % Check for strings in the first column and DimVars in the second.
    if ~iscellstr(dispUnits(:,1))
        error('First column of dispUnits cell array must contain strings.')
    end
    if ~all(cellfun('isclass', dispUnits(:,2), 'DimVar'))
        warning(...
            'Second column of dispUnits cell array should be all DimVars.')
    end
end

%% Build DimVar.

if nSIunits > nProvidedUnits
    nDimensions = nSIunits;
    
    dimensionNames = siDimensionNames;
    dimensionNames(1:nProvidedUnits) = providedNames;
    
    multipliers = siMultipliers;
    multipliers(1:nProvidedUnits) = providedMultipliers;
else
    nDimensions = nProvidedUnits;
    dimensionNames = providedNames;
    multipliers = providedMultipliers;
end

clear class
for nd = 1:nDimensions
    u.(dimensionNames{nd}) = DimVar(dimensionNames, dimensionNames{nd},...
        dispUnits);
end

%------ Get back into SI units for defining all subsequent units.------
for i = 1:nSIunits
    u.(siUnitSystem{i,1}) = u.(dimensionNames{i})*multipliers(i);
end 

%% Derived units:
% This list should mirror that in u.m.
% These are only used for custom display of derived units.

%---- length ----

u.km = 1e3*u.m;
u.dm = 1e-1*u.m;
u.cm = 1e-2*u.m;
u.mm = 1e-3*u.m;
u.um = 1e-6*u.m;
u.micron = u.um;
u.nm = 1e-9*u.m;
u.ang = 1e-10*u.m;            % ångström
u.in = 2.54*u.cm;
u.mil = 1e-3*u.in;
u.ft = 12*u.in;
u.kft = 1e3*u.ft;
u.yd = 3*u.ft;
u.mi = 5280*u.ft;
u.nmi = 1852*u.m;             % nautical mile
u.NM = u.nmi;                 % nautical mile
u.a0 = .529e-10*u.m;          % Bohr radius

%---- area ----

u.ha = 10000*u.m^2;
u.hectare = u.ha;
u.ac = 43560*u.ft^2;
u.acre = u.ac;

%---- volume ----

u.cc = (u.cm)^3;
u.L = 1000*u.cc;
u.mL = u.cc;
u.cuin = 16.387064*u.mL;      % cubic inch
u.gal = 231*u.cuin;           % US gallon
u.quart = u.gal/4;
u.pint = u.quart/2;
u.cup = u.pint/2;
u.floz = u.cup/8;             % fluid ounce
u.Tbls = u.floz/2;
u.tsp = u.Tbls/3;

%---- acceleration ----

u.g0 = 9.80665*u.m/u.s^2;     % standard gravity
u.gn = u.g0;                  % standard gravity
u.Gal = u.cm/u.s^2;           % en.wikipedia.org/wiki/Gal_(unit)

%---- force ----

u.N = u.kg*u.m/u.s^2;
u.kN = 1000*u.N;
u.dyn = 1e-5*u.N;
u.lbf = 4.4482216152605*u.N;  % pound force
u.kgf = u.kg*u.g0;            % kilogram force
u.kp = u.kgf;                 % kilopond
u.p = u.kp/1000;              % pond
u.sn = u.kN;                  % Sthene

%---- mass ----

u.gram = 1e-3*u.kg;
u.g = u.gram;                 % gram
u.mg = 1e-3*u.gram;
u.lbm = 0.45359237*u.kg;      % pound mass
u.lb = u.lbm;                 % pound mass
u.st = 14*u.lbm;              % stone
u.stone = u.st;
u.slug = u.lbf/(u.ft/u.s^2);
u.oz = (1/16)*u.lbm;
u.amu = 1.660539040e-27*u.kg; % atomic mass unit
u.Da = u.amu;                 % atomic mass unit
u.t = 1000*u.kg;
u.tonne = u.t;
u.mug = u.kgf/(u.m/u.s^2);    % metric slug
u.hyl = u.mug;
u.TMU = u.mug;                % technische Masseneinheit

%---- more force ----

u.pdl = u.lbm*u.ft/u.s^2;     % poundal
u.gramForce = u.gram*u.g0;
u.gf = u.gramForce;           % gram force
u.ozf = u.oz*u.g0;            % ounce force

%---- time ----

u.ms = 1e-3*u.s;
u.us = 1e-6*u.s;
u.ns = 1e-9*u.s;
u.ps = 1e-12*u.s;
u.min = 60*u.s;
u.hr = 60*u.min;
u.day = 24*u.hr;
u.week = 7*u.day;
u.fortnight = 2*u.week;
u.year = 365.25*u.day;        % Julian year; defines light-year
u.month = u.year/12;          % 1/12th Julian year

%---- frequency ----

% hertz; note: Incompatible with angle and angular velocity units.
u.Hz = 1/u.s;
u.kHz = 1e3 *u.Hz;
u.MHz = 1e6 *u.Hz;
u.GHz = 1e9 *u.Hz;

%---- energy ----

u.J = u.N*u.m;
u.MJ = 1e6*u.J;
u.kJ = 1e3*u.J;
u.mJ = 1e-3*u.J;
u.uJ = 1e-6*u.J;
u.nJ = 1e-9*u.J;
u.eV = 1.6022e-19*u.J;
u.BTU = 1.0550559e3*u.J;
u.kWh = 3.6e6*u.J;            % kilowatt-hour
u.Wh = 3.6e3*u.J;             % watt-hour
u.cal = 4.1868*u.J;
u.kCal = 1e3*u.cal;
u.erg = 1e-7*u.J;             % en.wikipedia.org/wiki/Erg
u.quad = 1e15*u.BTU;          % en.wikipedia.org/wiki/Quad_(unit)

%---- temperature ----
% For reference: °C = °K-273.15; °F = °R-459.67.

u.R = u.K*5/9;
u.mK = 1e-3*u.K;
u.uK = 1e-6*u.K;
u.nK = 1e-9*u.K;

%---- pressure ----

u.Pa = u.N/u.m^2;
u.mPa = u.Pa/1000;
u.kPa = u.Pa*1000;
u.MPa = u.kPa*1000;
u.torr = 133.322*u.Pa;
u.mtorr = 1e-3*u.torr;
u.bar = 1e5*u.Pa;
u.mbar = 1e-3*u.bar;
u.atm = 101325*u.Pa;
u.psi = u.lbf/u.in^2;
u.ksi = 1000*u.psi;
u.Msi = 1000*u.ksi;
u.psf = u.lbf/u.ft^2;
u.Ba = 0.1*u.Pa;              % Barye
u.pz = u.kPa;                 % pièze
u.mmHg = 133.322387415*u.Pa;
u.inHg = 25.4*u.mmHg;

%---- viscosity ----

u.St = u.cm^2/u.s;            % stokes (kinematic viscosity)
u.cSt = u.St/100;
u.P = u.Pa * u.s / 10;        % poise (dynamic viscosity)
u.cP = u.mPa * u.s;

%---- power ----

u.W = u.J/u.s;
u.MW = 1e6*u.W;
u.kW = 1e3*u.W;
u.mW = 1e-3*u.W;
u.uW = 1e-6*u.W;
u.nW = 1e-9*u.W;
u.pW = 1e-12*u.W;
u.hp = 550*u.ft*u.lbf/u.s;    % mechanical horsepower
u.hpE = 746*u.W;              % electrical horsepower
u.PS = 75*u.kg*u.g0*u.m/u.s;  % metric horsepower (DIN 66036 definition)

%---- current ----

u.mA = 1e-3*u.A;
u.uA = 1e-6*u.A;
u.nA = 1e-9*u.A;

%---- charge ----

u.C = u.A*u.s;                % coulomb
u.e = 1.6022e-19*u.C;         % elementary charge
u.mC = 1e-3*u.C;
u.uC = 1e-6*u.C;
u.nC = 1e-9*u.C;
u.pC = 1e-12*u.C;

u.mAh = u.mA*u.hr;            % milliamp-hour
u.Ah = u.A*u.hr;              % amp-hour

%---- voltage ----

u.V = 1*u.J/u.C;              % volt
u.kV = 1e3*u.V;
u.mV = 1e-3*u.V;
u.uV = 1e-6*u.V;

%---- resistance ----

u.Ohm = u.V/u.A;
u.MOhm = 1e6*u.Ohm;
u.kOhm = 1e3*u.Ohm;
u.mOhm = 1e-3*u.Ohm;
u.S = 1/u.Ohm;                % siemens

%---- capacitance ----

u.F = u.A*u.s/u.V;            % farad
u.mF = 1e-3*u.F;
u.uF = 1e-6*u.F;
u.nF = 1e-9*u.F;
u.pF = 1e-12*u.F;

%---- inductance ----

u.H = u.Ohm*u.s;              % henry
u.mH = 1e-3*u.H;

%---- EM ----

u.T = 1*u.N/(u.A*u.m);        % tesla
u.gauss = 1e-4*u.T;
u.Wb = u.V*u.s;               % weber
u.mWb = u.Wb/1000;
u.uWb = 1e-6*u.Wb;
u.nWb = 1e-9*u.Wb;

%---- fundamental constants ----
% See http://www.efunda.com/units/show_constants.cfm for more.

u.kB = 1.38e-23*u.J/u.K;              % Boltzmann constant
u.sigma_SB = 5.670e-8 * u.W/(u.m^2 * u.K^4); % Stefan–Boltzmann constant
u.h = 6.62607004e-34 * u.J*u.s;       % Planck constant
u.hbar = u.h/(2*pi);                  % Dirac constant
u.mu_B = 9.27400999e-24 * u.J/u.T;    % Bohr magneton
u.mu_N = 5.050783699e-27 * u.J/u.T;   % nuclear magneton
u.c = 299792458 * u.m/u.s;            % speed of light in vacuum
u.eps0 = 8.8541878176204e-12* u.C/(u.V*u.m); % vacuum permittivity
u.mu0 = 1.2566370614359e-6 * u.J/(u.m*u.A^2); % vacuum permeability

% specific gas constant for air (ESDU 77022 definition)
u.Rair = 287.05287*u.J/u.kg/u.K;

u.ly = u.c*u.year;            % light-year
u.lightYear = u.ly;

%---- non-dimensionals ----

u.percent = 0.01;             % %
u.pct = u.percent;
u.permil = 0.001;             % ‰
u.permill = u.permil;
u.permille = u.permil;
u.permyriad = 1e-4;           % ?
u.bp = u.permyriad;           % basis point
u.ppm = 1e-6;                 % part per million
u.ppb = 1e-9;                 % part per billion
u.ppt = 1e-12;                % part per trillion
u.ppq = 1e-15;                % part per quadrillion

%---- angles ----
% Note: angles are dimensionless

u.rad = 1;                    % radian
u.sr = 1;                     % steradian
u.turn = 2*pi*u.rad;
u.rev = u.turn;               % revolution = 2*pi radians
u.deg = u.turn/360;           % degree
u.arcminute = u.deg/60;
u.arcsecond = u.arcminute/60;
u.grad = u.turn/400;          % gradian

%---- rotational speeds ----

u.rpm = u.rev/u.min;

%---- speeds ----

u.mps = u.m/u.s;
u.fps = u.ft/u.s;
u.kt = u.nmi/u.hr;
u.kn = u.kt;
u.kts = u.kt;
u.knot = u.kt;
u.KTAS = u.kt;
u.nmph = u.kt;
u.kph = u.km/u.hr;
u.mph = u.mi/u.hr;
u.fpm = u.ft/u.min;

%---- volume flow rates ----

u.cfm = u.ft^3/u.min;         % cubic feet per minute
u.cfs = u.ft^3/u.s;           % cubic feet per second

%---- other derived SI ----

u.kat = u.mol/u.s;            % katal
u.lm = u.cd*u.sr;             % lumen
u.lx = u.lm/u.m^2;            % lux

%---- currency ----
% See also mathworks.com/matlabcentral/fileexchange/47255

u.dollar = u.USD;
u.cent = u.USD/100;

end

%   Original seed for this class by Rob deCarvalho.
%     http://www.mathworks.com/matlabcentral/fileexchange/authors/22148
%     http://www.mathworks.com/matlabcentral/fileexchange/10070