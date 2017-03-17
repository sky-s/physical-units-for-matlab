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

% References:
% http://physics.nist.gov/cuu/Constants/index.html
% http://www.translatorscafe.com/unit-converter
% http://en.wikipedia.org
% http://www.efunda.com/units/index.cfm


%% ----- Set up the fundamental dimensions over which to assign units -----
% See also www.mathworks.com/help/physmod/simscape/ug/unit-definitions.html,
% www.mathworks.com/help/symbolic/units-list.html.

defaultUnitSystem = 'SI';

% Check input
if nargin < 1
    unitSystem = defaultUnitSystem;
end

% Everything is based off of SI:
siUnitSystem =  {
    'm'        1
    'kg'       1
    's'        1
    'A'        1
    'K'        1
    'mol'      1
    'cd'       1
    'currency' 1
    };

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
            unitSystem = {
                'in'  1/0.0254
                'lbm'   2.2046
                's'    1
                'A'    1
                'R'    1.8};
        case 'IPSK'
            unitSystem = {
                'in'  1/0.0254
                'lbm'   2.2046};
        case {'FPS' 'EE' 'AE' 'IMPERIAL' 'AMERICAN'} 
            % EE = English Engineering; AE = Absolute English
            unitSystem = {
                'ft'  1/0.3048
                'lbm'   2.2046 
                's'    1
                'A'    1
                'R'    1.8};
        case 'FPSK'
            unitSystem = {
                'ft'  1/0.3048
                'lbm'   2.2046};
        case {'FSS' 'BG' 'GRAVITATIONAL FPS' 'TECHNICAL FPS'} 
            % BG = British Gravitational
            unitSystem = {
                'ft'  1/0.3048
                'slug' 0.3048/4.4482216152605 %exact
                's'    1
                'A'    1
                'R'    1.8};
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
                'currency'  1};
        case {'VERBOSE FPS' 'FPS VERBOSE' 'VERBOSE EE' 'EE VERBOSE'}
            unitSystem = {
                'foot'  1/0.3048
                'poundMass' 2.2046
                'second'    1
                'ampere'    1
                'Rankine'   1.8
                'mole'      1
                'candela'   1
                'currency'  1};
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
            varargin = {'N' 'Pa' 'J' 'W' 'C' 'V' 'F' 'Ohm' 'S' 'Wb' 'T' 'H'};
            
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
% This list should mirror that in the classdef u.m.
% These are only used for custom display of derived units.

%---- length ----

u.km = 1e3*u.m;               % kilometer
u.dm = 1e-1*u.m;              % decimeter
u.cm = 1e-2*u.m;              % centimeter
u.mm = 1e-3*u.m;              % millimeter
u.um = 1e-6*u.m;              % micrometer
u.micron = u.um;              % micron
u.nm = 1e-9*u.m;              % nanometer
u.pm = 1e-12*u.m;             % picometer
u.fm = 1e-15*u.m;             % femtometer
u.fermi = u.fm;               % fermi
u.Ao = 1e-10*u.m;             % ångström
u.ang = u.Ao;                 % ångström
u.a0 = 0.52917721067e-10*u.m; % Bohr radius
u.a_0 = u.a0;                 % Bohr radius
u.lP = 1.616229e-35*u.m;      % Planck length
u.xu = 1.0021e-13*u.m;        % x unit
u.xu_Cu = 1.00207697e-13*u.m; % x unit (copper)
u.xu_Mo = 1.00209952e-13*u.m; % x unit (molybdenum)
u.in = 2.54*u.cm;             % inch
u.mil = 1e-3*u.in;            % mil
u.line = u.in/10;             % line
u.hand = 4*u.in;              % hand
u.span = 9*u.in;              % span
u.smoot = 67*u.in;            % smoot
u.ft = 12*u.in;               % foot
u.ft_US = 1200/3937*u.m;      % US survey foot
u.kft = 1e3*u.ft;             % kilofoot
u.yd = 3*u.ft;                % yard
u.fathom = 6*u.ft;            % fathom
u.ftm = u.fathom;             % fathom
u.li = 0.66*u.ft;             % link
u.rod = 5.5*u.yd;             % rod
u.ch = 66*u.ft;               % chain
u.furlong = 220*u.yd;         % furlong
u.fur = u.furlong;            % furlong
u.mi = 5280*u.ft;             % mile
u.mi_US = 6336/3937*u.km;     % US survey mile
u.nmi = 1852*u.m;             % nautical mile
u.NM = u.nmi;                 % nautical mile
u.inm = u.nmi;                % nautical mile
u.nm_UK = 6080*u.ft;          % Imperial nautical mile
u.nmile = u.nm_UK;            % Imperial nautical mile
u.au = 149597870.7*u.km;      % astronomical unit
u.pc = 648000/pi*u.au;        % parsec

%---- reciprocal length ----

u.dpt = 1/u.m;                % diopter
u.R_inf = 1.0973731568508e7/u.m; % Rydberg constant

%---- area ----

u.square = 100*u.ft^2;        % square
u.ha = 10000*u.m^2;           % hectare
u.hectare = u.ha;             % hectare
u.a = 100*u.m^2;              % are
u.are = u.a;                  % are
u.ac = 43560*u.ft^2;          % acre
u.acre = u.ac;                % acre
u.ro = 1/4*u.acre;            % rood
u.twp = 36*u.mi^2;            % township
u.circ_mil = pi/4*u.mil^2;    % circular mil
u.circ_inch = pi/4*u.in^2;    % circular inch
u.b = 100*u.fm^2;             % barn
u.barn = u.b;                 % barn

%---- volume ----

u.cc = u.cm^3;                % cubic centimeter
u.L = 1000*u.cc;              % liter
u.l = u.L;                    % liter
u.dL = 100*u.cc;              % deciliter
u.dl = u.dL;                  % deciliter
u.cL = 10*u.cc;               % centiliter
u.cl = u.cL;                  % centiliter
u.mL = u.cc;                  % milliliter
u.ml = u.mL;                  % milliliter
u.uL = u.mm^3;                % microliter
u.ul = u.uL;                  % microliter
u.kL = u.m^3;                 % kiloliter
u.kl = u.kL;                  % kiloliter
u.cuin = 16.387064*u.mL;      % cubic inch
u.FBM = u.ft^2*u.in;          % board foot
u.gal = 231*u.cuin;           % US gallon
u.gal_UK = 4.54609*u.l;       % UK imperial gallon
u.igal = u.gal_UK;            % UK imperial gallon
u.quart = u.gal/4;            % US quart
u.qt_UK = u.gal_UK/4;         % British imperial quart
u.liq_qt = u.quart;           % US quart
u.pint = u.quart/2;           % US pint
u.pint_UK = u.qt_UK/2;        % British imperial pint
u.liq_pt = u.pint;            % US pint
u.cup = u.pint/2;             % US cup
u.floz = u.cup/8;             % US fluid ounce
u.floz_UK = u.gal_UK/160;     % British imperial fluid ounce
u.Tbls = u.floz/2;            % US tablespoon
u.tsp = u.Tbls/3;             % US teaspoon
u.acft = u.acre*u.ft;         % acre-foot
u.acin = u.acre*u.in;         % acre-inch
u.barrel = 7056*u.in^3;       % US customary dry barrel
u.bbl = u.barrel;             % US customary dry barrel
u.fldr = u.floz/8;            % US customary fluid dram
u.fldr_UK = u.floz_UK/8;      % British imperial fluid drachm (dram)
u.minim = u.fldr/60;          % US customary minim
u.minim_UK = u.fldr_UK/60;    % British imperial minim
u.gill = 4*u.floz;            % US customary fluid gill
u.gill_UK = u.gal_UK/32;      % British imperial gill

%---- acceleration ----

u.g0 = 9.80665*u.m/u.s^2;     % standard gravity
u.gn = u.g0;                  % standard gravity
u.g_n = u.g0;                 % standard gravity
u.gee = u.g0;                 % standard gravity
u.Gal = u.cm/u.s^2;           % gal

%---- force ----

u.N = u.kg*u.m/u.s^2;         % newton
u.kN = 1000*u.N;              % kilonewton
u.mN = 1e-3*u.N;              % millinewton
u.dyn = 1e-5*u.N;             % dyne
u.lbf = 4.4482216152605*u.N;  % pound force
u.kip = 1000*u.lbf;           % kip
u.kgf = u.kg*u.g0;            % kilogram force
u.kp = u.kgf;                 % kilopond
u.p = u.kp/1000;              % pond
u.sn = u.kN;                  % sthène

%---- mass ----

u.gram = 1e-3*u.kg;           % gram
u.g = u.gram;                 % gram
u.mg = 1e-3*u.gram;           % milligram
u.ug = 1e-6*u.gram;           % microgram
u.Mg = 1e6*u.gram;            % Megagram/metric ton
u.t = 1000*u.kg;              % metric ton
u.tonne = u.t;                % metric ton
u.Mt = 1e6*u.t;               % metric megaton
u.lbm = 0.45359237*u.kg;      % pound mass
u.lb = u.lbm;                 % pound mass
u.tn = 2000*u.lbm;            % US customary short ton
u.ton_UK = 2240*u.lbm;        % British imperial ton
u.st = 14*u.lbm;              % stone
u.stone = u.st;               % stone
u.cwt = 100*u.lbm;            % US customary short hundredweight
u.cwt_UK = 8*u.stone;         % British imperial short hundredweight
u.quarter = u.cwt_UK/4;       % British imperial quarter
u.slug = u.lbf/(u.ft/u.s^2);  % slug
u.oz = u.lbm/16;              % ounce
u.dr = u.oz/16;               % dram
u.gr = u.lbm/7000;            % grain
u.ct = 200*u.mg;              % carat
u.amu = 1.660539040e-27*u.kg; % atomic mass unit
u.Da = u.amu;                 % atomic mass unit
u.mu = u.amu;                 % atomic mass unit
u.mP = 2.176470e-8*u.kg;      % Planck mass
u.m_e = 9.10938356e-31*u.kg;  % electron mass
u.mug = u.kgf/(u.m/u.s^2);    % metric slug
u.hyl = u.mug;                % hyl
u.TMU = u.mug;                % technische Masseneinheit

%---- more force ----

u.pdl = u.lbm*u.ft/u.s^2;     % poundal
u.gramForce = u.gram*u.g0;    % gram force
u.gf = u.gramForce;           % gram force
u.ozf = u.oz*u.g0;            % ounce force
u.tonf = u.tn*u.g0;           % short ton force

%---- mass per length ----

u.den = u.gram/(9*u.km);      % denier
u.tex = u.gram/u.km;          % tex
u.dtex = u.tex/10;            % decitex

%---- time ----

u.ms = 1e-3*u.s;              % millisecond
u.us = 1e-6*u.s;              % microsecond
u.ns = 1e-9*u.s;              % nanosecond
u.ps = 1e-12*u.s;             % picosecond
u.tP = 5.39116e-44*u.s;       % Planck time
u.min = 60*u.s;               % minute
u.h = 60*u.min;               % hour
u.hr = u.h;                   % hour
u.d = 24*u.hr;                % day
u.day = u.d;                  % day
u.week = 7*u.day;             % week
u.fortnight = 2*u.week;       % fortnight
u.month_30 = 30*u.day;        % 30-day month
u.yr = 365.25*u.day;          % julian year
u.y = u.yr;                   % julian year
u.year = u.yr;                % julian year
u.year_julian = u.year;       % julian year
u.year_360 = 360*u.day;       % 360-day year
u.year_Tropical = 365.24219*u.day; % tropical year
u.year_Gregorian = 365.2425*u.day; % gregorian year
u.month = u.yr/12;            % 1/12th julian year

%---- frequency ----

u.Hz = 1/u.s;   % hertz (NB: incompatible with angle and angular velocity)
u.kHz = 1e3*u.Hz;             % kilohertz
u.MHz = 1e6*u.Hz;             % megahertz
u.GHz = 1e9*u.Hz;             % gigahertz

%---- energy ----

u.Nm = u.N*u.m;               % newton meter
u.J = u.Nm;                   % joule
u.MJ = 1e6*u.J;               % megajoule
u.kJ = 1e3*u.J;               % kilojoule
u.mJ = 1e-3*u.J;              % millijoule
u.uJ = 1e-6*u.J;              % microjoule
u.nJ = 1e-9*u.J;              % nanojoule
u.eV = 1.6021766208e-19*u.J;  % electronvolt
u.BTU = 1055.06*u.J;          % British thermal unit (ISO)
u.Btu = u.BTU;                % British thermal unit (ISO)
u.Btu_IT = 1055.0559*u.J;     % British thermal unit (International Table)
u.Btu_th = 1054.3503*u.J;     % British thermal unit (thermochemical)
u.kpm = u.kp*u.m;             % kilopond meter
u.Ws = u.J;                   % watt-second
u.kWh = 3.6e6*u.J;            % kilowatt-hour
u.Wh = 3.6e3*u.J;             % watt-hour
u.cal = 4.1868*u.J;           % calorie (International Table)
u.cal_IT = u.cal;             % calorie (International Table)
u.cal_4 = 4.204*u.J;          % calorie (4°C)
u.cal_15 = 4.1855*u.J;        % calorie (15°C)
u.cal_20 = 4.182*u.J;         % calorie (20°C)
u.cal_mean = 4.190*u.J; 	  % calorie (mean)
u.cal_th = 4.184*u.J;         % calorie (thermochemical)
u.kcal = 1e3*u.cal;           % kilocalorie
u.kcal_IT = 1e3*u.cal_IT;     % kilocalorie (International Table)
u.Cal = u.kcal;               % large calorie / food calorie
u.kcal_4 = 1e3*u.cal_4;       % kilocalorie (4°C)
u.kcal_15 = 1e3*u.cal_15;     % kilocalorie (15°C)
u.kcal_20 = 1e3*u.cal_20;     % kilocalorie (20°C)
u.kcal_mean = 1e3*u.cal_mean; % kilocalorie (mean)
u.kcal_th = 1e3*u.cal_th;     % kilocalorie (thermochemical)
u.erg = 1e-7*u.J;             % en.wikipedia.org/wiki/Erg
u.E_h = 4.359744650e-18*u.J;  % hartree energy
u.Ha = u.E_h;                 % hartree
u.thm = 1e5*u.BTU;            % therm
u.therm = u.thm;              % therm
u.quad = 1e15*u.BTU;          % quad

%---- temperature ----
% For reference: °C = °K-273.15; °F = °R-459.67.

u.R = u.K*5/9;                % rankine (°F = °R-459.67)
u.mK = 1e-3*u.K;              % millikelvin
u.uK = 1e-6*u.K;              % microkelvin
u.nK = 1e-9*u.K;              % nanokelvin
u.deltaK = u.K;               % kelvin (relative temperature)
u.deltadegC = u.K;            % celsius (relative, °C = °K-273.15)
u.deltadegR = u.R;            % rankine (relative temperature)
u.deltadegF = u.R;            % fahrenheit (relative, °F = °R-459.67)
u.TP = 1.416808e32*u.K;       % Planck temperature

%---- pressure ----

u.Pa = u.N/u.m^2;             % pascal
u.mPa = 1e-3*u.Pa;            % millipascal
u.kPa = 1e3*u.Pa;             % kilopascal
u.MPa = 1e6*u.Pa;             % megapascal
u.GPa = 1e9*u.Pa;             % gigapascal
u.torr = 133.322*u.Pa;        % torr
u.Torr = u.torr;              % torr
u.mtorr = 1e-3*u.torr;        % millitorr
u.bar = 1e5*u.Pa;             % bar
u.mbar = 1e-3*u.bar;          % millibar
u.kbar = 1e3*u.bar;           % kilobar
u.atm = 101325*u.Pa;          % standard atmosphere
u.at = u.kgf/u.cm^2;          % technical atmosphere
u.psi = u.lbf/u.in^2;         % pound force per square inch
u.ksi = 1e3*u.psi;            % kip per square inch
u.Msi = 1e6*u.psi;            % million pound force per square inch
u.psf = u.lbf/u.ft^2;         % pound force per square foot
u.ksf = u.kip/u.ft^2;         % kip per square foot
u.Ba = 0.1*u.Pa;              % barye
u.pz = u.kPa;                 % pièze
u.mmHg = 13.5951*u.kgf/u.m^2; % millimeter of mercury
u.cmHg = 10*u.mmHg;           % centimeter of mercury
u.mHg = 1e3*u.mmHg;           % meter of mercury
u.inHg = 2.54*u.cmHg;         % inch of mercury
u.ftHg = 12*u.inHg;           % foot of mercury
u.mmH20 = u.kgf/u.m^2;        % millimeter of water (density 1 g/cc)
u.mmAq = u.mmH20;             % millimeter of water
u.cmH20 = 10*u.mmH20;         % centimeter of water
u.cmAq = u.cmH20;             % centimeter of water
u.mH20 = 1e3*u.mmH20;         % meter of water
u.mAq = u.mH20;               % meter of water
u.inH20 = 2.54*u.cmH20;       % inch of water
u.inAq = u.inH20;             % inch of water
u.wc = u.inH20;               % inch water column
u.ftH20 = 12*u.inH20;         % foot of water
u.ftAq = u.ftH20;             % foot of water

%---- viscosity ----

u.St = u.cm^2/u.s;            % stokes
u.cSt = u.St/100;             % centistokes
u.newt = u.in^2/u.s;          % newt
u.P = u.Pa*u.s / 10;          % poise
u.cP = u.mPa*u.s;             % centipoise
u.reyn = u.lbf*u.s/u.in^2;    % reyn

%---- power ----

u.W = u.J/u.s;                % watt
u.MW = 1e6*u.W;               % megawatt
u.GW = 1e9*u.W;               % gigawatt
u.kW = 1e3*u.W;               % kilowatt
u.mW = 1e-3*u.W;              % milliwatt
u.uW = 1e-6*u.W;              % microwatt
u.nW = 1e-9*u.W;              % nanowatt
u.pW = 1e-12*u.W;             % picowatt
u.hp = 550*u.ft*u.lbf/u.s;    % mechanical horsepower (550 ft-lbf/s)
u.HP_I = u.hp;                % mechanical horsepower (550 ft-lbf/s)
u.hpE = 746*u.W;              % electrical horsepower
u.HP_E = u.hpE;               % electrical horsepower
u.PS = 75*u.kg*u.g0*u.m/u.s;  % metric horsepower (DIN 66036)
u.HP = u.PS;                  % metric horsepower (DIN 66036)
u.HP_DIN = u.PS;              % metric horsepower (DIN 66036)

%---- current ----

u.mA = 1e-3*u.A;              % milliampere
u.uA = 1e-6*u.A;              % microampere
u.nA = 1e-9*u.A;              % nanoampere
u.pA = 1e-12*u.A;             % picoampere
u.kA = 1e3*u.A;               % kiloampere
u.abA = 10*u.A;               % abampere
u.Bi = u.abA;                 % biot

%---- charge ----

u.C = u.A*u.s;                % coulomb
u.e = 1.6021766208e-19*u.C;   % elementary charge
u.mC = 1e-3*u.C;              % millicoulomb
u.uC = 1e-6*u.C;              % microcoulomb
u.nC = 1e-9*u.C;              % nanocoulomb
u.pC = 1e-12*u.C;             % picocoulomb
u.abC = 10*u.C;               % abcoulomb
u.aC = u.abC;                 % abcoulomb
u.statC = u.dyn^(1/2)*u.cm;   % statcoulomb
u.Fr = u.statC;               % franklin
u.esu = u.statC;              % electrostatic unit of charge
u.mAh = u.mA*u.hr;            % milliamp-hour
u.Ah = u.A*u.hr;              % amp-hour

%---- voltage ----

u.V = 1*u.J/u.C;              % volt
u.kV = 1e3*u.V;               % kilovolt
u.MV = 1e6*u.V;               % megavolt
u.GV = 1e9*u.V;               % gigavolt
u.mV = 1e-3*u.V;              % millivolt
u.uV = 1e-6*u.V;              % microvolt

%---- resistance/conductance ----

u.Ohm = u.V/u.A;              % ohm
u.GOhm = 1e9*u.Ohm;           % gigaohm
u.MOhm = 1e6*u.Ohm;           % megaohm
u.kOhm = 1e3*u.Ohm;           % kiloohm
u.mOhm = 1e-3*u.Ohm;          % milliohm
u.uOhm = 1e-6*u.Ohm;          % microohm
u.nOhm = 1e-9*u.Ohm;          % nanoohm
u.abOhm = u.nOhm;             % abohm
u.Z0 = 376.730313461*u.Ohm;   % characteristic impedance of vacuum
u.R_K = 25812.8074555*u.Ohm;  % von Klitzing constant
u.R_K_90 = 25812.807*u.Ohm;   % von Klitzing constant (conventional value)
u.S = 1/u.Ohm;                % siemens
u.mS = 1e-3*u.S;              % millisiemens
u.uS = 1e-6*u.S;              % microsiemens
u.nS = 1e-9*u.S;              % nanosiemens
u.G0 = 7.7480917310e-5*u.S;   % conductance quantum 


%---- capacitance ----

u.F = u.A*u.s/u.V;            % farad
u.mF = 1e-3*u.F;              % millifarad
u.uF = 1e-6*u.F;              % microfarad
u.nF = 1e-9*u.F;              % nanofarad
u.pF = 1e-12*u.F;             % picofarad

%---- inductance ----

u.H = u.Ohm*u.s;              % henry
u.mH = 1e-3*u.H;              % millihenry
u.uH = 1e-6*u.H;              % microhenry
u.nH = 1e-9*u.H;              % nanohenry
u.abH = u.nH;                 % abhenry
u.kH = 1e3*u.H;               % kilohenry
u.MH = 1e6*u.H;               % megahenry
u.GH = 1e9*u.H;               % gigahenry

%---- EM ----

u.T = 1*u.N/(u.A*u.m);        % tesla
u.Gs = 1e-4*u.T;              % gauss
u.Wb = u.V*u.s;               % weber
u.Mx = 1e-8*u.Wb;             % maxwell
u.mWb = u.Wb/1000;            % milliweber
u.uWb = 1e-6*u.Wb;            % microweber
u.nWb = 1e-9*u.Wb;            % nanoweber
u.Oe = 250/pi*u.A/u.m;        % oersted
u.Gb = 2.5/pi*u.A;            % gilbert

%---- non-dimensionals ----

u.percent = 0.01;             % %
u.pct = u.percent;            % %
u.permil = 0.001;             % ‰
u.permill = u.permil;         % ‰
u.permille = u.permil;        % ‰
u.permyriad = 1e-4;           % permyriad
u.bp = u.permyriad;           % basis point
u.ppm = 1e-6;                 % part per million
u.ppb = 1e-9;                 % part per billion
u.ppt = 1e-12;                % part per trillion
u.ppq = 1e-15;                % part per quadrillion

%---- angles ----
% Note: angles are dimensionless

u.rad = 1;                    % radian
u.sr = 1;                     % steradian
u.turn = 2*pi*u.rad;          % turn
u.rev = u.turn;               % revolution = 2*pi radians
u.deg = u.turn/360;           % degree
u.arcminute = u.deg/60;       % arcminute
u.arcmin = u.arcminute;       % arcminute
u.arcsecond = u.arcminute/60; % arcsecond
u.arcsec = u.arcsecond;       % arcsecond
u.grad = u.turn/400;          % gradian

%---- rotational speed ----

u.rpm = u.rev/u.min;          % revolution per minute
u.rps = u.rev/u.s;            % revolution per second

%---- velocity ----

u.mps = u.m/u.s;              % meter per second
u.kyne = u.cm/u.s;            % kyne
u.Kyne = u.kyne;              % kyne
u.fps = u.ft/u.s;             % foot per second
u.fpm = u.ft/u.min;           % foot per minute
u.kt = u.nmi/u.hr;            % knot
u.kn = u.kt;                  % knot
u.kts = u.kt;                 % knot
u.knot = u.kt;                % knot
u.knot_UK = u.nm_UK/u.hr;     % British imperial knot
u.KTAS = u.kt;                % knot
u.nmph = u.kt;                % knot
u.kph = u.km/u.hr;            % kilometer per hour
u.kmh = u.kph;                % kilometer per hour
u.mph = u.mi/u.hr;            % mile per hour

%---- volume flow rate ----

u.cfm = u.ft^3/u.min;         % cubic foot per minute
u.cfs = u.ft^3/u.s;           % cubic foot per second
u.gpm = u.gal/u.min;          % US customary gallon per minute
u.gpm_UK = u.gal_UK/u.min;    % British imperial gallon per minute
u.lpm = u.l/u.min;            % liter per minute

%---- fuel economy ----

u.l_100km = u.l/(100*u.km);   % liter per 100 km
u.mpg = u.mi/u.gal;           % mile per gallon

%---- Luminance etc. ----

u.asb = u.cd/u.m^2;           % apostilb
u.sb = u.cd/u.cm^2;           % stilb
u.ph = 1e4*u.cd*u.sr/u.m^2;   % phot
u.cp = 0.981*u.cd;            % candlepower
u.lm = u.cd*u.sr;             % lumen
u.lx = u.lm/u.m^2;            % lux
u.nx = 1e-3*u.lx;             % nox

%---- other derived SI ----

u.kat = u.mol/u.s;            % katal
u.M = u.mol/u.L;              % molar
u.molarity = u.M;             % molarity
u.Nms = u.N*u.m*u.s;          % newton-meter-second

%---- radiation ----

u.Gy = u.J/u.kg;              % gray
u.Sv = u.J/u.kg;              % sievert
u.Rad =  u.Gy/100;            % absorbed radiation dose
u.rem = u.Sv/100;             % roentgen equivalent man
u.roentgen = 2.58e-4*u.C/u.kg;% roentgen
u.Ly = u.cal_th/u.cm^2;       % langley
u.lan = u.Ly;                 % langley
u.Bq = 1/u.s;                 % becquerel
u.Ci = 3.7e10*u.Bq;           % curie

%---- constants ----

u.kB = 1.38064852e-23*u.J/u.K;        % Boltzmann constant
u.k_B = u.kB;                         % Boltzmann constant
u.sigma_SB = 5.670367e-8*u.W/(u.m^2*u.K^4); % Stefan–Boltzmann constant
u.h_c = 6.626070040e-34*u.J*u.s;      % Planck constant
u.h_bar = u.h/(2*pi);                 % Dirac constant
u.mu_B = 9.274009994e-24*u.J/u.T;     % Bohr magneton
u.mu_N = 5.050783699e-27*u.J/u.T;     % nuclear magneton
u.c = 299792458*u.m/u.s;              % speed of light in vacuum
u.c_0 = u.c;                          % speed of light in vacuum
u.ly = u.c*u.year;                    % light-year
u.lightYear = u.ly;                   % light-year
u.mu0 = pi*4e-7*u.N/u.A^2;            % vacuum permeability
u.eps0 = u.c^-2/u.mu0;                % vacuum permittivity
u.G = 6.67408e-11*u.m^3/u.kg/u.s^2;   % gravitational constant
u.N_A = 6.022140857e23/u.mol;         % Avogadro constant
u.NA = u.N_A;                         % Avogadro constant
u.NAh = u.N_A*u.h_c;                  % molar Planck constant
u.M_u = u.g/u.mol;                    % molar mass constant
u.K_J = 483597.8525e9*u.Hz/u.V;       % Josephson constant
u.K_J_90 = 483597.9*u.Hz/u.V;         % Josephson constant (conv. value)
u.F_c = 96485.33289*u.C/u.mol;        % Faraday constant
u.alpha = 7.2973525664e-3;            % fine-structure constant 
u.c1 = 3.741771790e-16*u.W/u.m^2;     % first radiation constant
u.c2 = 1.43877736e-2*u.m*u.K;         % second radiation constant
u.b_prime = 5.8789238e10*u.Hz/u.K;    % Wien frequency displ. law const.
u.b_c = 2.8977729e-3*u.m*u.K;         % Wien wavelength displ. law const.
u.R_air = 287.05287*u.J/u.kg/u.K;     % spec. gas const., air (ESDU 77022)
u.R_bar = 8.3144598*u.J/u.mol/u.K;    % molar gas constant

%---- currency ----
% For display purposes - not for exchange rates.
% See also mathworks.com/matlabcentral/fileexchange/47255

u.cent = u.currency/100;      % cent (currency)
u.Cent = u.cent;              % cent (currency)
u.pip = u.cent/100;           % pip (currency)
u.USD = u.currency;           % currency
u.EUR = u.currency;           % currency
u.GBP = u.currency;           % currency
u.JPY = u.currency;           % currency
u.AUD = u.currency;           % currency
u.CAD = u.currency;           % currency
u.CHF = u.currency;           % currency
u.CNY = u.currency;           % currency
u.dollar = u.currency;        % currency
u.franc = u.currency;         % currency

%---- used by symunits but not here ----
% gg -gauge
% land - league
% ha_US - US survey hectare
% molecule
% bps - bit per second
% Bd - baud
% B - byte
% bit
% HP_UK - British imperial horsepower
% PS_SAE - net horsepower (SAE J1349)
% PS_DIN - horsepower (DIN 70020)
% dry volumes


end

%   Original seed for this class by Rob deCarvalho.
%     http://www.mathworks.com/matlabcentral/fileexchange/authors/22148
%     http://www.mathworks.com/matlabcentral/fileexchange/10070
