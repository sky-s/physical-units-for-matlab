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
    'bit'      1
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
                'bit'       1
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
                'bit'       1
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

u.meter = u.m;
u.km = 1e3*u.m;               % kilometer
u.kilometer = u.km;
u.dm = 1e-1*u.m;              % decimeter
u.decimeter = u.dm;
u.cm = 1e-2*u.m;              % centimeter
u.centimeter = u.cm;
u.mm = 1e-3*u.m;              % millimeter
u.millimeter = u.mm;
u.um = 1e-6*u.m;              % micrometer
u.micrometer = u.um;
u.micron = u.um;              % micron
u.nm = 1e-9*u.m;              % nanometer
u.nanometer = u.nm;
u.pm = 1e-12*u.m;             % picometer
u.picometer = u.pm;
u.fm = 1e-15*u.m;             % femtometer
u.femtometer = u.fm;
u.fermi = u.fm;               % fermi
u.Ao = 1e-10*u.m;             % ångström
u.ang = u.Ao;                 % ångström
u.angstrom = u.ang;
u.angstroem = u.ang;
u.a0 = 0.52917721067e-10*u.m; % Bohr radius
u.a_0 = u.a0;                 % Bohr radius
u.BohrRadius = u.a0;
u.lP = 1.616229e-35*u.m;      % Planck length
u.PlanckLength = u.lP;
u.xu = 1.0021e-13*u.m;        % x unit
u.xUnit = u.xu;
u.xu_Cu = 1.00207697e-13*u.m; % x unit (copper)
u.xUnit_copper = u.xu_Cu;
u.xu_Mo = 1.00209952e-13*u.m; % x unit (molybdenum)
u.xUnit_molybdenum = u.xu_Mo;
u.in = 2.54*u.cm;             % inch
u.inch = u.in;
u.mil = 1e-3*u.in;            % mil
u.line = u.in/10;             % line
u.hand = 4*u.in;              % hand
u.span = 9*u.in;              % span
u.smoot = 67*u.in;            % smoot
u.ft = 12*u.in;               % foot
u.foot = u.ft;
u.ft_US = 1200/3937*u.m;      % US survey foot
u.foot_US = u.ft_US;          % US survey foot
u.kft = 1e3*u.ft;             % kilofoot
u.kilofoot = u.kft;
u.FL = 100*u.ft;              % flight level
u.flightLevel = u.FL;
u.yd = 3*u.ft;                % yard
u.yard = u.yd;
u.ftm = 6*u.ft;               % fathom
u.fathom = u.ftm;
u.li = 0.66*u.ft;             % link
u.link = u.li;
u.rod = 5.5*u.yd;             % rod
u.ch = 66*u.ft;               % chain
u.chain = u.ch;
u.fur = 220*u.yd;             % furlong
u.furlong = u.fur;
u.mi = 5280*u.ft;             % mile
u.mile = u.mi;
u.mi_US = 6336/3937*u.km;     % US survey mile
u.mile_US = u.mi_US;          % US survey mile
u.nmi = 1852*u.m;             % nautical mile
u.NM = u.nmi;                 % nautical mile
u.inm = u.nmi;                % nautical mile
u.nauticalMile = u.nmi;
u.nm_UK = 6080*u.ft;          % Imperial nautical mile
u.nmile = u.nm_UK;            % Imperial nautical mile
u.dataMile = 6000*u.ft;
u.au = 149597870.7*u.km;      % astronomical unit
u.astronomicalUnit = u.au;
u.pc = 648000/pi*u.au;        % parsec
u.parsec = u.pc;

%---- reciprocal length ----

u.dpt = 1/u.m;                % diopter
u.diopter = u.dpt;
u.R_inf = 1.0973731568508e7/u.m; % Rydberg constant
u.RydbergConstant = u.R_inf;

%---- area ----

u.square = 100*u.ft^2;        % square
u.ha = 10000*u.m^2;           % hectare
u.hectare = u.ha;
u.a = 100*u.m^2;              % are
u.are = u.a;
u.ac = 43560*u.ft^2;          % acre
u.acre = u.ac;
u.ro = 1/4*u.acre;            % rood
u.rood = u.ro;
u.twp = 36*u.mi^2;            % township
u.township = u.twp;
u.circ_mil = pi/4*u.mil^2;    % circular mil
u.circularMil = u.circ_mil;
u.circ_inch = pi/4*u.in^2;    % circular inch
u.circularInch = u.circ_inch;
u.b = 100*u.fm^2;             % barn
u.barn = u.b;

%---- volume ----

u.cc = u.cm^3;                % cubic centimeter
u.cubicCentimeter = u.cc;
u.L = 1000*u.cc;              % liter
u.l = u.L;                    % liter
u.liter = u.L;
u.dL = 100*u.cc;              % deciliter
u.dl = u.dL;                  % deciliter
u.deciliter = u.dl;
u.cL = 10*u.cc;               % centiliter
u.cl = u.cL;                  % centiliter
u.centiliter = u.cl;
u.mL = u.cc;                  % milliliter
u.ml = u.mL;                  % milliliter
u.milliliter = u.ml;
u.uL = u.mm^3;                % microliter
u.ul = u.uL;                  % microliter
u.microliter = u.ul;
u.kL = u.m^3;                 % kiloliter
u.kl = u.kL;                  % kiloliter
u.kiloliter = u.kl;
u.cuin = 16.387064*u.mL;      % cubic inch
u.cubicInch = u.cuin;
u.FBM = u.ft^2*u.in;          % board foot
u.boardFoot = u.FBM;
u.gal = 231*u.cuin;           % gallon (US)
u.gallon = u.gal;
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
u.fluidOunce = u.floz;        % US fluid ounce
u.floz_UK = u.gal_UK/160;     % British imperial fluid ounce
u.Tbls = u.floz/2;            % US tablespoon
u.tablespoon = u.Tbls;        % US tablespoon
u.tsp = u.Tbls/3;             % US teaspoon
u.teaspoon = u.tsp;           % US teaspoon
u.acft = u.acre*u.ft;         % acre-foot
u.acre_foot = u.acft;
u.acin = u.acre*u.in;         % acre-inch
u.acre_inch = u.acin;
u.bbl = 7056*u.in^3;          % US customary dry barrel
u.barrel = u.bbl;
u.fldr = u.floz/8;            % US customary fluid dram
u.fluidDram = u.fldr;
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
u.standardGravity = u.g0;
u.Gal = u.cm/u.s^2;           % gal

%---- force ----

u.N = u.kg*u.m/u.s^2;         % newton
u.newton = u.N;
u.kN = 1e3*u.N;               % kilonewton
u.kilonewton = u.kN;
u.MN = 1e6*u.N;               % meganewton
u.meganewton = u.MN;
u.mN = 1e-3*u.N;              % millinewton
u.millinewton = u.mN;
u.uN = 1e-6*u.N;              % micronewton
u.micronewton = u.uN;
u.dyn = 1e-5*u.N;             % dyne
u.dyne = u.dyn;
u.lbf = 4.4482216152605*u.N;  % pound force
u.poundForce = u.lbf;
u.kip = 1000*u.lbf;           % kip
u.kilopoundForce = u.kip;
u.kgf = u.kg*u.g0;            % kilogram force
u.kilogramForce = u.kgf;
u.kp = u.kgf;                 % kilopond
u.kilopond = u.kp;
u.p = u.kp/1000;              % pond
u.pond = u.p;
u.sn = u.kN;                  % sthène
u.sthene = u.sn;

%---- mass ----

u.kilogram = u.kg;
u.kilo = u.kg;                % kilogram
u.g = 1e-3*u.kg;              % gram
u.gram = u.g;
u.mg = 1e-3*u.gram;           % milligram
u.milligram = u.mg;
u.ug = 1e-6*u.gram;           % microgram
u.microgram = u.ug;
u.Mg = 1e6*u.gram;            % Megagram/metric tonne
u.Megagram = u.Mg;
u.t = 1000*u.kg;              % metric tonne
u.tonne = u.t;                % metric ton
u.Mt = 1e6*u.t;               % metric megatonne
u.megatonne = u.Mt;
u.lbm = 0.45359237*u.kg;      % pound mass
u.poundMass = u.lbm;
u.lb = u.lbm;                 % pound mass
u.pound = u.lb;
u.tn = 2000*u.lbm;            % US customary short ton
u.ton = u.tn;                 % US customary short ton
u.ton_UK = 2240*u.lbm;        % British imperial ton
u.st = 14*u.lbm;              % stone
u.stone = u.st;
u.cwt = 100*u.lbm;            % US customary short hundredweight
u.hundredweight = u.cwt;
u.cwt_UK = 8*u.stone;         % British imperial short hundredweight
u.quarter = u.cwt_UK/4;       % British imperial quarter
u.slug = u.lbf/(u.ft/u.s^2);  % slug
u.oz = u.lbm/16;              % ounce
u.ounce = u.oz;
u.dr = u.oz/16;               % dram
u.dram = u.dr;
u.gr = u.lbm/7000;            % grain
u.grain = u.gr;
u.ct = 200*u.mg;              % carat
u.carat = u.ct;
u.amu = 1.660539040e-27*u.kg; % atomic mass unit
u.atomicMassUnit = u.amu;
u.Da = u.amu;                 % atomic mass unit
u.dalton = u.Da;
u.mu = u.amu;                 % atomic mass unit
u.mP = 2.176470e-8*u.kg;      % Planck mass
u.PlanckMass = u.mP;
u.m_e = 9.10938356e-31*u.kg;  % electron mass
u.electronMass = u.m_e;
u.mug = u.kgf/(u.m/u.s^2);    % metric slug
u.metricSlug = u.mug;
u.hyl = u.mug;                % hyl
u.TMU = u.mug;                % technische Masseneinheit
u.technischeMasseneinheit = u.TMU;

%---- more force ----

u.pdl = u.lbm*u.ft/u.s^2;     % poundal
u.poundal = u.pdl;
u.gf = u.gram*u.g0;           % gram force
u.gramForce = u.gf;
u.ozf = u.oz*u.g0;            % ounce force
u.ounceForce = u.ozf;
u.tonf = u.tn*u.g0;           % short ton force
u.tonForce = u.tonf;

%---- mass per length ----

u.den = u.gram/(9*u.km);      % denier
u.denier = u.den;
u.tex = u.gram/u.km;          % tex
u.dtex = u.tex/10;            % decitex
u.decitex = u.dtex;

%---- time ----

u.second = u.s;
u.ms = 1e-3*u.s;              % millisecond
u.millisecond = u.ms;
u.us = 1e-6*u.s;              % microsecond
u.microsecond = u.us;
u.ns = 1e-9*u.s;              % nanosecond
u.nanosecond = u.ns;
u.ps = 1e-12*u.s;             % picosecond
u.picosecond = u.ps;
u.fs = 1e-15*u.s;             % femtosecond
u.femtosecond = u.fs;
u.tP = 5.39116e-44*u.s;       % Planck time
u.PlanckTime = u.tP;
u.min = 60*u.s;               % minute
u.minute = u.min;
u.h = 60*u.min;               % hour
u.hr = u.h;                   % hour
u.hour = u.hr;
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
u.hertz = u.Hz;
u.kHz = 1e3*u.Hz;             % kilohertz
u.kilohertz = u.kHz;
u.MHz = 1e6*u.Hz;             % megahertz
u.megahertz = u.MHz;
u.GHz = 1e9*u.Hz;             % gigahertz
u.gigahertz = u.GHz;
u.THz = 1e12*u.Hz;            % terahertz
u.terahertz = u.THz;
u.Bd = 1/u.s;                 % baud
u.baud = u.Bd;

%---- energy ----

u.Nm = u.N*u.m;               % newton-meter
u.newton_meter = u.Nm;
u.J = u.Nm;                   % joule
u.joule = u.J;
u.kJ = 1e3*u.J;               % kilojoule
u.kilojoule = u.kJ;
u.MJ = 1e6*u.J;               % megajoule
u.megajoule = u.MJ;
u.GJ = 1e9*u.J;               % gigajoule
u.gigajoule = u.GJ;
u.mJ = 1e-3*u.J;              % millijoule
u.millijoule = u.mJ;
u.uJ = 1e-6*u.J;              % microjoule
u.microjoule = u.uJ;
u.nJ = 1e-9*u.J;              % nanojoule
u.nanojoule = u.nJ;
u.eV = 1.6021766208e-19*u.J;  % electronvolt
u.electronvolt = u.eV;
u.BTU = 1055.06*u.J;          % British thermal unit (ISO)
u.Btu = u.BTU;                % British thermal unit (ISO)
u.britishThermalUnit = u.Btu;
u.Btu_IT = 1055.0559*u.J;     % British thermal unit (International Table)
u.Btu_th = 1054.3503*u.J;     % British thermal unit (thermochemical)
u.kpm = u.kp*u.m;             % kilopond-meter
u.kilopond_meter = u.kpm;
u.Ws = u.J;                   % watt-second
u.watt_second = u.Ws;
u.kWh = 3.6e6*u.J;            % kilowatt-hour
u.kilowatt_hour = u.kWh;
u.Wh = 3.6e3*u.J;             % watt-hour
u.watt_hour = u.Wh;
u.cal = 4.1868*u.J;           % calorie (International Table)
u.calorie = u.cal;
u.cal_IT = u.cal;             % calorie (International Table)
u.cal_4 = 4.204*u.J;          % calorie (4°C)
u.cal_15 = 4.1855*u.J;        % calorie (15°C)
u.cal_20 = 4.182*u.J;         % calorie (20°C)
u.cal_mean = 4.190*u.J; 	  % calorie (mean)
u.cal_th = 4.184*u.J;         % calorie (thermochemical)
u.kcal = 1e3*u.cal;           % kilocalorie
u.kilocalorie = u.kcal;
u.kcal_IT = 1e3*u.cal_IT;     % kilocalorie (International Table)
u.Cal = u.kcal;               % large calorie / food calorie
u.foodCalorie = u.Cal;
u.largeCalorie = u.Cal;
u.kcal_4 = 1e3*u.cal_4;       % kilocalorie (4°C)
u.kcal_15 = 1e3*u.cal_15;     % kilocalorie (15°C)
u.kcal_20 = 1e3*u.cal_20;     % kilocalorie (20°C)
u.kcal_mean = 1e3*u.cal_mean; % kilocalorie (mean)
u.kcal_th = 1e3*u.cal_th;     % kilocalorie (thermochemical)
u.erg = 1e-7*u.J;             % en.wikipedia.org/wiki/Erg
u.E_h = 4.359744650e-18*u.J;  % Hartree energy
u.Ha = u.E_h;                 % hartree
u.hartree = u.Ha;
u.thm = 1e5*u.BTU;            % therm
u.therm = u.thm;
u.quad = 1e15*u.BTU;          % quad

%---- temperature ----
% For reference: °C = °K-273.15; °F = °R-459.67.

u.kelvin = u.K;
u.R = u.K*5/9;                % rankine (°F = °R-459.67)
u.rankine = u.R;
u.mK = 1e-3*u.K;              % millikelvin
u.millikelvin = u.mK;
u.uK = 1e-6*u.K;              % microkelvin
u.microkelvin = u.uK;
u.nK = 1e-9*u.K;              % nanokelvin
u.nanokelvin = u.nK;
u.deltaK = u.K;               % kelvin (relative temperature)
u.deltadegC = u.K;            % celsius (relative, °C = °K-273.15)
u.deltadegR = u.R;            % rankine (relative temperature)
u.deltadegF = u.R;            % fahrenheit (relative, °F = °R-459.67)
u.TP = 1.416808e32*u.K;       % Planck temperature
u.PlanckTemperature = u.TP;

%---- pressure ----

u.Pa = u.N/u.m^2;             % pascal
u.pascal = u.Pa;
u.mPa = 1e-3*u.Pa;            % millipascal
u.millipascal = u.mPa;
u.uPa = 1e-6*u.Pa;            % micropascal
u.micropascal = u.uPa;
u.kPa = 1e3*u.Pa;             % kilopascal
u.kilopascal = u.kPa;
u.MPa = 1e6*u.Pa;             % megapascal
u.megapascal = u.MPa;
u.GPa = 1e9*u.Pa;             % gigapascal
u.gigapascal = u.GPa;
u.torr = 133.322*u.Pa;        % torr
u.Torr = u.torr;              % torr
u.mtorr = 1e-3*u.torr;        % millitorr
u.millitorr = u.mtorr;
u.bar = 1e5*u.Pa;             % bar
u.mbar = 1e-3*u.bar;          % millibar
u.millibar = u.mbar;
u.kbar = 1e3*u.bar;           % kilobar
u.kilobar = u.kbar;
u.atm = 101325*u.Pa;          % standard atmosphere
u.atmosphere = u.atm;
u.standardAtmosphere = u.atm;
u.at = u.kgf/u.cm^2;          % technical atmosphere
u.technicalAtmosphere = u.at;
u.psi = u.lbf/u.in^2;         % pound force per square inch
u.poundPerSquareInch = u.psi;
u.ksi = 1e3*u.psi;            % kip per square inch
u.kipPerSquareInch = u.ksi;
u.Msi = 1e6*u.psi;            % million pound force per square inch
u.megapoundPerSquareInch = u.Msi;
u.psf = u.lbf/u.ft^2;         % pound force per square foot
u.poundPerSquareFoot = u.psf;
u.ksf = u.kip/u.ft^2;         % kip per square foot
u.kipPerSquareFoot = u.ksf;
u.Ba = 0.1*u.Pa;              % barye
u.barye = u.Ba;
u.pz = u.kPa;                 % pièze
u.pieze = u.pz;
u.mmHg = 13.5951*u.kgf/u.m^2; % millimeter of mercury
u.millimeterMercury = u.mmHg;
u.cmHg = 10*u.mmHg;           % centimeter of mercury
u.centimeterMercury = u.cmHg;
u.mHg = 1e3*u.mmHg;           % meter of mercury
u.meterMercury = u.mHg;
u.inHg = 2.54*u.cmHg;         % inch of mercury
u.inchMercury = u.inHg;
u.ftHg = 12*u.inHg;           % foot of mercury
u.footMercury = u.ftHg;
u.mmH20 = u.kgf/u.m^2;        % millimeter of water (density 1 g/cc)
u.mmAq = u.mmH20;             % millimeter of water
u.millimeterWater = u.mmH20;
u.cmH20 = 10*u.mmH20;         % centimeter of water
u.cmAq = u.cmH20;             % centimeter of water
u.centimeterWater = u.cmH20;
u.mH20 = 1e3*u.mmH20;         % meter of water
u.mAq = u.mH20;               % meter of water
u.meterWater = u.mH20;
u.inH20 = 2.54*u.cmH20;       % inch of water
u.inAq = u.inH20;             % inch of water
u.inchWater = u.inH20;
u.wc = u.inH20;               % inch water column
u.inchWaterColumn = u.wc;
u.ftH20 = 12*u.inH20;         % foot of water
u.ftAq = u.ftH20;             % foot of water
u.footWater = u.ftH20;

%---- viscosity ----

u.St = u.cm^2/u.s;            % stokes
u.stokes = u.St;
u.cSt = u.St/100;             % centistokes
u.centistokes = u.cSt;
u.newt = u.in^2/u.s;          % newt
u.P = u.Pa*u.s / 10;          % poise
u.poise = u.P;
u.cP = u.mPa*u.s;             % centipoise
u.centipoise = u.cP;
u.reyn = u.lbf*u.s/u.in^2;    % reyn

%---- power ----

u.W = u.J/u.s;                % watt
u.watt = u.W;
u.kW = 1e3*u.W;               % kilowatt
u.kilowatt = u.kW;
u.MW = 1e6*u.W;               % megawatt
u.megawatt = u.MW;
u.GW = 1e9*u.W;               % gigawatt
u.gigawatt = u.GW;
u.TW = 1e12*u.W;              % terawatt
u.terawatt = u.TW;
u.mW = 1e-3*u.W;              % milliwatt
u.milliwatt = u.mW;
u.uW = 1e-6*u.W;              % microwatt
u.microwatt = u.uW;
u.nW = 1e-9*u.W;              % nanowatt
u.nanowatt = u.nW;
u.pW = 1e-12*u.W;             % picowatt
u.picowatt = u.pW;
u.hp = 550*u.ft*u.lbf/u.s;    % mechanical horsepower (550 ft-lbf/s)
u.horsepower = u.hp;
u.HP_I = u.hp;                % mechanical horsepower (550 ft-lbf/s)
u.hpE = 746*u.W;              % electrical horsepower
u.HP_E = u.hpE;               % electrical horsepower
u.electricalHorsepower = u.hp;
u.PS = 75*u.kg*u.g0*u.m/u.s;  % metric horsepower (DIN 66036)
u.HP = u.PS;                  % metric horsepower (DIN 66036)
u.HP_DIN = u.PS;              % metric horsepower (DIN 66036)
u.metricHorsepower = u.PS;

%---- current ----

u.amp = u.A;                  % ampere
u.ampere = u.A;
u.mA = 1e-3*u.A;              % milliampere
u.milliampere = u.mA;
u.uA = 1e-6*u.A;              % microampere
u.microampere = u.uA;
u.nA = 1e-9*u.A;              % nanoampere
u.nanoampere = u.nA;
u.pA = 1e-12*u.A;             % picoampere
u.picoampere = u.pA;
u.kA = 1e3*u.A;               % kiloampere
u.kiloampere = u.kA;
u.abA = 10*u.A;               % abampere
u.abampere = u.abA;
u.Bi = u.abA;                 % biot
u.biot = u.Bi;

%---- charge ----

u.C = u.A*u.s;                % coulomb
u.coulomb = u.C;
u.mC = 1e-3*u.C;              % millicoulomb
u.millicoulomb = u.mC;
u.uC = 1e-6*u.C;              % microcoulomb
u.microcoulomb = u.uC;
u.nC = 1e-9*u.C;              % nanocoulomb
u.nanocoulomb = u.nC;
u.pC = 1e-12*u.C;             % picocoulomb
u.picocoulomb = u.pC;
u.abC = 10*u.C;               % abcoulomb
u.aC = u.abC;                 % abcoulomb
u.abcoulomb = u.abC;
u.statC = u.dyn^(1/2)*u.cm;   % statcoulomb
u.statcoulomb = u.statC;
u.Fr = u.statC;               % franklin
u.franklin = u.Fr;
u.esu = u.statC;              % electrostatic unit of charge
u.electrostaticUnitCharge = u.esu;
u.e = 1.6021766208e-19*u.C;   % elementary charge
u.elementaryCharge = u.e;
u.mAh = u.mA*u.hr;            % milliamp-hour
u.milliamp_hour = u.mAh;
u.Ah = u.A*u.hr;              % amp-hour
u.amp_hour = u.Ah;

%---- voltage ----

u.V = 1*u.J/u.C;              % volt
u.volt = u.V;
u.kV = 1e3*u.V;               % kilovolt
u.kilovolt = u.kV;
u.MV = 1e6*u.V;               % megavolt
u.megavolt = u.MV;
u.GV = 1e9*u.V;               % gigavolt
u.gigavolt = u.GV;
u.mV = 1e-3*u.V;              % millivolt
u.millivolt = u.mV;
u.uV = 1e-6*u.V;              % microvolt
u.microvolt = u.uV;

%---- resistance/conductance ----

u.Ohm = u.V/u.A;              % ohm
u.GOhm = 1e9*u.Ohm;           % gigaohm
u.gigaohm = u.GOhm;
u.MOhm = 1e6*u.Ohm;           % megaohm
u.megaohm = u.MOhm;
u.kOhm = 1e3*u.Ohm;           % kiloohm
u.kiloohm = u.kOhm;
u.mOhm = 1e-3*u.Ohm;          % milliohm
u.milliohm = u.mOhm;
u.uOhm = 1e-6*u.Ohm;          % microohm
u.microohm = u.uOhm;
u.nOhm = 1e-9*u.Ohm;          % nanoohm
u.nanoohm = u.nOhm;
u.abOhm = u.nOhm;             % abohm
u.Z0 = 376.730313461*u.Ohm;   % characteristic impedance of vacuum
u.impedanceOfVacuum = u.Z0;
u.R_K = 25812.8074555*u.Ohm;  % von Klitzing constant
u.vonKlitzingConstant = u.R_K;
u.R_K_90 = 25812.807*u.Ohm;   % von Klitzing constant (conventional value)
u.vonKlitzingConstant_conv = u.R_K_90;
u.S = 1/u.Ohm;                % siemens
u.siemens = u.S;
u.mS = 1e-3*u.S;              % millisiemens
u.millisiemens = u.mS;
u.uS = 1e-6*u.S;              % microsiemens
u.microsiemens = u.uS;
u.nS = 1e-9*u.S;              % nanosiemens
u.nanosiemens = u.nS;
u.G0 = 7.7480917310e-5*u.S;   % conductance quantum
u.conductanceQuantum = u.G0;

%---- capacitance ----

u.F = u.A*u.s/u.V;            % farad
u.farad = u.F;
u.mF = 1e-3*u.F;              % millifarad
u.millifarad = u.mF;
u.uF = 1e-6*u.F;              % microfarad
u.microfarad = u.uF;
u.nF = 1e-9*u.F;              % nanofarad
u.nanofarad = u.nF;
u.pF = 1e-12*u.F;             % picofarad
u.picofarad = u.pF;

%---- inductance ----

u.H = u.Ohm*u.s;              % henry
u.henry = u.H;
u.mH = 1e-3*u.H;              % millihenry
u.millihenry = u.mH;
u.uH = 1e-6*u.H;              % microhenry
u.microhenry = u.uH;
u.nH = 1e-9*u.H;              % nanohenry
u.nanohenry = u.nH;
u.abH = u.nH;                 % abhenry
u.abhenry = u.abH;
u.kH = 1e3*u.H;               % kilohenry
u.kilohenry = u.kH;
u.MH = 1e6*u.H;               % megahenry
u.megahenry = u.MH;
u.GH = 1e9*u.H;               % gigahenry
u.gigahenry = u.GH;

%---- EM ----

u.T = 1*u.N/(u.A*u.m);        % tesla
u.tesla = u.T;
u.Gs = 1e-4*u.T;              % gauss
u.gauss = u.Gs;
u.Wb = u.V*u.s;               % weber
u.weber = u.Wb;
u.Mx = 1e-8*u.Wb;             % maxwell
u.maxwell = u.Mx;
u.mWb = u.Wb/1000;            % milliweber
u.milliweber = u.mWb;
u.uWb = 1e-6*u.Wb;            % microweber
u.microweber = u.uWb;
u.nWb = 1e-9*u.Wb;            % nanoweber
u.nanoweber = u.nWb;
u.Oe = 250/pi*u.A/u.m;        % oersted
u.oersted = u.Oe;
u.Gb = 2.5/pi*u.A;            % gilbert
u.gilbert = u.Gb;

%---- non-dimensionals ----

u.percent = 0.01;             % %
u.pct = u.percent;            % %
u.permil = 0.001;             % ‰
u.permill = u.permil;         % ‰
u.permille = u.permil;        % ‰
u.permyriad = 1e-4;           % permyriad
u.bp = u.permyriad;           % basis point
u.basisPoint = u.bp;
u.ppm = 1e-6;                 % part per million
u.partPerMillion = u.ppm;
u.ppb = 1e-9;                 % part per billion
u.partPerBillion = u.ppb;
u.ppt = 1e-12;                % part per trillion
u.partPerTrillion = u.ppt;
u.ppq = 1e-15;                % part per quadrillion
u.partPerQuadrillion = u.ppq;

%---- angles ----
% Note: angles are dimensionless

u.rad = 1;                    % radian
u.radian = u.rad;
u.sr = 1;                     % steradian
u.steradian = u.sr;
u.turn = 2*pi*u.rad;          % turn
u.rev = u.turn;               % revolution = 2*pi radians
u.revolution = u.rev;
u.deg = u.turn/360;           % degree
u.degree = u.deg;
u.arcmin = u.deg/60;          % arcminute
u.arcminute = u.arcmin;
u.arcsec = u.arcmin/60;       % arcsecond
u.arcsecond = u.arcsec;
u.grad = u.turn/400;          % gradian
u.gradian = u.grad;

%---- rotational speed ----

u.rpm = u.rev/u.min;          % revolution per minute
u.revolutionPerMinute = u.rpm;
u.rps = u.rev/u.s;            % revolution per second
u.revolutionPerSecond = u.rps;

%---- velocity ----

u.mps = u.m/u.s;              % meter per second
u.meterPerSecond = u.mps;
u.kyne = u.cm/u.s;            % kyne
u.Kyne = u.kyne;              % kyne
u.fps = u.ft/u.s;             % foot per second
u.footPerSecond = u.fps;
u.fpm = u.ft/u.min;           % foot per minute
u.footPerMinute = u.fpm;
u.kt = u.nmi/u.hr;            % knot
u.kn = u.kt;                  % knot
u.kts = u.kt;                 % knot
u.knot = u.kt;
u.knot_UK = u.nm_UK/u.hr;     % British imperial knot
u.KTAS = u.kt;                % knot
u.nmph = u.kt;                % nautical mile per hour
u.nauticalMilePerHour = u.nmph;
u.kph = u.km/u.hr;            % kilometer per hour
u.kmh = u.kph;                % kilometer per hour
u.kilometerPerHour = u.kmh;
u.mph = u.mi/u.hr;            % mile per hour
u.milePerHour = u.mph;

%---- volume flow rate ----

u.cfm = u.ft^3/u.min;         % cubic foot per minute
u.cubicFootPerMinute = u.cfm;
u.cfs = u.ft^3/u.s;           % cubic foot per second
u.cubicFootPerSecond = u.cfs;
u.gpm = u.gal/u.min;          % US customary gallon per minute
u.gallonPerMinute = u.gpm;
u.gpm_UK = u.gal_UK/u.min;    % British imperial gallon per minute
u.lpm = u.l/u.min;            % liter per minute
u.literPerMinute = u.lpm;

%---- fuel economy ----

u.l_100km = u.l/(100*u.km);   % liter per 100 km
u.literPer100kilometer = u.l_100km;
u.mpg = u.mi/u.gal;           % mile per gallon
u.milePerGallon = u.mpg;

%---- Luminance etc. ----

u.candela = u.cd;
u.asb = u.cd/u.m^2;           % apostilb
u.apostilb = u.asb;
u.sb = u.cd/u.cm^2;           % stilb
u.stilb = u.sb;
u.ph = 1e4*u.cd*u.sr/u.m^2;   % phot
u.phot = u.ph;
u.cp = 0.981*u.cd;            % candlepower
u.candlepower = u.cp;
u.lm = u.cd*u.sr;             % lumen
u.lumen = u.lm;
u.lx = u.lm/u.m^2;            % lux
u.lux = u.lx;
u.nx = 1e-3*u.lx;             % nox
u.nox = u.nx;

%---- other derived SI ----

u.mole = u.mol;
u.kat = u.mol/u.s;            % katal
u.katal = u.kat;
u.M = u.mol/u.L;              % molar
u.molar = u.M;
u.molarity = u.M;             % molarity
u.Nms = u.N*u.m*u.s;          % newton-meter-second
u.newton_meter_second = u.Nms;

%---- radiation ----

u.Gy = u.J/u.kg;              % gray
u.gray = u.Gy;
u.Sv = u.J/u.kg;              % sievert
u.sievert = u.Sv;
u.Rad =  u.Gy/100;            % absorbed radiation dose
u.rem = u.Sv/100;             % roentgen equivalent man
u.roentgenEquivalentMan = u.rem;
u.roentgen = 2.58e-4*u.C/u.kg;% roentgen
u.Ly = u.cal_th/u.cm^2;       % langley
u.lan = u.Ly;                 % langley
u.langley = u.lan;
u.Bq = 1/u.s;                 % becquerel
u.becquerel = u.Bq;
u.Ci = 3.7e10*u.Bq;           % curie
u.curie = u.Ci;

%---- constants ----

u.k_B = 1.38064852e-23*u.J/u.K;       % Boltzmann constant
u.BoltzmannConstant = u.k_B;
u.sigma_SB = 5.670367e-8*u.W/(u.m^2*u.K^4); % Stefan–Boltzmann constant
u.Stefan_BoltzmannConstant = u.sigma_SB;
u.h_c = 6.626070040e-34*u.J*u.s;      % Planck constant
u.PlanckConstant = u.h_c;
u.h_bar = u.h/(2*pi);                 % Dirac constant
u.DiracConstant = u.h_bar;
u.mu_B = 9.274009994e-24*u.J/u.T;     % Bohr magneton
u.BohrMagneton = u.mu_B;
u.mu_N = 5.050783699e-27*u.J/u.T;     % nuclear magneton
u.nuclearMagneton = u.mu_N;
u.c = 299792458*u.m/u.s;              % speed of light in vacuum
u.c_0 = u.c;                          % speed of light in vacuum
u.lightSpeed = u.c;
u.speedOfLight = u.c;
u.ly = u.c*u.year;                    % light-year
u.lightYear = u.ly;                   % light-year
u.mu0 = pi*4e-7*u.N/u.A^2;            % vacuum permeability
u.vacuumPermeability = u.mu0;
u.eps0 = u.c^-2/u.mu0;                % vacuum permittivity
u.vacuumPermittivity = u.eps0;
u.G = 6.67408e-11*u.m^3/u.kg/u.s^2;   % gravitational constant
u.gravitationalConstant = u.G;
u.N_A = 6.022140857e23/u.mol;         % Avogadro constant
u.NA = u.N_A;                         % Avogadro constant
u.AvogadroConstant = u.N_A;
u.NAh = u.N_A*u.h_c;                  % molar Planck constant
u.molarPlanckConstant = u.NAh;
u.M_u = u.g/u.mol;                    % molar mass constant
u.molarMassConstant = u.M_u;
u.K_J = 483597.8525e9*u.Hz/u.V;       % Josephson constant
u.JosephsonConstant = u.K_J;
u.K_J_90 = 483597.9*u.Hz/u.V;         % Josephson constant (conv. value)
u.JosephsonConstant_conv = u.K_J_90;
u.F_c = 96485.33289*u.C/u.mol;        % Faraday constant
u.FaradayConstant = u.F_c;
u.alpha = 7.2973525664e-3;            % fine-structure constant
u.fine_structureConstant = u.alpha;
u.SommerfeldConstant = u.alpha;
u.c1 = 3.741771790e-16*u.W/u.m^2;     % first radiation constant
u.firstRadiationConstant = u.c1;
u.c2 = 1.43877736e-2*u.m*u.K;         % second radiation constant
u.secondRadiationConstant = u.c2;
u.b_prime = 5.8789238e10*u.Hz/u.K;    % Wien frequency displ. law const.
u.WienFrequencyDisplacementLawConstant = u.b_prime;
u.b_c = 2.8977729e-3*u.m*u.K;         % Wien wavelength displ. law const.
u.WienWavelengthDisplacementLawConstant = u.b_c;
u.R_air = 287.05287*u.J/u.kg/u.K;     % spec. gas const., air (ESDU 77022)
u.specificGasConstant_air = u.R_air;
u.R_bar = 8.3144598*u.J/u.mol/u.K;    % molar gas constant
u.molarGasConstant = u.R_bar;
u.radarStatuteMile = 2*u.mi/u.c;
u.radarNauticalMile = 2*u.NM/u.c;
u.radarDataMile = 2*u.dataMile/u.c;
u.radarKilometer = 2*u.km/u.c;

%---- digital information ----

u.nibble = 4*u.bit;                   % nibble
u.B = 8*u.bit;                        % byte
u.byte = u.B;                         % byte
u.octet = u.B;                        % octet
u.kB = 1e3*u.B;                       % kilobyte
u.kilobyte = u.kB;
u.MB = 1e6*u.B;                       % megabyte
u.megabyte = u.MB;
u.GB = 1e9*u.B;                       % gigabyte
u.gigabyte = u.GB;
u.TB = 1e12*u.B;                      % terabyte
u.terabyte = u.TB;
u.PB = 1e15*u.B;                      % petabyte
u.petabyte = u.PB;
u.EB = 1e18*u.B;                      % exabyte
u.exabyte = u.EB;
u.Kibit = 2^10*u.bit;                 % kibibit
u.kibibit = u.Kibit;
u.KiB = 2^10*u.B;                     % kibibyte
u.KB = u.KiB;                         % kibibyte
u.kibibyte = u.KB;
u.Mibit = 2^20*u.bit;                 % mebibit
u.mebibit = u.Mibit;
u.MiB = 2^20*u.B;                     % mebibyte
u.mebibyte = u.MiB;
u.Gibit = 2^30*u.bit;                 % gibibit
u.gibibit = u.Gibit;
u.GiB = 2^30*u.B;                     % gibibyte
u.gibibyte = u.GiB;
u.Tibit = 2^40*u.bit;                 % tebibit
u.tebibit = u.Tibit;
u.TiB = 2^40*u.B;                     % tebibyte
u.tebibyte = u.TiB;
u.Pibit = 2^50*u.bit;                 % pebibit
u.pebibit = u.Pibit;
u.PiB = 2^50*u.B;                     % pebibyte
u.pebibyte = u.PiB;
u.Eibit = 2^60*u.bit;                 % exbibit
u.exbibit = u.Eibit;
u.EiB = 2^60*u.B;                     % exbibyte
u.exbibyte = u.EiB;
u.bps = u.bit/u.s;                    % bit per second
u.bitPerSecond = u.bps;
u.kbps = 1e3*u.bps;                   % kilobit per second
u.kilobitPerSecond = u.kbps;
u.Mbps = 1e6*u.bps;                   % megabit per second
u.megabitPerSecond = u.Mbps;
u.Gbps = 1e9*u.bps;                   % gigabit per second
u.gigabitPerSecond = u.Gbps;
u.Tbps = 1e12*u.bps;                  % terabit per second
u.terabitPerSecond = u.Tbps;

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

%---- used by symunit but not here ----
% gg -gauge
% land - league
% ha_US - US survey hectare
% molecule
% HP_UK - British imperial horsepower
% PS_SAE - net horsepower (SAE J1349)
% PS_DIN - horsepower (DIN 70020)
% dry volumes


end

%   Original seed for this class by Rob deCarvalho.
%     http://www.mathworks.com/matlabcentral/fileexchange/authors/22148
%     http://www.mathworks.com/matlabcentral/fileexchange/10070
