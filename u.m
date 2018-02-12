classdef u < handle
% u  Physical units.
%
%   If the Physical Units Toolbox is on your MATLAB path, there is nothing to
%   initialize, add to your workspace, or pass to functions - simply
%   multiply/divide by u.(unitName) to attach physical units to a variable. For
%   example, to define a speed using a supported unit: carSpeed = 100 * u.kph.
%   Or, define a speed with an unsupported unit as a combination of supported
%   units: snailSpeed = 20 * u.m/u.week.
%
%   Calling u by itself will display all available units in u.
%
%   Variables with physical units attached are of the class DimVar
%   ("dimenensioned variable"). Math operations performed on dimensioned
%   variables will automatically perform dimensional analysis and can create new
%   units or cancel units and return a normal variable.
%
%   Variables with units listed in displayUnits display, plot, etc. in terms of
%   those units. A variable with units not in the list will display as a
%   combination of fundamental base units (mass, length, time, temperature,
%   ...). To customize the list, see displayUnits. For more advanced
%   customization of the base units themselves, see baseUnitSystem.
%
%   Customization is set by calls to displayUnits and/or baseUnitSystem (either
%   function files or variables in the base workspace). Tailor preferences for a
%   specific project by defining these variables at the top of a script (before
%   any units are called) or placing unique versions of the files in a project's
%   directory. Be sure to clear the class when changing projects or else the old
%   customizations will remain in effect.
%
%   Some MATLAB functions won't accept variables with physical units. Most of
%   the time displayingvalue, which returns value in terms of preferred display
%   units, will be the appropriate tool, but there is also double and u2num.
%
%   Example 1: Shaft power.
%       rotationSpeed = 2500 * u.rpm;
%       torque = 95 * str2u('ft-lbf');  % Use alternate string-based definition.
%       power = rotationSpeed * torque; % Returns variable with units of power.
%       horsePower = power / u.hp;      % Convert/cancel units.
%
%   Example 2: Unit conversion.
%       100 * u.acre/u.ha;  % Convert 100 acres to hectares.
%       u.st/u.kg;          % Return conversion factor for stone to kilos.
%
%   See also displayUnits, clear, str2u, symunit, displayingvalue,
%     dispdisp - http://www.mathworks.com/matlabcentral/fileexchange/48637.

%   Copyright Sky Sartorius
%   www.mathworks.com/matlabcentral/fileexchange/authors/101715
%   github.com/sky-s/physical-units-for-matlab

properties (Hidden, Constant = true)
    %% User-defined base and display:
    % Establishes base unit system and preferences based on baseUnitSystem and
    % displayUnits.
    
    baseUnitSystem =    evalin('base','baseUnitSystem')
    dispUnits =         evalin('base','displayUnits')

    coreUnits = buildCoreUnits(u.baseUnitSystem);
end
properties (Constant = true)

    %% Core units:
    baseNames = u.baseUnitSystem(:,1)'

    m           = u.coreUnits.m         % meter
    kg          = u.coreUnits.kg        % kilogram
    s           = u.coreUnits.s         % second
    A           = u.coreUnits.A         % ampere
    K           = u.coreUnits.K         % kelvin (°C = °K-273.15)
    mol         = u.coreUnits.mol       % mole
    cd          = u.coreUnits.cd        % candela
    bit         = u.coreUnits.bit       % bit
    currency    = u.coreUnits.currency  % currency

    %% Derived units list:
    % References:
    % http://physics.nist.gov/cuu/Constants/index.html
    % http://www.translatorscafe.com/unit-converter
    % http://en.wikipedia.org
    % http://www.efunda.com/units/index.cfm
    
    %---- length ----

    meter = u.m;
    km = 1e3*u.m;               % kilometer
    kilometer = u.km;
    dm = 1e-1*u.m;              % decimeter
    decimeter = u.dm;
    cm = 1e-2*u.m;              % centimeter
    centimeter = u.cm;
    mm = 1e-3*u.m;              % millimeter
    millimeter = u.mm;
    um = 1e-6*u.m;              % micrometer
    micrometer = u.um;
    micron = u.um;              % micron
    nm = 1e-9*u.m;              % nanometer
    nanometer = u.nm;
    pm = 1e-12*u.m;             % picometer
    picometer = u.pm;
    fm = 1e-15*u.m;             % femtometer
    femtometer = u.fm;
    fermi = u.fm;               % fermi
    Ao = 1e-10*u.m;             % ångström
    ang = u.Ao;                 % ångström
    angstrom = u.ang;
    angstroem = u.ang;
    a0 = 0.52917721067e-10*u.m; % Bohr radius
    a_0 = u.a0;                 % Bohr radius
    BohrRadius = u.a0;
    lP = 1.616229e-35*u.m;      % Planck length
    PlanckLength = u.lP;
    xu = 1.0021e-13*u.m;        % x unit
    xUnit = u.xu;
    xu_Cu = 1.00207697e-13*u.m; % x unit (copper)
    xUnit_copper = u.xu_Cu;
    xu_Mo = 1.00209952e-13*u.m; % x unit (molybdenum)
    xUnit_molybdenum = u.xu_Mo;
    in = 2.54*u.cm;             % inch
    inch = u.in;
    mil = 1e-3*u.in;            % mil
    line = u.in/10;             % line
    hand = 4*u.in;              % hand
    span = 9*u.in;              % span
    smoot = 67*u.in;            % smoot
    ft = 12*u.in;               % foot
    foot = u.ft;
    ft_US = 1200/3937*u.m;      % US survey foot
    foot_US = u.ft_US;          % US survey foot
    kft = 1e3*u.ft;             % kilofoot
    kilofoot = u.kft;
    FL = 100*u.ft;              % flight level
    flightLevel = u.FL;
    yd = 3*u.ft;                % yard
    yard = u.yd;
    ftm = 6*u.ft;               % fathom
    fathom = u.ftm;
    li = 0.66*u.ft;             % link
    link = u.li;
    rod = 5.5*u.yd;             % rod
    ch = 66*u.ft;               % chain
    chain = u.ch;
    fur = 220*u.yd;             % furlong
    furlong = u.fur;
    mi = 5280*u.ft;             % mile
    mile = u.mi;
    mi_US = 6336/3937*u.km;     % US survey mile
    mile_US = u.mi_US;          % US survey mile
    nmi = 1852*u.m;             % nautical mile
    NM = u.nmi;                 % nautical mile
    inm = u.nmi;                % nautical mile
    nauticalMile = u.nmi;
    nm_UK = 6080*u.ft;          % Imperial nautical mile
    nmile = u.nm_UK;            % Imperial nautical mile
    dataMile = 6000*u.ft;
    au = 149597870.7*u.km;      % astronomical unit
    astronomicalUnit = u.au;
    pc = 648000/pi*u.au;        % parsec
    parsec = u.pc;

    %---- reciprocal length ----

    dpt = 1/u.m;                % diopter
    diopter = u.dpt;
    R_inf = 1.0973731568508e7/u.m; % Rydberg constant
    RydbergConstant = u.R_inf;

    %---- area ----

    sqft = u.ft^2;              % square foot
    square = 100*u.sqft;        % square
    ha = 10000*u.m^2;           % hectare
    hectare = u.ha;
    a = 100*u.m^2;              % are
    are = u.a;
    ac = 43560*u.sqft;          % acre
    acre = u.ac;
    ro = 1/4*u.acre;            % rood
    rood = u.ro;
    twp = 36*u.mi^2;            % township
    township = u.twp;
    circ_mil = pi/4*u.mil^2;    % circular mil
    circularMil = u.circ_mil;
    circ_inch = pi/4*u.in^2;    % circular inch
    circularInch = u.circ_inch;
    b = 100*u.fm^2;             % barn
    barn = u.b;

    %---- volume ----

    cc = u.cm^3;                % cubic centimeter
    cubicCentimeter = u.cc;
    L = 1000*u.cc;              % liter
    l = u.L;                    % liter
    liter = u.L;
    dL = 100*u.cc;              % deciliter
    dl = u.dL;                  % deciliter
    deciliter = u.dl;
    cL = 10*u.cc;               % centiliter
    cl = u.cL;                  % centiliter
    centiliter = u.cl;
    mL = u.cc;                  % milliliter
    ml = u.mL;                  % milliliter
    milliliter = u.ml;
    uL = u.mm^3;                % microliter
    ul = u.uL;                  % microliter
    microliter = u.ul;
    kL = u.m^3;                 % kiloliter
    kl = u.kL;                  % kiloliter
    kiloliter = u.kl;
    cuin = 16.387064*u.mL;      % cubic inch
    cubicInch = u.cuin;
    FBM = u.ft^2*u.in;          % board foot
    boardFoot = u.FBM;
    gal = 231*u.cuin;           % gallon (US)
    gallon = u.gal;
    gal_UK = 4.54609*u.l;       % UK imperial gallon
    igal = u.gal_UK;            % UK imperial gallon
    quart = u.gal/4;            % US quart
    qt_UK = u.gal_UK/4;         % British imperial quart
    liq_qt = u.quart;           % US quart
    pint = u.quart/2;           % US pint
    pint_UK = u.qt_UK/2;        % British imperial pint
    liq_pt = u.pint;            % US pint
    cup = u.pint/2;             % US cup
    floz = u.cup/8;             % US fluid ounce
    fluidOunce = u.floz;        % US fluid ounce
    floz_UK = u.gal_UK/160;     % British imperial fluid ounce
    Tbls = u.floz/2;            % US tablespoon
    tablespoon = u.Tbls;        % US tablespoon
    tsp = u.Tbls/3;             % US teaspoon
    teaspoon = u.tsp;           % US teaspoon
    acft = u.acre*u.ft;         % acre-foot
    acre_foot = u.acft;
    acin = u.acre*u.in;         % acre-inch
    acre_inch = u.acin;
    bbl = 7056*u.in^3;          % US customary dry barrel
    barrel = u.bbl;
    fldr = u.floz/8;            % US customary fluid dram
    fluidDram = u.fldr;
    fldr_UK = u.floz_UK/8;      % British imperial fluid drachm (dram)
    minim = u.fldr/60;          % US customary minim
    minim_UK = u.fldr_UK/60;    % British imperial minim
    gill = 4*u.floz;            % US customary fluid gill
    gill_UK = u.gal_UK/32;      % British imperial gill

    %---- acceleration ----

    g0 = 9.80665*u.m/u.s^2;     % standard gravity
    gn = u.g0;                  % standard gravity
    g_n = u.g0;                 % standard gravity
    gee = u.g0;                 % standard gravity
    standardGravity = u.g0;
    Gal = u.cm/u.s^2;           % gal

    %---- force ----

    N = u.kg*u.m/u.s^2;         % newton
    newton = u.N;
    kN = 1e3*u.N;               % kilonewton
    kilonewton = u.kN;
    MN = 1e6*u.N;               % meganewton
    meganewton = u.MN;
    mN = 1e-3*u.N;              % millinewton
    millinewton = u.mN;
    uN = 1e-6*u.N;              % micronewton
    micronewton = u.uN;
    dyn = 1e-5*u.N;             % dyne
    dyne = u.dyn;
    lbf = 4.4482216152605*u.N;  % pound force
    poundForce = u.lbf;
    kip = 1000*u.lbf;           % kip
    kilopoundForce = u.kip;
    kgf = u.kg*u.g0;            % kilogram force
    kilogramForce = u.kgf;
    kp = u.kgf;                 % kilopond
    kilopond = u.kp;
    p = u.kp/1000;              % pond
    pond = u.p;
    sn = u.kN;                  % sthène
    sthene = u.sn;

    %---- mass ----

    kilogram = u.kg;
    kilo = u.kg;                % kilogram
    g = 1e-3*u.kg;              % gram
    gram = u.g;
    mg = 1e-3*u.gram;           % milligram
    milligram = u.mg;
    ug = 1e-6*u.gram;           % microgram
    microgram = u.ug;
    Mg = 1e6*u.gram;            % Megagram/metric tonne
    Megagram = u.Mg;
    t = 1000*u.kg;              % metric tonne
    tonne = u.t;                % metric ton
    Mt = 1e6*u.t;               % metric megatonne
    megatonne = u.Mt;
    lbm = 0.45359237*u.kg;      % pound mass
    poundMass = u.lbm;
    lb = u.lbm;                 % pound mass
    pound = u.lb;
    tn = 2000*u.lbm;            % US customary short ton
    ton = u.tn;                 % US customary short ton
    ton_UK = 2240*u.lbm;        % British imperial ton
    st = 14*u.lbm;              % stone
    stone = u.st;
    cwt = 100*u.lbm;            % US customary short hundredweight
    hundredweight = u.cwt;
    cwt_UK = 8*u.stone;         % British imperial short hundredweight
    quarter = u.cwt_UK/4;       % British imperial quarter
    slug = u.lbf/(u.ft/u.s^2);  % slug
    slinch = u.lbf/(u.in/u.s^2);
    blob = u.slinch;
    oz = u.lbm/16;              % ounce
    ounce = u.oz;
    dr = u.oz/16;               % dram
    dram = u.dr;
    gr = u.lbm/7000;            % grain
    grain = u.gr;
    ct = 200*u.mg;              % carat
    carat = u.ct;
    amu = 1.660539040e-27*u.kg; % atomic mass unit
    atomicMassUnit = u.amu;
    Da = u.amu;                 % atomic mass unit
    dalton = u.Da;
    mu = u.amu;                 % atomic mass unit
    mP = 2.176470e-8*u.kg;      % Planck mass
    PlanckMass = u.mP;
    m_e = 9.10938356e-31*u.kg;  % electron mass
    electronMass = u.m_e;
    mug = u.kgf/(u.m/u.s^2);    % metric slug
    metricSlug = u.mug;
    hyl = u.mug;                % hyl
    par = u.mug;                % par
    TMU = u.mug;                % technische Masseneinheit
    technischeMasseneinheit = u.TMU;
    glug = u.g*u.g0/(u.cm/u.s^2);

    %---- more force ----

    pdl = u.lbm*u.ft/u.s^2;     % poundal
    poundal = u.pdl;
    gf = u.gram*u.g0;           % gram force
    gramForce = u.gf;
    ozf = u.oz*u.g0;            % ounce force
    ounceForce = u.ozf;
    tonf = u.tn*u.g0;           % short ton force
    tonForce = u.tonf;

    %---- mass per length ----

    den = u.gram/(9*u.km);      % denier
    denier = u.den;
    tex = u.gram/u.km;          % tex
    dtex = u.tex/10;            % decitex
    decitex = u.dtex;

    %---- time ----

    second = u.s;
    ms = 1e-3*u.s;              % millisecond
    millisecond = u.ms;
    us = 1e-6*u.s;              % microsecond
    microsecond = u.us;
    ns = 1e-9*u.s;              % nanosecond
    nanosecond = u.ns;
    ps = 1e-12*u.s;             % picosecond
    picosecond = u.ps;
    fs = 1e-15*u.s;             % femtosecond
    femtosecond = u.fs;
    tP = 5.39116e-44*u.s;       % Planck time
    PlanckTime = u.tP;
    min = 60*u.s;               % minute
    minute = u.min;
    h = 60*u.min;               % hour
    hr = u.h;                   % hour
    hour = u.hr;
    d = 24*u.hr;                % day
    day = u.d;                  % day
    week = 7*u.day;             % week
    fortnight = 2*u.week;       % fortnight
    month_30 = 30*u.day;        % 30-day month
    yr = 365.25*u.day;          % julian year
    y = u.yr;                   % julian year
    year = u.yr;                % julian year
    year_julian = u.year;       % julian year
    year_360 = 360*u.day;       % 360-day year
    year_Tropical = 365.24219*u.day; % tropical year
    year_Gregorian = 365.2425*u.day; % gregorian year
    month = u.yr/12;            % 1/12th julian year
    flick = u.s/705600000;

    %---- frequency ----

    Hz = 1/u.s;   % hertz (NB: incompatible with angle and angular velocity)
    hertz = u.Hz;
    kHz = 1e3*u.Hz;             % kilohertz
    kilohertz = u.kHz;
    MHz = 1e6*u.Hz;             % megahertz
    megahertz = u.MHz;
    GHz = 1e9*u.Hz;             % gigahertz
    gigahertz = u.GHz;
    THz = 1e12*u.Hz;            % terahertz
    terahertz = u.THz;
    Bd = 1/u.s;                 % baud
    baud = u.Bd;

    %---- energy ----

    Nm = u.N*u.m;               % newton-meter
    newton_meter = u.Nm;
    J = u.Nm;                   % joule
    joule = u.J;
    kJ = 1e3*u.J;               % kilojoule
    kilojoule = u.kJ;
    MJ = 1e6*u.J;               % megajoule
    megajoule = u.MJ;
    GJ = 1e9*u.J;               % gigajoule
    gigajoule = u.GJ;
    mJ = 1e-3*u.J;              % millijoule
    millijoule = u.mJ;
    uJ = 1e-6*u.J;              % microjoule
    microjoule = u.uJ;
    nJ = 1e-9*u.J;              % nanojoule
    nanojoule = u.nJ;
    eV = 1.6021766208e-19*u.J;  % electronvolt
    electronvolt = u.eV;
    BTU = 1055.06*u.J;          % British thermal unit (ISO)
    Btu = u.BTU;                % British thermal unit (ISO)
    britishThermalUnit = u.Btu;
    Btu_IT = 1055.0559*u.J;     % British thermal unit (International Table)
    Btu_th = 1054.3503*u.J;     % British thermal unit (thermochemical)
    kpm = u.kp*u.m;             % kilopond-meter
    kilopond_meter = u.kpm;
    Ws = u.J;                   % watt-second
    watt_second = u.Ws;
    kWh = 3.6e6*u.J;            % kilowatt-hour
    kilowatt_hour = u.kWh;
    Wh = 3.6e3*u.J;             % watt-hour
    watt_hour = u.Wh;
    cal = 4.1868*u.J;           % calorie (International Table)
    calorie = u.cal;
    cal_IT = u.cal;             % calorie (International Table)
    cal_4 = 4.204*u.J;          % calorie (4°C)
    cal_15 = 4.1855*u.J;        % calorie (15°C)
    cal_20 = 4.182*u.J;         % calorie (20°C)
    cal_mean = 4.190*u.J; 	  % calorie (mean)
    cal_th = 4.184*u.J;         % calorie (thermochemical)
    kcal = 1e3*u.cal;           % kilocalorie
    kilocalorie = u.kcal;
    kcal_IT = 1e3*u.cal_IT;     % kilocalorie (International Table)
    Cal = u.kcal;               % large calorie / food calorie
    foodCalorie = u.Cal;
    largeCalorie = u.Cal;
    kcal_4 = 1e3*u.cal_4;       % kilocalorie (4°C)
    kcal_15 = 1e3*u.cal_15;     % kilocalorie (15°C)
    kcal_20 = 1e3*u.cal_20;     % kilocalorie (20°C)
    kcal_mean = 1e3*u.cal_mean; % kilocalorie (mean)
    kcal_th = 1e3*u.cal_th;     % kilocalorie (thermochemical)
    erg = 1e-7*u.J;             % en.wikipedia.org/wiki/Erg
    E_h = 4.359744650e-18*u.J;  % Hartree energy
    Ha = u.E_h;                 % hartree
    hartree = u.Ha;
    thm = 1e5*u.BTU;            % therm
    therm = u.thm;
    quad = 1e15*u.BTU;          % quad

    %---- temperature ----
    % For reference: °C = °K-273.15; °F = °R-459.67.

    kelvin = u.K;
    R = u.K*5/9;                % rankine (°F = °R-459.67)
    rankine = u.R;
    mK = 1e-3*u.K;              % millikelvin
    millikelvin = u.mK;
    uK = 1e-6*u.K;              % microkelvin
    microkelvin = u.uK;
    nK = 1e-9*u.K;              % nanokelvin
    nanokelvin = u.nK;
    deltaK = u.K;               % kelvin (relative temperature)
    deltadegC = u.K;            % celsius (relative, °C = °K-273.15)
    deltadegR = u.R;            % rankine (relative temperature)
    deltadegF = u.R;            % fahrenheit (relative, °F = °R-459.67)
    TP = 1.416808e32*u.K;       % Planck temperature
    PlanckTemperature = u.TP;

    %---- pressure ----

    Pa = u.N/u.m^2;             % pascal
    pascal = u.Pa;
    mPa = 1e-3*u.Pa;            % millipascal
    millipascal = u.mPa;
    uPa = 1e-6*u.Pa;            % micropascal
    micropascal = u.uPa;
    kPa = 1e3*u.Pa;             % kilopascal
    kilopascal = u.kPa;
    MPa = 1e6*u.Pa;             % megapascal
    megapascal = u.MPa;
    GPa = 1e9*u.Pa;             % gigapascal
    gigapascal = u.GPa;
    torr = 133.322*u.Pa;        % torr
    Torr = u.torr;              % torr
    mtorr = 1e-3*u.torr;        % millitorr
    millitorr = u.mtorr;
    bar = 1e5*u.Pa;             % bar
    mbar = 1e-3*u.bar;          % millibar
    millibar = u.mbar;
    kbar = 1e3*u.bar;           % kilobar
    kilobar = u.kbar;
    atm = 101325*u.Pa;          % standard atmosphere
    atmosphere = u.atm;
    standardAtmosphere = u.atm;
    at = u.kgf/u.cm^2;          % technical atmosphere
    technicalAtmosphere = u.at;
    psi = u.lbf/u.in^2;         % pound force per square inch
    poundPerSquareInch = u.psi;
    ksi = 1e3*u.psi;            % kip per square inch
    kipPerSquareInch = u.ksi;
    Msi = 1e6*u.psi;            % million pound force per square inch
    megapoundPerSquareInch = u.Msi;
    psf = u.lbf/u.ft^2;         % pound force per square foot
    poundPerSquareFoot = u.psf;
    ksf = u.kip/u.ft^2;         % kip per square foot
    kipPerSquareFoot = u.ksf;
    Ba = 0.1*u.Pa;              % barye
    barye = u.Ba;
    pz = u.kPa;                 % pièze
    pieze = u.pz;
    mmHg = 13.5951*u.kgf/u.m^2; % millimeter of mercury
    millimeterMercury = u.mmHg;
    cmHg = 10*u.mmHg;           % centimeter of mercury
    centimeterMercury = u.cmHg;
    mHg = 1e3*u.mmHg;           % meter of mercury
    meterMercury = u.mHg;
    inHg = 2.54*u.cmHg;         % inch of mercury
    inchMercury = u.inHg;
    ftHg = 12*u.inHg;           % foot of mercury
    footMercury = u.ftHg;
    mmH20 = u.kgf/u.m^2;        % millimeter of water (density 1 g/cc)
    mmAq = u.mmH20;             % millimeter of water
    millimeterWater = u.mmH20;
    cmH20 = 10*u.mmH20;         % centimeter of water
    cmAq = u.cmH20;             % centimeter of water
    centimeterWater = u.cmH20;
    mH20 = 1e3*u.mmH20;         % meter of water
    mAq = u.mH20;               % meter of water
    meterWater = u.mH20;
    inH20 = 2.54*u.cmH20;       % inch of water
    inAq = u.inH20;             % inch of water
    inchWater = u.inH20;
    wc = u.inH20;               % inch water column
    inchWaterColumn = u.wc;
    ftH20 = 12*u.inH20;         % foot of water
    ftAq = u.ftH20;             % foot of water
    footWater = u.ftH20;

    %---- viscosity ----

    St = u.cm^2/u.s;            % stokes
    stokes = u.St;
    cSt = u.St/100;             % centistokes
    centistokes = u.cSt;
    newt = u.in^2/u.s;          % newt
    P = u.Pa*u.s / 10;          % poise
    poise = u.P;
    cP = u.mPa*u.s;             % centipoise
    centipoise = u.cP;
    reyn = u.lbf*u.s/u.in^2;    % reyn

    %---- power ----

    W = u.J/u.s;                % watt
    watt = u.W;
    kW = 1e3*u.W;               % kilowatt
    kilowatt = u.kW;
    MW = 1e6*u.W;               % megawatt
    megawatt = u.MW;
    GW = 1e9*u.W;               % gigawatt
    gigawatt = u.GW;
    TW = 1e12*u.W;              % terawatt
    terawatt = u.TW;
    mW = 1e-3*u.W;              % milliwatt
    milliwatt = u.mW;
    uW = 1e-6*u.W;              % microwatt
    microwatt = u.uW;
    nW = 1e-9*u.W;              % nanowatt
    nanowatt = u.nW;
    pW = 1e-12*u.W;             % picowatt
    picowatt = u.pW;
    hp = 550*u.ft*u.lbf/u.s;    % mechanical horsepower (550 ft-lbf/s)
    horsepower = u.hp;
    HP_I = u.hp;                % mechanical horsepower (550 ft-lbf/s)
    hpE = 746*u.W;              % electrical horsepower
    HP_E = u.hpE;               % electrical horsepower
    electricalHorsepower = u.hp;
    PS = 75*u.kg*u.g0*u.m/u.s;  % metric horsepower (DIN 66036)
    HP = u.PS;                  % metric horsepower (DIN 66036)
    HP_DIN = u.PS;              % metric horsepower (DIN 66036)
    metricHorsepower = u.PS;

    %---- current ----

    amp = u.A;                  % ampere
    ampere = u.A;
    mA = 1e-3*u.A;              % milliampere
    milliampere = u.mA;
    uA = 1e-6*u.A;              % microampere
    microampere = u.uA;
    nA = 1e-9*u.A;              % nanoampere
    nanoampere = u.nA;
    pA = 1e-12*u.A;             % picoampere
    picoampere = u.pA;
    kA = 1e3*u.A;               % kiloampere
    kiloampere = u.kA;
    abA = 10*u.A;               % abampere
    abampere = u.abA;
    Bi = u.abA;                 % biot
    biot = u.Bi;

    %---- charge ----

    C = u.A*u.s;                % coulomb
    coulomb = u.C;
    mC = 1e-3*u.C;              % millicoulomb
    millicoulomb = u.mC;
    uC = 1e-6*u.C;              % microcoulomb
    microcoulomb = u.uC;
    nC = 1e-9*u.C;              % nanocoulomb
    nanocoulomb = u.nC;
    pC = 1e-12*u.C;             % picocoulomb
    picocoulomb = u.pC;
    abC = 10*u.C;               % abcoulomb
    aC = u.abC;                 % abcoulomb
    abcoulomb = u.abC;
    statC = u.dyn^(1/2)*u.cm;   % statcoulomb
    statcoulomb = u.statC;
    Fr = u.statC;               % franklin
    franklin = u.Fr;
    esu = u.statC;              % electrostatic unit of charge
    electrostaticUnitCharge = u.esu;
    e = 1.6021766208e-19*u.C;   % elementary charge
    elementaryCharge = u.e;
    mAh = u.mA*u.hr;            % milliamp-hour
    milliamp_hour = u.mAh;
    Ah = u.A*u.hr;              % amp-hour
    amp_hour = u.Ah;

    %---- voltage ----

    V = 1*u.J/u.C;              % volt
    volt = u.V;
    kV = 1e3*u.V;               % kilovolt
    kilovolt = u.kV;
    MV = 1e6*u.V;               % megavolt
    megavolt = u.MV;
    GV = 1e9*u.V;               % gigavolt
    gigavolt = u.GV;
    mV = 1e-3*u.V;              % millivolt
    millivolt = u.mV;
    uV = 1e-6*u.V;              % microvolt
    microvolt = u.uV;

    %---- resistance/conductance ----

    Ohm = u.V/u.A;              % ohm
    GOhm = 1e9*u.Ohm;           % gigaohm
    gigaohm = u.GOhm;
    MOhm = 1e6*u.Ohm;           % megaohm
    megaohm = u.MOhm;
    kOhm = 1e3*u.Ohm;           % kiloohm
    kiloohm = u.kOhm;
    mOhm = 1e-3*u.Ohm;          % milliohm
    milliohm = u.mOhm;
    uOhm = 1e-6*u.Ohm;          % microohm
    microohm = u.uOhm;
    nOhm = 1e-9*u.Ohm;          % nanoohm
    nanoohm = u.nOhm;
    abOhm = u.nOhm;             % abohm
    Z0 = 376.730313461*u.Ohm;   % characteristic impedance of vacuum
    impedanceOfVacuum = u.Z0;
    R_K = 25812.8074555*u.Ohm;  % von Klitzing constant
    vonKlitzingConstant = u.R_K;
    R_K_90 = 25812.807*u.Ohm;   % von Klitzing constant (conventional value)
    vonKlitzingConstant_conv = u.R_K_90;
    S = 1/u.Ohm;                % siemens
    siemens = u.S;
    mS = 1e-3*u.S;              % millisiemens
    millisiemens = u.mS;
    uS = 1e-6*u.S;              % microsiemens
    microsiemens = u.uS;
    nS = 1e-9*u.S;              % nanosiemens
    nanosiemens = u.nS;
    G0 = 7.7480917310e-5*u.S;   % conductance quantum
    conductanceQuantum = u.G0;

    %---- capacitance ----

    F = u.A*u.s/u.V;            % farad
    farad = u.F;
    mF = 1e-3*u.F;              % millifarad
    millifarad = u.mF;
    uF = 1e-6*u.F;              % microfarad
    microfarad = u.uF;
    nF = 1e-9*u.F;              % nanofarad
    nanofarad = u.nF;
    pF = 1e-12*u.F;             % picofarad
    picofarad = u.pF;

    %---- inductance ----

    H = u.Ohm*u.s;              % henry
    henry = u.H;
    mH = 1e-3*u.H;              % millihenry
    millihenry = u.mH;
    uH = 1e-6*u.H;              % microhenry
    microhenry = u.uH;
    nH = 1e-9*u.H;              % nanohenry
    nanohenry = u.nH;
    abH = u.nH;                 % abhenry
    abhenry = u.abH;
    kH = 1e3*u.H;               % kilohenry
    kilohenry = u.kH;
    MH = 1e6*u.H;               % megahenry
    megahenry = u.MH;
    GH = 1e9*u.H;               % gigahenry
    gigahenry = u.GH;

    %---- EM ----

    T = 1*u.N/(u.A*u.m);        % tesla
    tesla = u.T;
    Gs = 1e-4*u.T;              % gauss
    gauss = u.Gs;
    Wb = u.V*u.s;               % weber
    weber = u.Wb;
    Mx = 1e-8*u.Wb;             % maxwell
    maxwell = u.Mx;
    mWb = u.Wb/1000;            % milliweber
    milliweber = u.mWb;
    uWb = 1e-6*u.Wb;            % microweber
    microweber = u.uWb;
    nWb = 1e-9*u.Wb;            % nanoweber
    nanoweber = u.nWb;
    Oe = 250/pi*u.A/u.m;        % oersted
    oersted = u.Oe;
    Gb = 2.5/pi*u.A;            % gilbert
    gilbert = u.Gb;

    %---- non-dimensionals ----

    percent = 0.01;             % %
    pct = u.percent;            % %
    permil = 0.001;             % ‰
    permill = u.permil;         % ‰
    permille = u.permil;        % ‰
    permyriad = 1e-4;           % permyriad
    bp = u.permyriad;           % basis point
    basisPoint = u.bp;
    ppm = 1e-6;                 % part per million
    partPerMillion = u.ppm;
    ppb = 1e-9;                 % part per billion
    partPerBillion = u.ppb;
    ppt = 1e-12;                % part per trillion
    partPerTrillion = u.ppt;
    ppq = 1e-15;                % part per quadrillion
    partPerQuadrillion = u.ppq;

    %---- angles ----
    % Note: angles are dimensionless

    rad = 1;                    % radian
    radian = u.rad;
    sr = 1;                     % steradian
    steradian = u.sr;
    turn = 2*pi*u.rad;          % turn
    rev = u.turn;               % revolution = 2*pi radians
    revolution = u.rev;
    deg = u.turn/360;           % degree
    degree = u.deg;
    arcmin = u.deg/60;          % arcminute
    arcminute = u.arcmin;
    arcsec = u.arcmin/60;       % arcsecond
    arcsecond = u.arcsec;
    grad = u.turn/400;          % gradian
    gradian = u.grad;

    %---- rotational speed ----

    rpm = u.rev/u.min;          % revolution per minute
    revolutionPerMinute = u.rpm;
    rps = u.rev/u.s;            % revolution per second
    revolutionPerSecond = u.rps;

    %---- velocity ----

    mps = u.m/u.s;              % meter per second
    meterPerSecond = u.mps;
    kyne = u.cm/u.s;            % kyne
    Kyne = u.kyne;              % kyne
    fps = u.ft/u.s;             % foot per second
    footPerSecond = u.fps;
    fpm = u.ft/u.min;           % foot per minute
    footPerMinute = u.fpm;
    kt = u.nmi/u.hr;            % knot
    kn = u.kt;                  % knot
    kts = u.kt;                 % knot
    knot = u.kt;
    knot_UK = u.nm_UK/u.hr;     % British imperial knot
    KTAS = u.kt;                % knot
    nmph = u.kt;                % nautical mile per hour
    nauticalMilePerHour = u.nmph;
    kph = u.km/u.hr;            % kilometer per hour
    kmh = u.kph;                % kilometer per hour
    kilometerPerHour = u.kmh;
    mph = u.mi/u.hr;            % mile per hour
    milePerHour = u.mph;

    %---- volume flow rate ----

    cfm = u.ft^3/u.min;         % cubic foot per minute
    cubicFootPerMinute = u.cfm;
    cfs = u.ft^3/u.s;           % cubic foot per second
    cubicFootPerSecond = u.cfs;
    gpm = u.gal/u.min;          % US customary gallon per minute
    gallonPerMinute = u.gpm;
    gph = u.gal/u.hr;           % US customary gallon per hour
    gallonPerHour = u.gph;
    gpm_UK = u.gal_UK/u.min;    % British imperial gallon per minute
    lpm = u.l/u.min;            % liter per minute
    literPerMinute = u.lpm;

    %---- fuel economy ----

    l_100km = u.l/(100*u.km);   % liter per 100 km
    literPer100kilometer = u.l_100km;
    mpg = u.mi/u.gal;           % mile per gallon
    milePerGallon = u.mpg;

    %---- Luminance etc. ----

    candela = u.cd;
    asb = u.cd/u.m^2;           % apostilb
    apostilb = u.asb;
    sb = u.cd/u.cm^2;           % stilb
    stilb = u.sb;
    ph = 1e4*u.cd*u.sr/u.m^2;   % phot
    phot = u.ph;
    cp = 0.981*u.cd;            % candlepower
    candlepower = u.cp;
    lm = u.cd*u.sr;             % lumen
    lumen = u.lm;
    lx = u.lm/u.m^2;            % lux
    lux = u.lx;
    nx = 1e-3*u.lx;             % nox
    nox = u.nx;

    %---- other derived SI ----

    mole = u.mol;
    kat = u.mol/u.s;            % katal
    katal = u.kat;
    M = u.mol/u.L;              % molar
    molar = u.M;
    molarity = u.M;             % molarity
    Nms = u.N*u.m*u.s;          % newton-meter-second
    newton_meter_second = u.Nms;

    %---- radiation ----

    Gy = u.J/u.kg;              % gray
    gray = u.Gy;
    Sv = u.J/u.kg;              % sievert
    sievert = u.Sv;
    Rad =  u.Gy/100;            % absorbed radiation dose
    rem = u.Sv/100;             % roentgen equivalent man
    roentgenEquivalentMan = u.rem;
    roentgen = 2.58e-4*u.C/u.kg;% roentgen
    Ly = u.cal_th/u.cm^2;       % langley
    lan = u.Ly;                 % langley
    langley = u.lan;
    Bq = 1/u.s;                 % becquerel
    becquerel = u.Bq;
    Ci = 3.7e10*u.Bq;           % curie
    curie = u.Ci;

    %---- constants ----

    k_B = 1.38064852e-23*u.J/u.K;       % Boltzmann constant
    BoltzmannConstant = u.k_B;
    sigma_SB = 5.670367e-8*u.W/(u.m^2*u.K^4); % Stefan–Boltzmann constant
    Stefan_BoltzmannConstant = u.sigma_SB;
    h_c = 6.626070040e-34*u.J*u.s;      % Planck constant
    PlanckConstant = u.h_c;
    h_bar = u.h_c/(2*pi);               % Dirac constant
    DiracConstant = u.h_bar;
    mu_B = 9.274009994e-24*u.J/u.T;     % Bohr magneton
    BohrMagneton = u.mu_B;
    mu_N = 5.050783699e-27*u.J/u.T;     % nuclear magneton
    nuclearMagneton = u.mu_N;
    c = 299792458*u.m/u.s;              % speed of light in vacuum
    c_0 = u.c;                          % speed of light in vacuum
    lightSpeed = u.c;
    speedOfLight = u.c;
    ly = u.c*u.year;                    % light-year
    lightYear = u.ly;                   % light-year
    mu0 = pi*4e-7*u.N/u.A^2;            % vacuum permeability
    vacuumPermeability = u.mu0;
    eps0 = u.c^-2/u.mu0;                % vacuum permittivity
    vacuumPermittivity = u.eps0;
    G = 6.67408e-11*u.m^3/u.kg/u.s^2;   % gravitational constant
    gravitationalConstant = u.G;
    N_A = 6.022140857e23/u.mol;         % Avogadro constant
    NA = u.N_A;                         % Avogadro constant
    AvogadroConstant = u.N_A;
    NAh = u.N_A*u.h_c;                  % molar Planck constant
    molarPlanckConstant = u.NAh;
    M_u = u.g/u.mol;                    % molar mass constant
    molarMassConstant = u.M_u;
    K_J = 483597.8525e9*u.Hz/u.V;       % Josephson constant
    JosephsonConstant = u.K_J;
    K_J_90 = 483597.9*u.Hz/u.V;         % Josephson constant (conv. value)
    JosephsonConstant_conv = u.K_J_90;
    F_c = 96485.33289*u.C/u.mol;        % Faraday constant
    FaradayConstant = u.F_c;
    alpha = 7.2973525664e-3;            % fine-structure constant
    fine_structureConstant = u.alpha;
    SommerfeldConstant = u.alpha;
    c1 = 3.741771790e-16*u.W/u.m^2;     % first radiation constant
    firstRadiationConstant = u.c1;
    c2 = 1.43877736e-2*u.m*u.K;         % second radiation constant
    secondRadiationConstant = u.c2;
    b_prime = 5.8789238e10*u.Hz/u.K;    % Wien frequency displ. law const.
    WienFrequencyDisplacementLawConstant = u.b_prime;
    b_c = 2.8977729e-3*u.m*u.K;         % Wien wavelength displ. law const.
    WienWavelengthDisplacementLawConstant = u.b_c;
    R_air = 287.05287*u.J/u.kg/u.K;     % spec. gas const., air (ESDU 77022)
    specificGasConstant_air = u.R_air;
    R_bar = 8.3144598*u.J/u.mol/u.K;    % molar gas constant
    molarGasConstant = u.R_bar;
    radarStatuteMile = 2*u.mi/u.c;
    radarNauticalMile = 2*u.NM/u.c;
    radarDataMile = 2*u.dataMile/u.c;
    radarKilometer = 2*u.km/u.c;

    %---- digital information ----

    nibble = 4*u.bit;                   % nibble
    B = 8*u.bit;                        % byte
    byte = u.B;                         % byte
    octet = u.B;                        % octet
    kB = 1e3*u.B;                       % kilobyte
    kilobyte = u.kB;
    MB = 1e6*u.B;                       % megabyte
    megabyte = u.MB;
    GB = 1e9*u.B;                       % gigabyte
    gigabyte = u.GB;
    TB = 1e12*u.B;                      % terabyte
    terabyte = u.TB;
    PB = 1e15*u.B;                      % petabyte
    petabyte = u.PB;
    EB = 1e18*u.B;                      % exabyte
    exabyte = u.EB;
    Kibit = 2^10*u.bit;                 % kibibit
    kibibit = u.Kibit;
    KiB = 2^10*u.B;                     % kibibyte
    KB = u.KiB;                         % kibibyte
    kibibyte = u.KB;
    Mibit = 2^20*u.bit;                 % mebibit
    mebibit = u.Mibit;
    MiB = 2^20*u.B;                     % mebibyte
    mebibyte = u.MiB;
    Gibit = 2^30*u.bit;                 % gibibit
    gibibit = u.Gibit;
    GiB = 2^30*u.B;                     % gibibyte
    gibibyte = u.GiB;
    Tibit = 2^40*u.bit;                 % tebibit
    tebibit = u.Tibit;
    TiB = 2^40*u.B;                     % tebibyte
    tebibyte = u.TiB;
    Pibit = 2^50*u.bit;                 % pebibit
    pebibit = u.Pibit;
    PiB = 2^50*u.B;                     % pebibyte
    pebibyte = u.PiB;
    Eibit = 2^60*u.bit;                 % exbibit
    exbibit = u.Eibit;
    EiB = 2^60*u.B;                     % exbibyte
    exbibyte = u.EiB;
    bps = u.bit/u.s;                    % bit per second
    bitPerSecond = u.bps;
    kbps = 1e3*u.bps;                   % kilobit per second
    kilobitPerSecond = u.kbps;
    Mbps = 1e6*u.bps;                   % megabit per second
    megabitPerSecond = u.Mbps;
    Gbps = 1e9*u.bps;                   % gigabit per second
    gigabitPerSecond = u.Gbps;
    Tbps = 1e12*u.bps;                  % terabit per second
    terabitPerSecond = u.Tbps;

    %---- currency ----
    % For display purposes - not for exchange rates.
    % See also mathworks.com/matlabcentral/fileexchange/47255

    cent = u.currency/100;      % cent (currency)
    Cent = u.cent;              % cent (currency)
    pip = u.cent/100;           % pip (currency)
    USD = u.currency;           % currency
    EUR = u.currency;           % currency
    GBP = u.currency;           % currency
    JPY = u.currency;           % currency
    AUD = u.currency;           % currency
    CAD = u.currency;           % currency
    CHF = u.currency;           % currency
    CNY = u.currency;           % currency
    dollar = u.currency;        % currency
    franc = u.currency;         % currency

    %---- used by Matlab's symunit but not here ----
    % gg -gauge
    % land - league
    % ha_US - US survey hectare
    % molecule
    % HP_UK - British imperial horsepower
    % PS_SAE - net horsepower (SAE J1349)
    % PS_DIN - horsepower (DIN 70020)
    % dry volumes
end

%% METHODS
methods
    %% Plotting and display:
    function disp(o)
        try    
            dispdisp(o); % mathworks.com/matlabcentral/fileexchange/48637
        catch
            % builtin('disp',o);
            
            url = 'http://www.mathworks.com/matlabcentral/fileexchange/48637/';
            dlCmd = sprintf('matlab:unzip(websave(tempname,''%s%s''),pwd);u',...
                url,'?download=true');
            
            warning('The function <a href="%s">%s</a> %s\n%s',...
                'www.mathworks.com/matlabcentral/fileexchange/48637',...
                'dispdisp',...
                'is recommended for display of physical units.',...
                ['<a href="' dlCmd ...
                '">Direct download of dispdisp into current directory</a>']);
        end
    end
end
end

%% Processing base units.
function U = buildCoreUnits(baseUnitSystem)
coreBaseNames = {'m' 'kg' 's' 'A' 'K' 'mol' 'cd' 'bit' 'currency'};

if ischar(baseUnitSystem) && strcmpi('none',baseUnitSystem)
    % Use normal variables - not DimVars - if baseUnitSystem is 'none'.
    U = cell2struct(num2cell(ones(size(coreBaseNames))),coreBaseNames,2);
    return
end

validateattributes(baseUnitSystem,{'cell'},{'size',[9,2]},'u','baseUnitSystem');

if ~iscellstr(u.baseUnitSystem(:,1))
    error('First column of baseUnitSystem must be type char.')
end

baseValues = baseUnitSystem(:,2);
if ~all(cellfun('isclass', baseValues, 'double'))
    error('Second column of baseUnitSystem must contain doubles.')
end

expos = eye(numel(coreBaseNames));
for i = 1:numel(baseValues)
    U.(coreBaseNames{i}) = DimVar(expos(i,:),baseValues{i});
end
end

%%
%   Original inspiration for this tool by Rob deCarvalho.
%     http://www.mathworks.com/matlabcentral/fileexchange/authors/22148
%     http://www.mathworks.com/matlabcentral/fileexchange/10070