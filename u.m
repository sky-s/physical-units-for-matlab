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
%   When displaying variables with units or using them in plot, etc., the units
%   used for display will be, in order if available and valid: 1) per-variable
%   preferred display units, 2) units listed in displayUnits, 3) a combination
%   of fundamental base units (mass, length, time, temperature, ...). To set (or
%   clear) custom display units for a given variable, see the scd function. To
%   customize the displayUnits list, see displayUnits. For more advanced
%   customization of the base units themselves, see baseUnitSystem.
%
%   Display customization is set by calls to displayUnits and/or baseUnitSystem
%   (either function files or variables in the base workspace). Tailor
%   preferences for a specific project by defining these variables at the top of
%   a script (before any units are called) or placing unique versions of the
%   files in a project's directory. Be sure to clear the class when changing
%   projects or else the old customizations will remain in effect.
%
%   Some MATLAB functions won't accept variables with physical units (DimVars).
%   Most of the time displayingvalue, which returns value in terms of preferred
%   display units, will be the appropriate tool, but there is also double, which
%   returns the value in terms of base units, and u2num.
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
%   Example 3: Custom display.
%       fieldSize = 3*u.sqkm
%       % Muliplies and divides remove any per-variable custom display units.
%       % It's nice to display in the units that make sense for that variable.
%       rate = 3.7*u.acre/u.day
%       rate = scd(rate,'sqm/hr')
%       timeNeeded = fieldSize/rate
%       timeNeeded = scd(timeNeeded,'month')
% 
%   See also displayUnits, baseUnitSystem, scd, clear, displayingvalue,
%   DimVar.double, u2num, str2u, symunit,
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

    m           = scd(u.coreUnits.m,'m') % meter
    kg          = scd(u.coreUnits.kg,'kg') % kilogram
    s           = scd(u.coreUnits.s,'s') % second
    A           = scd(u.coreUnits.A,'A') % ampere
    K           = scd(u.coreUnits.K,'K') % kelvin (°C = °K-273.15)
    mol         = scd(u.coreUnits.mol,'mol') % mole
    cd          = scd(u.coreUnits.cd,'cd') % candela
    bit         = scd(u.coreUnits.bit,'bit') % bit
    currency    = scd(u.coreUnits.currency,'currency') % currency
    unit        = scd(u.coreUnits.unit,'unit') % user unit

    %% Derived units list:
    % References:
    % http://physics.nist.gov/cuu/Constants/index.html
    % http://www.translatorscafe.com/unit-converter
    % http://en.wikipedia.org
    % http://www.efunda.com/units/index.cfm
    
    %---- length ----

    meter = scd(u.m,'meter') 
    km = scd(1e3*u.m,'km') % kilometer
    kilometer = scd(u.km,'kilometer') 
    dm = scd(1e-1*u.m,'dm') % decimeter
    decimeter = scd(u.dm,'decimeter') 
    cm = scd(1e-2*u.m,'cm') % centimeter
    centimeter = scd(u.cm,'centimeter') 
    mm = scd(1e-3*u.m,'mm') % millimeter
    millimeter = scd(u.mm,'millimeter') 
    um = scd(1e-6*u.m,'um') % micrometer
    micrometer = scd(u.um,'micrometer') 
    micron = scd(u.um,'micron') % micron
    nm = scd(1e-9*u.m,'nm') % nanometer
    nanometer = scd(u.nm,'nanometer') 
    pm = scd(1e-12*u.m,'pm') % picometer
    picometer = scd(u.pm,'picometer') 
    fm = scd(1e-15*u.m,'fm') % femtometer
    femtometer = scd(u.fm,'femtometer') 
    fermi = scd(u.fm,'fermi') % fermi
    Ao = scd(1e-10*u.m,'Ao') % ångström
    ang = scd(u.Ao,'ang') % ångström
    angstrom = scd(u.ang,'angstrom') 
    angstroem = scd(u.ang,'angstroem') 
    a0 = scd(0.52917721067e-10*u.m,'a0') % Bohr radius
    a_0 = scd(u.a0,'a_0') % Bohr radius
    BohrRadius = scd(u.a0,'BohrRadius') 
    lP = scd(1.616229e-35*u.m,'lP') % Planck length
    PlanckLength = scd(u.lP,'PlanckLength') 
    xu = scd(1.0021e-13*u.m,'xu') % x unit
    xUnit = scd(u.xu,'xUnit') 
    xu_Cu = scd(1.00207697e-13*u.m,'xu_Cu') % x unit (copper)
    xUnit_copper = scd(u.xu_Cu,'xUnit_copper') 
    xu_Mo = scd(1.00209952e-13*u.m,'xu_Mo') % x unit (molybdenum)
    xUnit_molybdenum = scd(u.xu_Mo,'xUnit_molybdenum') 
    in = scd(2.54*u.cm,'in') % inch
    inch = scd(u.in,'inch') 
    mil = scd(1e-3*u.in,'mil') % mil
    line = scd(u.in/10,'line') % line
    hand = scd(4*u.in,'hand') % hand
    span = scd(9*u.in,'span') % span
    smoot = scd(67*u.in,'smoot') % smoot
    ft = scd(12*u.in,'ft') % foot
    foot = scd(u.ft,'foot') 
    ft_US = scd(1200/3937*u.m,'ft_US') % US survey foot
    foot_US = scd(u.ft_US,'foot_US') % US survey foot
    kft = scd(1e3*u.ft,'kft') % kilofoot
    kilofoot = scd(u.kft,'kilofoot') 
    FL = scd(100*u.ft,'FL') % flight level
    flightLevel = scd(u.FL,'flightLevel') 
    yd = scd(3*u.ft,'yd') % yard
    yard = scd(u.yd,'yard') 
    ftm = scd(6*u.ft,'ftm') % fathom
    fathom = scd(u.ftm,'fathom') 
    li = scd(0.66*u.ft,'li') % link
    link = scd(u.li,'link') 
    rod = scd(5.5*u.yd,'rod') % rod
    ch = scd(66*u.ft,'ch') % chain
    chain = scd(u.ch,'chain') 
    fur = scd(220*u.yd,'fur') % furlong
    furlong = scd(u.fur,'furlong') 
    mi = scd(5280*u.ft,'mi') % mile
    mile = scd(u.mi,'mile') 
    mi_US = scd(6336/3937*u.km,'mi_US') % US survey mile
    mile_US = scd(u.mi_US,'mile_US') % US survey mile
    nmi = scd(1852*u.m,'nmi') % nautical mile
    NM = scd(u.nmi,'NM') % nautical mile
    inm = scd(u.nmi,'inm') % nautical mile
    nauticalMile = scd(u.nmi,'nauticalMile') 
    nm_UK = scd(6080*u.ft,'nm_UK') % Imperial nautical mile
    nmile = scd(u.nm_UK,'nmile') % Imperial nautical mile
    dataMile = scd(6000*u.ft,'dataMile') 
    au = scd(149597870.7*u.km,'au') % astronomical unit
    astronomicalUnit = scd(u.au,'astronomicalUnit') 
    pc = scd(648000/pi*u.au,'pc') % parsec
    parsec = scd(u.pc,'parsec') 

    %---- reciprocal length ----

    dpt = scd(1/u.m,'dpt') % diopter
    diopter = scd(u.dpt,'diopter') 
    R_inf = scd(1.0973731568508e7/u.m,'R_inf') % Rydberg constant
    RydbergConstant = scd(u.R_inf,'RydbergConstant') 

    %---- area ----

    ft2 = scd(u.ft^2,'ft²') % square foot
    sqft = scd(u.ft2,'sqft') % square foot
    square = scd(100*u.sqft,'square') % square
    ha = scd(10000*u.m^2,'ha') % hectare
    hectare = scd(u.ha,'hectare') 
    a = scd(100*u.m^2,'a') % are
    are = scd(u.a,'are') 
    ac = scd(43560*u.sqft,'ac') % acre
    acre = scd(u.ac,'acre') 
    ro = scd(1/4*u.acre,'ro') % rood
    rood = scd(u.ro,'rood') 
    twp = scd(36*u.mi^2,'twp') % township
    township = scd(u.twp,'township') 
    circ_mil = scd(pi/4*u.mil^2,'circ_mil') % circular mil
    circularMil = scd(u.circ_mil,'circularMil') 
    circ_inch = scd(pi/4*u.in^2,'circ_inch') % circular inch
    circularInch = scd(u.circ_inch,'circularInch') 
    b = scd(100*u.fm^2,'b') % barn
    barn = scd(u.b,'barn') 
    sqin = scd(u.in^2,'sqin') % square inch
    squareInch = scd(u.sqin,'squareInch')
    sqmil = scd(u.mil^2,'sqmil') % square mil
    squareMil = scd(u.sqmil,'squareMil')
    sqmi = scd(u.mi^2,'sqmi') % square mile
    squareMile = scd(u.sqmi,'squareMile')
    sqnmi = scd(u.nmi^2,'sqnmi') % square nautical mile
    squareNauticalMile = scd(u.sqnmi,'squareNauticalMile')
    m2 = scd(u.m^2,'m²') % square meter
    sqm = scd(u.m^2,'sqm') % square meter
    squareMeter = scd(u.sqm,'squareMeter')
    sqkm = scd(u.km^2,'sqkm') % square kilometer
    squareKilometer = scd(u.sqkm,'squareKilometer')
    sqcm = scd(u.cm^2,'sqcm') % square centimeter
    squareCentimeter = scd(u.sqcm,'squareCentimeter')
    sqmm = scd(u.mm^2,'sqmm') % square millimeter
    squareMillimeter = scd(u.sqmm,'squareMillimeter')
    sqdm = scd(u.dm^2,'sqdm') % square decimeter
    squareDecimeter = scd(u.sqdm,'squareDecimeter')

    %---- volume ----

    m3 = scd(u.m^3,'m³') % cubic meter
    cc = scd(u.cm^3,'cc') % cubic centimeter
    cubicCentimeter = scd(u.cc,'cubicCentimeter') 
    L = scd(1000*u.cc,'L') % liter
    l = scd(u.L,'l') % liter
    liter = scd(u.L,'liter') 
    dL = scd(100*u.cc,'dL') % deciliter
    dl = scd(u.dL,'dl') % deciliter
    deciliter = scd(u.dl,'deciliter') 
    cL = scd(10*u.cc,'cL') % centiliter
    cl = scd(u.cL,'cl') % centiliter
    centiliter = scd(u.cl,'centiliter') 
    mL = scd(u.cc,'mL') % milliliter
    ml = scd(u.mL,'ml') % milliliter
    milliliter = scd(u.ml,'milliliter') 
    uL = scd(u.mm^3,'uL') % microliter
    ul = scd(u.uL,'ul') % microliter
    microliter = scd(u.ul,'microliter') 
    kL = scd(u.m^3,'kL') % kiloliter
    kl = scd(u.kL,'kl') % kiloliter
    kiloliter = scd(u.kl,'kiloliter') 
    cuin = scd(16.387064*u.mL,'cuin') % cubic inch
    cubicInch = scd(u.cuin,'cubicInch')
    ft3 = scd(u.ft^3,'ft³') % cubic foot
    FBM = scd(u.sqft*u.in,'FBM') % board foot
    boardFoot = scd(u.FBM,'boardFoot') 
    gal = scd(231*u.cuin,'gal') % gallon (US)
    gallon = scd(u.gal,'gallon') 
    gal_UK = scd(4.54609*u.l,'gal_UK') % UK imperial gallon
    igal = scd(u.gal_UK,'igal') % UK imperial gallon
    quart = scd(u.gal/4,'quart') % US quart
    qt_UK = scd(u.gal_UK/4,'qt_UK') % British imperial quart
    liq_qt = scd(u.quart,'liq_qt') % US quart
    pint = scd(u.quart/2,'pint') % US pint
    pint_UK = scd(u.qt_UK/2,'pint_UK') % British imperial pint
    liq_pt = scd(u.pint,'liq_pt') % US pint
    cup = scd(u.pint/2,'cup') % US cup
    floz = scd(u.cup/8,'floz') % US fluid ounce
    fluidOunce = scd(u.floz,'fluidOunce') % US fluid ounce
    floz_UK = scd(u.gal_UK/160,'floz_UK') % British imperial fluid ounce
    Tbls = scd(u.floz/2,'Tbls') % US tablespoon
    tablespoon = scd(u.Tbls,'tablespoon') % US tablespoon
    tsp = scd(u.Tbls/3,'tsp') % US teaspoon
    teaspoon = scd(u.tsp,'teaspoon') % US teaspoon
    acft = scd(u.acre*u.ft,'acft') % acre-foot
    acre_foot = scd(u.acft,'acre_foot') 
    acin = scd(u.acre*u.in,'acin') % acre-inch
    acre_inch = scd(u.acin,'acre_inch') 
    bbl = scd(7056*u.in^3,'bbl') % US customary dry barrel
    barrel = scd(u.bbl,'barrel') 
    fldr = scd(u.floz/8,'fldr') % US customary fluid dram
    fluidDram = scd(u.fldr,'fluidDram') 
    fldr_UK = scd(u.floz_UK/8,'fldr_UK') % British imperial fluid drachm (dram)
    minim = scd(u.fldr/60,'minim') % US customary minim
    minim_UK = scd(u.fldr_UK/60,'minim_UK') % British imperial minim
    gill = scd(4*u.floz,'gill') % US customary fluid gill
    gill_UK = scd(u.gal_UK/32,'gill_UK') % British imperial gill

    %---- acceleration ----

    g0 = scd(9.80665*u.m/u.s^2,'g0') % standard gravity
    g_0 = scd(u.g0,'g_0') % standard gravity
    gn = scd(u.g0,'gn') % standard gravity
    g_n = scd(u.g0,'g_n') % standard gravity
    gee = scd(u.g0,'gee') % standard gravity
    standardGravity = scd(u.g0,'standardGravity') 
    Gal = scd(u.cm/u.s^2,'Gal') % gal

    %---- force ----

    N = scd(u.kg*u.m/u.s^2,'N') % newton
    newton = scd(u.N,'newton') 
    kN = scd(1e3*u.N,'kN') % kilonewton
    kilonewton = scd(u.kN,'kilonewton') 
    MN = scd(1e6*u.N,'MN') % meganewton
    meganewton = scd(u.MN,'meganewton') 
    mN = scd(1e-3*u.N,'mN') % millinewton
    millinewton = scd(u.mN,'millinewton') 
    uN = scd(1e-6*u.N,'uN') % micronewton
    micronewton = scd(u.uN,'micronewton') 
    dyn = scd(1e-5*u.N,'dyn') % dyne
    dyne = scd(u.dyn,'dyne') 
    lbf = scd(4.4482216152605*u.N,'lbf') % pound force
    lb_f = scd(u.lbf,'lb_f') % pound force
    poundForce = scd(u.lbf,'poundForce') 
    kip = scd(1000*u.lbf,'kip') % kip
    kilopoundForce = scd(u.kip,'kilopoundForce') 
    kgf = scd(u.kg*u.g0,'kgf') % kilogram force
    kg_f = scd(u.kgf,'kg_f') % kilogram force
    kilogramForce = scd(u.kgf,'kilogramForce') 
    kp = scd(u.kgf,'kp') % kilopond
    kilopond = scd(u.kp,'kilopond') 
    p = scd(u.kp/1000,'p') % pond
    pond = scd(u.p,'pond') 
    sn = scd(u.kN,'sn') % sthène
    sthene = scd(u.sn,'sthene') 

    %---- mass ----

    kilogram = scd(u.kg,'kilogram') 
    kilo = scd(u.kg,'kilo') % kilogram
    g = scd(1e-3*u.kg,'g') % gram
    gram = scd(u.g,'gram') 
    mg = scd(1e-3*u.gram,'mg') % milligram
    milligram = scd(u.mg,'milligram') 
    ug = scd(1e-6*u.gram,'ug') % microgram
    microgram = scd(u.ug,'microgram') 
    Mg = scd(1e6*u.gram,'Mg') % Megagram/metric tonne
    Megagram = scd(u.Mg,'Megagram') 
    t = scd(1000*u.kg,'t') % metric tonne
    tonne = scd(u.t,'tonne') % metric ton
    Mt = scd(1e6*u.t,'Mt') % metric megatonne
    megatonne = scd(u.Mt,'megatonne') 
    lbm = scd(0.45359237*u.kg,'lbm') % pound mass
    lb_m = scd(u.lbm,'lb_m') % pound mass
    poundMass = scd(u.lbm,'poundMass') 
    lb = scd(u.lbm,'lb') % pound mass
    pound = scd(u.lb,'pound') 
    tn = scd(2000*u.lbm,'tn') % US customary short ton
    ton = scd(u.tn,'ton') % US customary short ton
    ton_UK = scd(2240*u.lbm,'ton_UK') % British imperial ton
    st = scd(14*u.lbm,'st') % stone
    stone = scd(u.st,'stone') 
    cwt = scd(100*u.lbm,'cwt') % US customary short hundredweight
    hundredweight = scd(u.cwt,'hundredweight') 
    cwt_UK = scd(8*u.stone,'cwt_UK') % British imperial short hundredweight
    quarter = scd(u.cwt_UK/4,'quarter') % British imperial quarter
    slug = scd(u.lbf/(u.ft/u.s^2),'slug') % slug
    slinch = scd(u.lbf/(u.in/u.s^2),'slinch') 
    blob = scd(u.slinch,'blob') 
    oz = scd(u.lbm/16,'oz') % ounce
    ounce = scd(u.oz,'ounce') 
    dr = scd(u.oz/16,'dr') % dram
    dram = scd(u.dr,'dram') 
    gr = scd(u.lbm/7000,'gr') % grain
    grain = scd(u.gr,'grain') 
    ct = scd(200*u.mg,'ct') % carat
    carat = scd(u.ct,'carat') 
    amu = scd(1.660539040e-27*u.kg,'amu') % atomic mass unit
    atomicMassUnit = scd(u.amu,'atomicMassUnit') 
    Da = scd(u.amu,'Da') % atomic mass unit
    dalton = scd(u.Da,'dalton') 
    mu = scd(u.amu,'mu') % atomic mass unit
    mP = scd(2.176470e-8*u.kg,'mP') % Planck mass
    PlanckMass = scd(u.mP,'PlanckMass') 
    m_e = scd(9.10938356e-31*u.kg,'m_e') % electron mass
    electronMass = scd(u.m_e,'electronMass') 
    mug = scd(u.kgf/(u.m/u.s^2),'mug') % metric slug
    metricSlug = scd(u.mug,'metricSlug') 
    hyl = scd(u.mug,'hyl') % hyl
    par = scd(u.mug,'par') % par
    TMU = scd(u.mug,'TMU') % technische Masseneinheit
    technischeMasseneinheit = scd(u.TMU,'technischeMasseneinheit') 
    glug = scd(u.g*u.g0/(u.cm/u.s^2),'glug') 

    %---- more force ----

    pdl = scd(u.lbm*u.ft/u.s^2,'pdl') % poundal
    poundal = scd(u.pdl,'poundal') 
    gf = scd(u.gram*u.g0,'gf') % gram force
    g_f = scd(u.gf,'g_f') %gram force
    gramForce = scd(u.gf,'gramForce') 
    ozf = scd(u.oz*u.g0,'ozf') % ounce force
    oz_f = scd(u.ozf,'oz_f') % ounce force
    ounceForce = scd(u.ozf,'ounceForce') 
    tonf = scd(u.tn*u.g0,'tonf') % short ton force
    ton_f = scd(u.tonf,'ton_f') % short ton force
    tonForce = scd(u.tonf,'tonForce') 

    %---- mass per length ----

    den = scd(u.gram/(9*u.km),'den') % denier
    denier = scd(u.den,'denier') 
    tex = scd(u.gram/u.km,'tex') % tex
    dtex = scd(u.tex/10,'dtex') % decitex
    decitex = scd(u.dtex,'decitex') 

    %---- time ----

    second = scd(u.s,'second') 
    sec = scd(u.s,'sec') % second
    ms = scd(1e-3*u.s,'ms') % millisecond
    millisecond = scd(u.ms,'millisecond') 
    us = scd(1e-6*u.s,'us') % microsecond
    microsecond = scd(u.us,'microsecond') 
    ns = scd(1e-9*u.s,'ns') % nanosecond
    nanosecond = scd(u.ns,'nanosecond') 
    ps = scd(1e-12*u.s,'ps') % picosecond
    picosecond = scd(u.ps,'picosecond') 
    fs = scd(1e-15*u.s,'fs') % femtosecond
    femtosecond = scd(u.fs,'femtosecond') 
    tP = scd(5.39116e-44*u.s,'tP') % Planck time
    PlanckTime = scd(u.tP,'PlanckTime') 
    min = scd(60*u.s,'min') % minute
    minute = scd(u.min,'minute') 
    h = scd(60*u.min,'h') % hour
    hr = scd(u.h,'hr') % hour
    hrs = scd(u.h,'hrs')
    hour = scd(u.hr,'hour') 
    d = scd(24*u.hr,'d') % day
    day = scd(u.d,'day') % day
    days = scd(u.d,'days')
    week = scd(7*u.day,'week') % week
    fortnight = scd(2*u.week,'fortnight') % fortnight
    month_30 = scd(30*u.day,'month_30') % 30-day month
    yr = scd(365.25*u.day,'yr') % julian year
    y = scd(u.yr,'y') % julian year
    year = scd(u.yr,'year') % julian year
    yrs = scd(u.yr,'yrs') % julian year
    year_julian = scd(u.year,'year_julian') % julian year
    year_360 = scd(360*u.day,'year_360') % 360-day year
    year_Tropical = scd(365.24219*u.day,'year_Tropical') % tropical year
    year_Gregorian = scd(365.2425*u.day,'year_Gregorian') % gregorian year
    month = scd(u.yr/12,'month') % 1/12th julian year
    flick = scd(u.s/705600000,'flick') 

    %---- frequency ----

    Hz = scd(1/u.s,'Hz') % hertz (NB: incompatible with angle and angular velocity)
    hertz = scd(u.Hz,'hertz') 
    kHz = scd(1e3*u.Hz,'kHz') % kilohertz
    kilohertz = scd(u.kHz,'kilohertz') 
    MHz = scd(1e6*u.Hz,'MHz') % megahertz
    megahertz = scd(u.MHz,'megahertz') 
    GHz = scd(1e9*u.Hz,'GHz') % gigahertz
    gigahertz = scd(u.GHz,'gigahertz') 
    THz = scd(1e12*u.Hz,'THz') % terahertz
    terahertz = scd(u.THz,'terahertz') 
    Bd = scd(1/u.s,'Bd') % baud
    baud = scd(u.Bd,'baud') 

    %---- energy ----

    Nm = scd(u.N*u.m,'Nm') % newton-meter
    newton_meter = scd(u.Nm,'newton_meter') 
    J = scd(u.Nm,'J') % joule
    joule = scd(u.J,'joule') 
    kJ = scd(1e3*u.J,'kJ') % kilojoule
    kilojoule = scd(u.kJ,'kilojoule') 
    MJ = scd(1e6*u.J,'MJ') % megajoule
    megajoule = scd(u.MJ,'megajoule') 
    GJ = scd(1e9*u.J,'GJ') % gigajoule
    gigajoule = scd(u.GJ,'gigajoule') 
    mJ = scd(1e-3*u.J,'mJ') % millijoule
    millijoule = scd(u.mJ,'millijoule') 
    uJ = scd(1e-6*u.J,'uJ') % microjoule
    microjoule = scd(u.uJ,'microjoule') 
    nJ = scd(1e-9*u.J,'nJ') % nanojoule
    nanojoule = scd(u.nJ,'nanojoule') 
    eV = scd(1.6021766208e-19*u.J,'eV') % electronvolt
    electronvolt = scd(u.eV,'electronvolt') 
    BTU = scd(1055.06*u.J,'BTU') % British thermal unit (ISO)
    Btu = scd(u.BTU,'Btu') % British thermal unit (ISO)
    britishThermalUnit = scd(u.Btu,'britishThermalUnit') 
    Btu_IT = scd(1055.0559*u.J,'Btu_IT') % British thermal unit (International Table)
    Btu_th = scd(1054.3503*u.J,'Btu_th') % British thermal unit (thermochemical)
    kpm = scd(u.kp*u.m,'kpm') % kilopond-meter
    kilopond_meter = scd(u.kpm,'kilopond_meter') 
    Ws = scd(u.J,'Ws') % watt-second
    watt_second = scd(u.Ws,'watt_second') 
    kWh = scd(3.6e6*u.J,'kWh') % kilowatt-hour
    kilowatt_hour = scd(u.kWh,'kilowatt_hour') 
    Wh = scd(3.6e3*u.J,'Wh') % watt-hour
    watt_hour = scd(u.Wh,'watt_hour') 
    cal = scd(4.1868*u.J,'cal') % calorie (International Table)
    calorie = scd(u.cal,'calorie') 
    cal_IT = scd(u.cal,'cal_IT') % calorie (International Table)
    cal_4 = scd(4.204*u.J,'cal_4') % calorie (4°C)
    cal_15 = scd(4.1855*u.J,'cal_15') % calorie (15°C)
    cal_20 = scd(4.182*u.J,'cal_20') % calorie (20°C)
    cal_mean = scd(4.190*u.J,'cal_mean') % calorie (mean)
    cal_th = scd(4.184*u.J,'cal_th') % calorie (thermochemical)
    tonTnt = scd(u.cal_th*1e9,'tonTnt') % ton of TNT
    kcal = scd(1e3*u.cal,'kcal') % kilocalorie
    kilocalorie = scd(u.kcal,'kilocalorie') 
    kcal_IT = scd(1e3*u.cal_IT,'kcal_IT') % kilocalorie (International Table)
    Cal = scd(u.kcal,'Cal') % large calorie / food calorie
    foodCalorie = scd(u.Cal,'foodCalorie') 
    largeCalorie = scd(u.Cal,'largeCalorie') 
    kcal_4 = scd(1e3*u.cal_4,'kcal_4') % kilocalorie (4°C)
    kcal_15 = scd(1e3*u.cal_15,'kcal_15') % kilocalorie (15°C)
    kcal_20 = scd(1e3*u.cal_20,'kcal_20') % kilocalorie (20°C)
    kcal_mean = scd(1e3*u.cal_mean,'kcal_mean') % kilocalorie (mean)
    kcal_th = scd(1e3*u.cal_th,'kcal_th') % kilocalorie (thermochemical)
    erg = scd(1e-7*u.J,'erg') % en.wikipedia.org/wiki/Erg
    E_h = scd(4.359744650e-18*u.J,'E_h') % Hartree energy
    Ha = scd(u.E_h,'Ha') % hartree
    hartree = scd(u.Ha,'hartree') 
    thm = scd(1e5*u.BTU,'thm') % therm
    therm = scd(u.thm,'therm') 
    quad = scd(1e15*u.BTU,'quad') % quad

    %---- temperature ----
    % For reference: °C = °K-273.15; °F = °R-459.67.

    kelvin = scd(u.K,'kelvin') 
    R = scd(u.K*5/9,'R') % rankine (°F = °R-459.67)
    rankine = scd(u.R,'rankine') 
    mK = scd(1e-3*u.K,'mK') % millikelvin
    millikelvin = scd(u.mK,'millikelvin') 
    uK = scd(1e-6*u.K,'uK') % microkelvin
    microkelvin = scd(u.uK,'microkelvin') 
    nK = scd(1e-9*u.K,'nK') % nanokelvin
    nanokelvin = scd(u.nK,'nanokelvin') 
    deltaK = scd(u.K,'deltaK') % kelvin (relative temperature)
    deltadegC = scd(u.K,'deltadegC') % celsius (relative, °C = °K-273.15)
    deltadegR = scd(u.R,'deltadegR') % rankine (relative temperature)
    deltadegF = scd(u.R,'deltadegF') % fahrenheit (relative, °F = °R-459.67)
    TP = scd(1.416808e32*u.K,'TP') % Planck temperature
    PlanckTemperature = scd(u.TP,'PlanckTemperature') 

    %---- pressure ----

    Pa = scd(u.N/u.sqm,'Pa') % pascal
    pascal = scd(u.Pa,'pascal') 
    mPa = scd(1e-3*u.Pa,'mPa') % millipascal
    millipascal = scd(u.mPa,'millipascal') 
    uPa = scd(1e-6*u.Pa,'uPa') % micropascal
    micropascal = scd(u.uPa,'micropascal') 
    kPa = scd(1e3*u.Pa,'kPa') % kilopascal
    kilopascal = scd(u.kPa,'kilopascal') 
    MPa = scd(1e6*u.Pa,'MPa') % megapascal
    megapascal = scd(u.MPa,'megapascal') 
    GPa = scd(1e9*u.Pa,'GPa') % gigapascal
    gigapascal = scd(u.GPa,'gigapascal') 
    torr = scd(133.322*u.Pa,'torr') % torr
    Torr = scd(u.torr,'Torr') % torr
    mtorr = scd(1e-3*u.torr,'mtorr') % millitorr
    millitorr = scd(u.mtorr,'millitorr') 
    bar = scd(1e5*u.Pa,'bar') % bar
    mbar = scd(1e-3*u.bar,'mbar') % millibar
    millibar = scd(u.mbar,'millibar') 
    kbar = scd(1e3*u.bar,'kbar') % kilobar
    kilobar = scd(u.kbar,'kilobar') 
    atm = scd(101325*u.Pa,'atm') % standard atmosphere
    atmosphere = scd(u.atm,'atmosphere') 
    standardAtmosphere = scd(u.atm,'standardAtmosphere') 
    at = scd(u.kgf/u.sqcm,'at') % technical atmosphere
    technicalAtmosphere = scd(u.at,'technicalAtmosphere') 
    psi = scd(u.lbf/u.sqin,'psi') % pound force per square inch
    poundPerSquareInch = scd(u.psi,'poundPerSquareInch') 
    ksi = scd(1e3*u.psi,'ksi') % kip per square inch
    kipPerSquareInch = scd(u.ksi,'kipPerSquareInch') 
    Msi = scd(1e6*u.psi,'Msi') % million pound force per square inch
    megapoundPerSquareInch = scd(u.Msi,'megapoundPerSquareInch') 
    psf = scd(u.lbf/u.sqft,'psf') % pound force per square foot
    poundPerSquareFoot = scd(u.psf,'poundPerSquareFoot') 
    ksf = scd(u.kip/u.sqft,'ksf') % kip per square foot
    kipPerSquareFoot = scd(u.ksf,'kipPerSquareFoot') 
    Ba = scd(0.1*u.Pa,'Ba') % barye
    barye = scd(u.Ba,'barye') 
    pz = scd(u.kPa,'pz') % pièze
    pieze = scd(u.pz,'pieze') 
    mmHg = scd(13.5951*u.kgf/u.sqm,'mmHg') % millimeter of mercury
    millimeterMercury = scd(u.mmHg,'millimeterMercury') 
    cmHg = scd(10*u.mmHg,'cmHg') % centimeter of mercury
    centimeterMercury = scd(u.cmHg,'centimeterMercury') 
    mHg = scd(1e3*u.mmHg,'mHg') % meter of mercury
    meterMercury = scd(u.mHg,'meterMercury') 
    inHg = scd(2.54*u.cmHg,'inHg') % inch of mercury
    inchMercury = scd(u.inHg,'inchMercury') 
    ftHg = scd(12*u.inHg,'ftHg') % foot of mercury
    footMercury = scd(u.ftHg,'footMercury') 
    mmH20 = scd(u.kgf/u.sqm,'mmH20') % millimeter of water (density 1 g/cc)
    mmAq = scd(u.mmH20,'mmAq') % millimeter of water
    millimeterWater = scd(u.mmH20,'millimeterWater') 
    cmH20 = scd(10*u.mmH20,'cmH20') % centimeter of water
    cmAq = scd(u.cmH20,'cmAq') % centimeter of water
    centimeterWater = scd(u.cmH20,'centimeterWater') 
    mH20 = scd(1e3*u.mmH20,'mH20') % meter of water
    mAq = scd(u.mH20,'mAq') % meter of water
    meterWater = scd(u.mH20,'meterWater') 
    inH20 = scd(2.54*u.cmH20,'inH20') % inch of water
    inAq = scd(u.inH20,'inAq') % inch of water
    inchWater = scd(u.inH20,'inchWater') 
    wc = scd(u.inH20,'wc') % inch water column
    inchWaterColumn = scd(u.wc,'inchWaterColumn') 
    ftH20 = scd(12*u.inH20,'ftH20') % foot of water
    ftAq = scd(u.ftH20,'ftAq') % foot of water
    footWater = scd(u.ftH20,'footWater') 

    %---- viscosity ----

    St = scd(u.sqcm/u.s,'St') % stokes
    stokes = scd(u.St,'stokes') 
    cSt = scd(u.St/100,'cSt') % centistokes
    centistokes = scd(u.cSt,'centistokes') 
    newt = scd(u.sqin/u.s,'newt') % newt
    P = scd(u.Pa*u.s / 10,'P') % poise
    poise = scd(u.P,'poise') 
    cP = scd(u.mPa*u.s,'cP') % centipoise
    centipoise = scd(u.cP,'centipoise') 
    reyn = scd(u.lbf*u.s/u.sqin,'reyn') % reyn

    %---- power ----

    W = scd(u.J/u.s,'W') % watt
    watt = scd(u.W,'watt') 
    kW = scd(1e3*u.W,'kW') % kilowatt
    kilowatt = scd(u.kW,'kilowatt') 
    MW = scd(1e6*u.W,'MW') % megawatt
    megawatt = scd(u.MW,'megawatt') 
    GW = scd(1e9*u.W,'GW') % gigawatt
    gigawatt = scd(u.GW,'gigawatt') 
    TW = scd(1e12*u.W,'TW') % terawatt
    terawatt = scd(u.TW,'terawatt') 
    mW = scd(1e-3*u.W,'mW') % milliwatt
    milliwatt = scd(u.mW,'milliwatt') 
    uW = scd(1e-6*u.W,'uW') % microwatt
    microwatt = scd(u.uW,'microwatt') 
    nW = scd(1e-9*u.W,'nW') % nanowatt
    nanowatt = scd(u.nW,'nanowatt') 
    pW = scd(1e-12*u.W,'pW') % picowatt
    picowatt = scd(u.pW,'picowatt') 
    hp = scd(550*u.ft*u.lbf/u.s,'hp') % mechanical horsepower (550 ft-lbf/s)
    horsepower = scd(u.hp,'horsepower') 
    HP_I = scd(u.hp,'HP_I') % mechanical horsepower (550 ft-lbf/s)
    hpE = scd(746*u.W,'hpE') % electrical horsepower
    HP_E = scd(u.hpE,'HP_E') % electrical horsepower
    electricalHorsepower = scd(u.hp,'electricalHorsepower') 
    PS = scd(75*u.kg*u.g0*u.m/u.s,'PS') % metric horsepower (DIN 66036)
    HP = scd(u.PS,'HP') % metric horsepower (DIN 66036)
    HP_DIN = scd(u.PS,'HP_DIN') % metric horsepower (DIN 66036)
    metricHorsepower = scd(u.PS,'metricHorsepower') 

    %---- current ----

    amp = scd(u.A,'amp') % ampere
    ampere = scd(u.A,'ampere') 
    mA = scd(1e-3*u.A,'mA') % milliampere
    milliampere = scd(u.mA,'milliampere') 
    uA = scd(1e-6*u.A,'uA') % microampere
    microampere = scd(u.uA,'microampere') 
    nA = scd(1e-9*u.A,'nA') % nanoampere
    nanoampere = scd(u.nA,'nanoampere') 
    pA = scd(1e-12*u.A,'pA') % picoampere
    picoampere = scd(u.pA,'picoampere') 
    kA = scd(1e3*u.A,'kA') % kiloampere
    kiloampere = scd(u.kA,'kiloampere') 
    abA = scd(10*u.A,'abA') % abampere
    abampere = scd(u.abA,'abampere') 
    Bi = scd(u.abA,'Bi') % biot
    biot = scd(u.Bi,'biot') 

    %---- charge ----

    C = scd(u.A*u.s,'C') % coulomb
    coulomb = scd(u.C,'coulomb') 
    mC = scd(1e-3*u.C,'mC') % millicoulomb
    millicoulomb = scd(u.mC,'millicoulomb') 
    uC = scd(1e-6*u.C,'uC') % microcoulomb
    microcoulomb = scd(u.uC,'microcoulomb') 
    nC = scd(1e-9*u.C,'nC') % nanocoulomb
    nanocoulomb = scd(u.nC,'nanocoulomb') 
    pC = scd(1e-12*u.C,'pC') % picocoulomb
    picocoulomb = scd(u.pC,'picocoulomb') 
    abC = scd(10*u.C,'abC') % abcoulomb
    aC = scd(u.abC,'aC') % abcoulomb
    abcoulomb = scd(u.abC,'abcoulomb') 
    statC = scd(u.dyn^(1/2)*u.cm,'statC') % statcoulomb
    statcoulomb = scd(u.statC,'statcoulomb') 
    Fr = scd(u.statC,'Fr') % franklin
    franklin = scd(u.Fr,'franklin') 
    esu = scd(u.statC,'esu') % electrostatic unit of charge
    electrostaticUnitCharge = scd(u.esu,'electrostaticUnitCharge') 
    e = scd(1.6021766208e-19*u.C,'e') % elementary charge
    elementaryCharge = scd(u.e,'elementaryCharge') 
    mAh = scd(u.mA*u.hr,'mAh') % milliamp-hour
    milliamp_hour = scd(u.mAh,'milliamp_hour') 
    Ah = scd(u.A*u.hr,'Ah') % amp-hour
    amp_hour = scd(u.Ah,'amp_hour') 

    %---- voltage ----

    V = scd(1*u.J/u.C,'V') % volt
    volt = scd(u.V,'volt') 
    kV = scd(1e3*u.V,'kV') % kilovolt
    kilovolt = scd(u.kV,'kilovolt') 
    MV = scd(1e6*u.V,'MV') % megavolt
    megavolt = scd(u.MV,'megavolt') 
    GV = scd(1e9*u.V,'GV') % gigavolt
    gigavolt = scd(u.GV,'gigavolt') 
    mV = scd(1e-3*u.V,'mV') % millivolt
    millivolt = scd(u.mV,'millivolt') 
    uV = scd(1e-6*u.V,'uV') % microvolt
    microvolt = scd(u.uV,'microvolt') 

    %---- resistance/conductance ----

    Ohm = scd(u.V/u.A,'Ohm') % ohm
    GOhm = scd(1e9*u.Ohm,'GOhm') % gigaohm
    gigaohm = scd(u.GOhm,'gigaohm') 
    MOhm = scd(1e6*u.Ohm,'MOhm') % megaohm
    megaohm = scd(u.MOhm,'megaohm') 
    kOhm = scd(1e3*u.Ohm,'kOhm') % kiloohm
    kiloohm = scd(u.kOhm,'kiloohm') 
    mOhm = scd(1e-3*u.Ohm,'mOhm') % milliohm
    milliohm = scd(u.mOhm,'milliohm') 
    uOhm = scd(1e-6*u.Ohm,'uOhm') % microohm
    microohm = scd(u.uOhm,'microohm') 
    nOhm = scd(1e-9*u.Ohm,'nOhm') % nanoohm
    nanoohm = scd(u.nOhm,'nanoohm') 
    abOhm = scd(u.nOhm,'abOhm') % abohm
    Z0 = scd(376.730313461*u.Ohm,'Z0') % characteristic impedance of vacuum
    impedanceOfVacuum = scd(u.Z0,'impedanceOfVacuum') 
    R_K = scd(25812.8074555*u.Ohm,'R_K') % von Klitzing constant
    vonKlitzingConstant = scd(u.R_K,'vonKlitzingConstant') 
    R_K_90 = scd(25812.807*u.Ohm,'R_K_90') % von Klitzing constant (conventional value)
    vonKlitzingConstant_conv = scd(u.R_K_90,'vonKlitzingConstant_conv') 
    S = scd(1/u.Ohm,'S') % siemens
    siemens = scd(u.S,'siemens') 
    mS = scd(1e-3*u.S,'mS') % millisiemens
    millisiemens = scd(u.mS,'millisiemens') 
    uS = scd(1e-6*u.S,'uS') % microsiemens
    microsiemens = scd(u.uS,'microsiemens') 
    nS = scd(1e-9*u.S,'nS') % nanosiemens
    nanosiemens = scd(u.nS,'nanosiemens') 
    G0 = scd(7.7480917310e-5*u.S,'G0') % conductance quantum
    conductanceQuantum = scd(u.G0,'conductanceQuantum') 

    %---- capacitance ----

    F = scd(u.A*u.s/u.V,'F') % farad
    farad = scd(u.F,'farad') 
    mF = scd(1e-3*u.F,'mF') % millifarad
    millifarad = scd(u.mF,'millifarad') 
    uF = scd(1e-6*u.F,'uF') % microfarad
    microfarad = scd(u.uF,'microfarad') 
    nF = scd(1e-9*u.F,'nF') % nanofarad
    nanofarad = scd(u.nF,'nanofarad') 
    pF = scd(1e-12*u.F,'pF') % picofarad
    picofarad = scd(u.pF,'picofarad') 

    %---- inductance ----

    H = scd(u.Ohm*u.s,'H') % henry
    henry = scd(u.H,'henry') 
    mH = scd(1e-3*u.H,'mH') % millihenry
    millihenry = scd(u.mH,'millihenry') 
    uH = scd(1e-6*u.H,'uH') % microhenry
    microhenry = scd(u.uH,'microhenry') 
    nH = scd(1e-9*u.H,'nH') % nanohenry
    nanohenry = scd(u.nH,'nanohenry') 
    abH = scd(u.nH,'abH') % abhenry
    abhenry = scd(u.abH,'abhenry') 
    kH = scd(1e3*u.H,'kH') % kilohenry
    kilohenry = scd(u.kH,'kilohenry') 
    MH = scd(1e6*u.H,'MH') % megahenry
    megahenry = scd(u.MH,'megahenry') 
    GH = scd(1e9*u.H,'GH') % gigahenry
    gigahenry = scd(u.GH,'gigahenry') 

    %---- EM ----

    T = scd(1*u.N/(u.A*u.m),'T') % tesla
    tesla = scd(u.T,'tesla') 
    Gs = scd(1e-4*u.T,'Gs') % gauss
    gauss = scd(u.Gs,'gauss') 
    Wb = scd(u.V*u.s,'Wb') % weber
    weber = scd(u.Wb,'weber') 
    Mx = scd(1e-8*u.Wb,'Mx') % maxwell
    maxwell = scd(u.Mx,'maxwell') 
    mWb = scd(u.Wb/1000,'mWb') % milliweber
    milliweber = scd(u.mWb,'milliweber') 
    uWb = scd(1e-6*u.Wb,'uWb') % microweber
    microweber = scd(u.uWb,'microweber') 
    nWb = scd(1e-9*u.Wb,'nWb') % nanoweber
    nanoweber = scd(u.nWb,'nanoweber') 
    Oe = scd(250/pi*u.A/u.m,'Oe') % oersted
    oersted = scd(u.Oe,'oersted') 
    Gb = scd(2.5/pi*u.A,'Gb') % gilbert
    gilbert = scd(u.Gb,'gilbert') 

    %---- non-dimensionals ----

    percent = 0.01 % %
    pct = u.percent % %
    permil = 0.001 % ‰
    permill = u.permil % ‰
    permille = u.permil % ‰
    permyriad = 1e-4 % permyriad
    bp = u.permyriad % basis point
    basisPoint = u.bp
    ppm = 1e-6 % part per million
    partPerMillion = u.ppm 
    ppb = 1e-9 % part per billion
    partPerBillion = u.ppb
    ppt = 1e-12 % part per trillion
    partPerTrillion = u.ppt
    ppq = 1e-15 % part per quadrillion
    partPerQuadrillion = u.ppq 
    
    %---- angles ----
    % Note: angles are dimensionless

    rad = 1 % radian
    radian = u.rad
    sr = 1 % steradian
    steradian = u.sr
    turn = 2*pi*u.rad % turn
    rev = u.turn % revolution = 2*pi radians
    revolution = u.rev
    deg = u.turn/360 % degree
    degree = u.deg
    arcmin = u.deg/60 % arcminute
    arcminute = u.arcmin
    arcsec = u.arcmin/60 % arcsecond
    arcsecond = u.arcsec
    grad = u.turn/400 % gradian
    gradian = u.grad
    
    %---- rotational speed ----

    rpm = scd(u.rev/u.min,'rpm') % revolution per minute
    revolutionPerMinute = scd(u.rpm,'revolutionPerMinute') 
    rps = scd(u.rev/u.s,'rps') % revolution per second
    revolutionPerSecond = scd(u.rps,'revolutionPerSecond') 

    %---- velocity ----

    mps = scd(u.m/u.s,'mps') % meter per second
    meterPerSecond = scd(u.mps,'meterPerSecond') 
    kyne = scd(u.cm/u.s,'kyne') % kyne
    Kyne = scd(u.kyne,'Kyne') % kyne
    fps = scd(u.ft/u.s,'fps') % foot per second
    footPerSecond = scd(u.fps,'footPerSecond') 
    fpm = scd(u.ft/u.min,'fpm') % foot per minute
    footPerMinute = scd(u.fpm,'footPerMinute') 
    kt = scd(u.nmi/u.hr,'kt') % knot
    kn = scd(u.kt,'kn') % knot
    kts = scd(u.kt,'kts') % knot
    knot = scd(u.kt,'knot') 
    knot_UK = scd(u.nm_UK/u.hr,'knot_UK') % British imperial knot
    KTAS = scd(u.kt,'KTAS') % knot
    nmph = scd(u.kt,'nmph') % nautical mile per hour
    nauticalMilePerHour = scd(u.nmph,'nauticalMilePerHour') 
    kph = scd(u.km/u.hr,'kph') % kilometer per hour
    kmh = scd(u.kph,'kmh') % kilometer per hour
    kps = scd(u.km/u.s,'kps') % kilometer per second
    kilometerPerHour = scd(u.kmh,'kilometerPerHour') 
    mph = scd(u.mi/u.hr,'mph') % mile per hour
    milePerHour = scd(u.mph,'milePerHour') 

    %---- volume flow rate ----

    cfm = scd(u.ft^3/u.min,'cfm') % cubic foot per minute
    cubicFootPerMinute = scd(u.cfm,'cubicFootPerMinute') 
    cfs = scd(u.ft^3/u.s,'cfs') % cubic foot per second
    cubicFootPerSecond = scd(u.cfs,'cubicFootPerSecond') 
    gpm = scd(u.gal/u.min,'gpm') % US customary gallon per minute
    gallonPerMinute = scd(u.gpm,'gallonPerMinute') 
    gph = scd(u.gal/u.hr,'gph') % US customary gallon per hour
    gallonPerHour = scd(u.gph,'gallonPerHour') 
    gpm_UK = scd(u.gal_UK/u.min,'gpm_UK') % British imperial gallon per minute
    lpm = scd(u.l/u.min,'lpm') % liter per minute
    literPerMinute = scd(u.lpm,'literPerMinute') 

    %---- fuel economy ----

    l_100km = scd(u.l/(100*u.km),'l_100km') % liter per 100 km
    literPer100kilometer = scd(u.l_100km,'literPer100kilometer') 
    mpg = scd(u.mi/u.gal,'mpg') % mile per gallon
    milePerGallon = scd(u.mpg,'milePerGallon') 

    %---- Luminance etc. ----

    candela = scd(u.cd,'candela') 
    asb = scd(u.cd/u.sqm,'asb') % apostilb
    apostilb = scd(u.asb,'apostilb') 
    sb = scd(u.cd/u.sqcm,'sb') % stilb
    stilb = scd(u.sb,'stilb') 
    ph = scd(1e4*u.cd*u.sr/u.sqm,'ph') % phot
    phot = scd(u.ph,'phot') 
    cp = scd(0.981*u.cd,'cp') % candlepower
    candlepower = scd(u.cp,'candlepower') 
    lm = scd(u.cd*u.sr,'lm') % lumen
    lumen = scd(u.lm,'lumen') 
    lx = scd(u.lm/u.sqm,'lx') % lux
    lux = scd(u.lx,'lux') 
    nx = scd(1e-3*u.lx,'nx') % nox
    nox = scd(u.nx,'nox') 

    %---- other derived SI ----

    mole = scd(u.mol,'mole') 
    kat = scd(u.mol/u.s,'kat') % katal
    katal = scd(u.kat,'katal') 
    M = scd(u.mol/u.L,'M') % molar
    molar = scd(u.M,'molar') 
    molarity = scd(u.M,'molarity') % molarity
    Nms = scd(u.N*u.m*u.s,'Nms') % newton-meter-second
    newton_meter_second = scd(u.Nms,'newton_meter_second') 

    %---- radiation ----

    Gy = scd(u.J/u.kg,'Gy') % gray
    gray = scd(u.Gy,'gray') 
    Sv = scd(u.J/u.kg,'Sv') % sievert
    sievert = scd(u.Sv,'sievert') 
    Rad = scd( u.Gy/100,'Rad') % absorbed radiation dose
    rem = scd(u.Sv/100,'rem') % roentgen equivalent man
    roentgenEquivalentMan = scd(u.rem,'roentgenEquivalentMan') 
    roentgen = scd(2.58e-4*u.C/u.kg,'roentgen') % roentgen
    Ly = scd(u.cal_th/u.sqcm,'Ly') % langley
    lan = scd(u.Ly,'lan') % langley
    langley = scd(u.lan,'langley') 
    Bq = scd(1/u.s,'Bq') % becquerel
    becquerel = scd(u.Bq,'becquerel') 
    Ci = scd(3.7e10*u.Bq,'Ci') % curie
    curie = scd(u.Ci,'curie') 

    %---- constants ----
    
    i = 1i
    j = 1j
    pi = pi % Archimedes' constant ?
    tau = 2*pi
    phi = (1 + sqrt(5))/2 % golden ratio
    EulersNumber = exp(1) % ("e" is reserved for elementary charge)

    k_B = scd(1.38064852e-23*u.J/u.K,'k_B') % Boltzmann constant
    BoltzmannConstant = scd(u.k_B,'BoltzmannConstant') 
    sigma_SB = scd(5.670367e-8*u.W/(u.sqm*u.K^4),'sigma_SB') % Stefan–Boltzmann constant
    Stefan_BoltzmannConstant = scd(u.sigma_SB,'Stefan_BoltzmannConstant') 
    h_c = scd(6.626070040e-34*u.J*u.s,'h_c') % Planck constant
    PlanckConstant = scd(u.h_c,'PlanckConstant') 
    h_bar = scd(u.h_c/(2*pi),'h_bar') % Dirac constant
    DiracConstant = scd(u.h_bar,'DiracConstant') 
    mu_B = scd(9.274009994e-24*u.J/u.T,'mu_B') % Bohr magneton
    BohrMagneton = scd(u.mu_B,'BohrMagneton') 
    mu_N = scd(5.050783699e-27*u.J/u.T,'mu_N') % nuclear magneton
    nuclearMagneton = scd(u.mu_N,'nuclearMagneton') 
    c = scd(299792458*u.m/u.s,'c') % speed of light in vacuum
    c_0 = scd(u.c,'c_0') % speed of light in vacuum
    lightSpeed = scd(u.c,'lightSpeed') 
    speedOfLight = scd(u.c,'speedOfLight') 
    ly = scd(u.c*u.year,'ly') % light-year
    lightYear = scd(u.ly,'lightYear') % light-year
    mu0 = scd(pi*4e-7*u.N/u.A^2,'mu0') % vacuum permeability
    vacuumPermeability = scd(u.mu0,'vacuumPermeability') 
    eps0 = scd(u.c^-2/u.mu0,'eps0') % vacuum permittivity
    vacuumPermittivity = scd(u.eps0,'vacuumPermittivity') 
    G = scd(6.67408e-11*u.m^3/u.kg/u.s^2,'G') % gravitational constant
    gravitationalConstant = scd(u.G,'gravitationalConstant') 
    N_A = scd(6.022140857e23/u.mol,'N_A') % Avogadro constant
    NA = scd(u.N_A,'NA') % Avogadro constant
    AvogadroConstant = scd(u.N_A,'AvogadroConstant') 
    NAh = scd(u.N_A*u.h_c,'NAh') % molar Planck constant
    molarPlanckConstant = scd(u.NAh,'molarPlanckConstant') 
    M_u = scd(u.g/u.mol,'M_u') % molar mass constant
    molarMassConstant = scd(u.M_u,'molarMassConstant') 
    K_J = scd(483597.8525e9*u.Hz/u.V,'K_J') % Josephson constant
    JosephsonConstant = scd(u.K_J,'JosephsonConstant') 
    K_J_90 = scd(483597.9*u.Hz/u.V,'K_J_90') % Josephson constant (conv. value)
    JosephsonConstant_conv = scd(u.K_J_90,'JosephsonConstant_conv') 
    F_c = scd(96485.33289*u.C/u.mol,'F_c') % Faraday constant
    FaradayConstant = scd(u.F_c,'FaradayConstant') 
    alpha = 7.2973525664e-3 % fine-structure constant
    fine_structureConstant = u.alpha
    SommerfeldConstant = u.alpha
    c1 = scd(3.741771790e-16*u.W/u.sqm,'c1') % first radiation constant
    firstRadiationConstant = scd(u.c1,'firstRadiationConstant') 
    c2 = scd(1.43877736e-2*u.m*u.K,'c2') % second radiation constant
    secondRadiationConstant = scd(u.c2,'secondRadiationConstant') 
    b_prime = scd(5.8789238e10*u.Hz/u.K,'b_prime') % Wien frequency displ. law const.
    WienFrequencyDisplacementLawConstant = scd(u.b_prime,'WienFrequencyDisplacementLawConstant') 
    b_c = scd(2.8977729e-3*u.m*u.K,'b_c') % Wien wavelength displ. law const.
    WienWavelengthDisplacementLawConstant = scd(u.b_c,'WienWavelengthDisplacementLawConstant') 
    R_air = scd(287.05287*u.J/u.kg/u.K,'R_air') % spec. gas const., air (ESDU 77022)
    specificGasConstant_air = scd(u.R_air,'specificGasConstant_air') 
    R_bar = scd(8.3144598*u.J/u.mol/u.K,'R_bar') % molar gas constant
    molarGasConstant = scd(u.R_bar,'molarGasConstant') 
    radarStatuteMile = scd(2*u.mi/u.c,'radarStatuteMile') 
    radarNauticalMile = scd(2*u.NM/u.c,'radarNauticalMile') 
    radarDataMile = scd(2*u.dataMile/u.c,'radarDataMile') 
    radarKilometer = scd(2*u.km/u.c,'radarKilometer') 

    %---- digital information ----

    nibble = scd(4*u.bit,'nibble') 
    B = scd(8*u.bit,'B') % byte
    byte = scd(u.B,'byte') 
    octet = scd(u.B,'octet') % octet
    kB = scd(1e3*u.B,'kB') % kilobyte
    kilobyte = scd(u.kB,'kilobyte') 
    MB = scd(1e6*u.B,'MB') % megabyte
    megabyte = scd(u.MB,'megabyte') 
    GB = scd(1e9*u.B,'GB') % gigabyte
    gigabyte = scd(u.GB,'gigabyte') 
    TB = scd(1e12*u.B,'TB') % terabyte
    terabyte = scd(u.TB,'terabyte') 
    PB = scd(1e15*u.B,'PB') % petabyte
    petabyte = scd(u.PB,'petabyte') 
    EB = scd(1e18*u.B,'EB') % exabyte
    exabyte = scd(u.EB,'exabyte') 
    Kibit = scd(2^10*u.bit,'Kibit') % kibibit
    kibibit = scd(u.Kibit,'kibibit') 
    KiB = scd(2^10*u.B,'KiB') % kibibyte
    KB = scd(u.KiB,'KB') % kibibyte
    kibibyte = scd(u.KB,'kibibyte') 
    Mibit = scd(2^20*u.bit,'Mibit') % mebibit
    mebibit = scd(u.Mibit,'mebibit') 
    MiB = scd(2^20*u.B,'MiB') % mebibyte
    mebibyte = scd(u.MiB,'mebibyte') 
    Gibit = scd(2^30*u.bit,'Gibit') % gibibit
    gibibit = scd(u.Gibit,'gibibit') 
    GiB = scd(2^30*u.B,'GiB') % gibibyte
    gibibyte = scd(u.GiB,'gibibyte') 
    Tibit = scd(2^40*u.bit,'Tibit') % tebibit
    tebibit = scd(u.Tibit,'tebibit') 
    TiB = scd(2^40*u.B,'TiB') % tebibyte
    tebibyte = scd(u.TiB,'tebibyte') 
    Pibit = scd(2^50*u.bit,'Pibit') % pebibit
    pebibit = scd(u.Pibit,'pebibit') 
    PiB = scd(2^50*u.B,'PiB') % pebibyte
    pebibyte = scd(u.PiB,'pebibyte') 
    Eibit = scd(2^60*u.bit,'Eibit') % exbibit
    exbibit = scd(u.Eibit,'exbibit') 
    EiB = scd(2^60*u.B,'EiB') % exbibyte
    exbibyte = scd(u.EiB,'exbibyte') 
    bps = scd(u.bit/u.s,'bps') % bit per second
    bitPerSecond = scd(u.bps,'bitPerSecond') 
    kbps = scd(1e3*u.bps,'kbps') % kilobit per second
    kilobitPerSecond = scd(u.kbps,'kilobitPerSecond') 
    Mbps = scd(1e6*u.bps,'Mbps') % megabit per second
    megabitPerSecond = scd(u.Mbps,'megabitPerSecond') 
    Gbps = scd(1e9*u.bps,'Gbps') % gigabit per second
    gigabitPerSecond = scd(u.Gbps,'gigabitPerSecond') 
    Tbps = scd(1e12*u.bps,'Tbps') % terabit per second
    terabitPerSecond = scd(u.Tbps,'terabitPerSecond') 

    %---- currency ----
    % For display purposes - not for exchange rates.
    % See also mathworks.com/matlabcentral/fileexchange/47255

    cent = scd(u.currency/100,'cent') % cent (currency)
    Cent = scd(u.cent,'Cent') % cent (currency)
    pip = scd(u.cent/100,'pip') % pip (currency)
    USD = scd(u.currency,'USD') % currency
    EUR = scd(u.currency,'EUR') % currency
    GBP = scd(u.currency,'GBP') % currency
    JPY = scd(u.currency,'JPY') % currency
    AUD = scd(u.currency,'AUD') % currency
    CAD = scd(u.currency,'CAD') % currency
    CHF = scd(u.currency,'CHF') % currency
    CNY = scd(u.currency,'CNY') % currency
    dollar = scd(u.currency,'dollar') % currency
    franc = scd(u.currency,'franc') % currency

    %---- used by Matlab's symunit but not here ----
    % gg - gauge
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
        f = fieldnames(o);
        for iField = 1:length(f)
            thisField = u.(f{iField});
            if isa(thisField,'DimVar')
                thisField = scd(thisField);
            end
            uDisplayStruct.(f{iField}) = thisField;
        end
                
        try    
            dispdisp(uDisplayStruct);
            % mathworks.com/matlabcentral/fileexchange/48637
        catch
            builtin('disp',uDisplayStruct);
            
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
coreBaseNames = {'m' 'kg' 's' 'A' 'K' 'mol' 'cd' 'bit' 'currency' 'unit'};

if ischar(baseUnitSystem) && strcmpi('none',baseUnitSystem)
    % Use normal variables - not DimVars - if baseUnitSystem is 'none'.
    U = cell2struct(num2cell(ones(size(coreBaseNames))),coreBaseNames,2);
    return
end

validateattributes(baseUnitSystem,{'cell'},{'size',[10,2]},...
    'u','baseUnitSystem');

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