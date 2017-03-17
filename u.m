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
    %   A variable with unusual units will display as a combination of fundamental
    %   base units (mass, length, time, temperature, ...). By default these are the
    %   SI base units of kg, m, s, A, K, mol, cd, and currency. However, the display
    %   of any unit, whether a base unit or a derived unit, can be customized. For
    %   example, instead of displaying power as [m^2][kg]/s^3, customizing base
    %   units could make it display in terms of [ft^2][lbm]/s^3. Or, if desired for
    %   a given project, all variables with units corresponding to "power" could
    %   instead be displayed in terms of a single derived unit, e.g. W or hp.
    %
    %   Display customization is set by whatever myUnits.m file is highest on the
    %   MATLAB search path, so a unique myUnits.m can be placed in a project's
    %   directory to tailor preferred display units for that project. Be sure to
    %   clear the class when changing projects or else the old myUnits.m will remain
    %   in effect.
    %
    %   Some MATLAB functions won't accept variables with physical units. See u2num.
    %
    %   Example 1: Shaft power.
    %       rotationSpeed = 2500 * u.rpm;
    %       torque = 95 * str2u('ft-lbf');  % Use alternate string-based definition.
    %       power = rotationSpeed * torque; % Returns variable with units of power.
    %       horsePower = power / u.hp;      % Convert/cancel units.
    %
    %   Example 2: Unit conversion.
    %       100 * u.acre/u.ha;  % Convert 100 acres to hectares.
    %       u.st/u.kg;          % Retern conversion factor for stone to kilos.
    %
    %   See also myUnits, clear, str2u, u2num, symunit,
    %     dispdisp - http://www.mathworks.com/matlabcentral/fileexchange/48637.
    
    %   Copyright Sky Sartorius
    %   www.mathworks.com/matlabcentral/fileexchange/authors/101715
    %   github.com/sky-s/physical-units-for-matlab
    
    properties (Constant = true)
        %% Initialization:
        
        % Establishes base unit system and preferences based on myUnits.m
        % highest on the MATLAB search path.
        U = myUnits()
        
        %% Base units:
        
        m           = u.U.m         % meter
        kg          = u.U.kg        % kilogram
        s           = u.U.s         % second
        A           = u.U.A         % ampere
        K           = u.U.K         % kelvin (°C = °K-273.15)
        mol         = u.U.mol       % mole
        cd          = u.U.cd        % candela
        bit         = u.U.bit       % bit
        currency    = u.U.currency  % currency
        
        %% Derived units:
        %---- length ----
        
        km = 1e3*u.m;               % kilometer
        dm = 1e-1*u.m;              % decimeter
        cm = 1e-2*u.m;              % centimeter
        mm = 1e-3*u.m;              % millimeter
        um = 1e-6*u.m;              % micrometer
        micron = u.um;              % micron
        nm = 1e-9*u.m;              % nanometer
        pm = 1e-12*u.m;             % picometer
        fm = 1e-15*u.m;             % femtometer
        fermi = u.fm;               % fermi
        Ao = 1e-10*u.m;             % ångström
        ang = u.Ao;                 % ångström
        a0 = 0.52917721067e-10*u.m; % Bohr radius
        a_0 = u.a0;                 % Bohr radius
        lP = 1.616229e-35*u.m;      % Planck length
        xu = 1.0021e-13*u.m;        % x unit
        xu_Cu = 1.00207697e-13*u.m; % x unit (copper)
        xu_Mo = 1.00209952e-13*u.m; % x unit (molybdenum)
        in = 2.54*u.cm;             % inch
        mil = 1e-3*u.in;            % mil
        line = u.in/10;             % line
        hand = 4*u.in;              % hand
        span = 9*u.in;              % span
        smoot = 67*u.in;            % smoot
        ft = 12*u.in;               % foot
        ft_US = 1200/3937*u.m;      % US survey foot
        kft = 1e3*u.ft;             % kilofoot
        yd = 3*u.ft;                % yard
        fathom = 6*u.ft;            % fathom
        ftm = u.fathom;             % fathom
        li = 0.66*u.ft;             % link
        rod = 5.5*u.yd;             % rod
        ch = 66*u.ft;               % chain
        furlong = 220*u.yd;         % furlong
        fur = u.furlong;            % furlong
        mi = 5280*u.ft;             % mile
        mi_US = 6336/3937*u.km;     % US survey mile
        nmi = 1852*u.m;             % nautical mile
        NM = u.nmi;                 % nautical mile
        inm = u.nmi;                % nautical mile
        nm_UK = 6080*u.ft;          % Imperial nautical mile
        nmile = u.nm_UK;            % Imperial nautical mile
        au = 149597870.7*u.km;      % astronomical unit
        pc = 648000/pi*u.au;        % parsec
        
        %---- reciprocal length ----
        
        dpt = 1/u.m;                % diopter
        R_inf = 1.0973731568508e7/u.m; % Rydberg constant
        
        %---- area ----
        
        square = 100*u.ft^2;        % square
        ha = 10000*u.m^2;           % hectare
        hectare = u.ha;             % hectare
        a = 100*u.m^2;              % are
        are = u.a;                  % are
        ac = 43560*u.ft^2;          % acre
        acre = u.ac;                % acre
        ro = 1/4*u.acre;            % rood
        twp = 36*u.mi^2;            % township
        circ_mil = pi/4*u.mil^2;    % circular mil
        circ_inch = pi/4*u.in^2;    % circular inch
        b = 100*u.fm^2;             % barn
        barn = u.b;                 % barn
        
        %---- volume ----
        
        cc = u.cm^3;                % cubic centimeter
        L = 1000*u.cc;              % liter
        l = u.L;                    % liter
        dL = 100*u.cc;              % deciliter
        dl = u.dL;                  % deciliter
        cL = 10*u.cc;               % centiliter
        cl = u.cL;                  % centiliter
        mL = u.cc;                  % milliliter
        ml = u.mL;                  % milliliter
        uL = u.mm^3;                % microliter
        ul = u.uL;                  % microliter
        kL = u.m^3;                 % kiloliter
        kl = u.kL;                  % kiloliter
        cuin = 16.387064*u.mL;      % cubic inch
        FBM = u.ft^2*u.in;          % board foot
        gal = 231*u.cuin;           % US gallon
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
        floz_UK = u.gal_UK/160;     % British imperial fluid ounce
        Tbls = u.floz/2;            % US tablespoon
        tsp = u.Tbls/3;             % US teaspoon
        acft = u.acre*u.ft;         % acre-foot
        acin = u.acre*u.in;         % acre-inch
        barrel = 7056*u.in^3;       % US customary dry barrel
        bbl = u.barrel;             % US customary dry barrel
        fldr = u.floz/8;            % US customary fluid dram
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
        Gal = u.cm/u.s^2;           % gal
        
        %---- force ----
        
        N = u.kg*u.m/u.s^2;         % newton
        kN = 1000*u.N;              % kilonewton
        mN = 1e-3*u.N;              % millinewton
        dyn = 1e-5*u.N;             % dyne
        lbf = 4.4482216152605*u.N;  % pound force
        kip = 1000*u.lbf;           % kip
        kgf = u.kg*u.g0;            % kilogram force
        kp = u.kgf;                 % kilopond
        p = u.kp/1000;              % pond
        sn = u.kN;                  % sthène
        
        %---- mass ----
        
        gram = 1e-3*u.kg;           % gram
        g = u.gram;                 % gram
        mg = 1e-3*u.gram;           % milligram
        ug = 1e-6*u.gram;           % microgram
        Mg = 1e6*u.gram;            % Megagram/metric ton
        t = 1000*u.kg;              % metric ton
        tonne = u.t;                % metric ton
        Mt = 1e6*u.t;               % metric megaton
        lbm = 0.45359237*u.kg;      % pound mass
        lb = u.lbm;                 % pound mass
        tn = 2000*u.lbm;            % US customary short ton
        ton_UK = 2240*u.lbm;        % British imperial ton
        st = 14*u.lbm;              % stone
        stone = u.st;               % stone
        cwt = 100*u.lbm;            % US customary short hundredweight
        cwt_UK = 8*u.stone;         % British imperial short hundredweight
        quarter = u.cwt_UK/4;       % British imperial quarter
        slug = u.lbf/(u.ft/u.s^2);  % slug
        oz = u.lbm/16;              % ounce
        dr = u.oz/16;               % dram
        gr = u.lbm/7000;            % grain
        ct = 200*u.mg;              % carat
        amu = 1.660539040e-27*u.kg; % atomic mass unit
        Da = u.amu;                 % atomic mass unit
        mu = u.amu;                 % atomic mass unit
        mP = 2.176470e-8*u.kg;      % Planck mass
        m_e = 9.10938356e-31*u.kg;  % electron mass
        mug = u.kgf/(u.m/u.s^2);    % metric slug
        hyl = u.mug;                % hyl
        TMU = u.mug;                % technische Masseneinheit
        
        %---- more force ----
        
        pdl = u.lbm*u.ft/u.s^2;     % poundal
        gramForce = u.gram*u.g0;    % gram force
        gf = u.gramForce;           % gram force
        ozf = u.oz*u.g0;            % ounce force
        tonf = u.tn*u.g0;           % short ton force
        
        %---- mass per length ----
        
        den = u.gram/(9*u.km);      % denier
        tex = u.gram/u.km;          % tex
        dtex = u.tex/10;            % decitex
        
        %---- time ----
        
        ms = 1e-3*u.s;              % millisecond
        us = 1e-6*u.s;              % microsecond
        ns = 1e-9*u.s;              % nanosecond
        ps = 1e-12*u.s;             % picosecond
        tP = 5.39116e-44*u.s;       % Planck time
        min = 60*u.s;               % minute
        h = 60*u.min;               % hour
        hr = u.h;                   % hour
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
        
        %---- frequency ----
        
        Hz = 1/u.s;   % hertz (NB: incompatible with angle and angular velocity)
        kHz = 1e3*u.Hz;             % kilohertz
        MHz = 1e6*u.Hz;             % megahertz
        GHz = 1e9*u.Hz;             % gigahertz
        Bd = 1/u.s;                 % baud
        
        %---- energy ----
        
        Nm = u.N*u.m;               % newton meter
        J = u.Nm;                   % joule
        MJ = 1e6*u.J;               % megajoule
        kJ = 1e3*u.J;               % kilojoule
        mJ = 1e-3*u.J;              % millijoule
        uJ = 1e-6*u.J;              % microjoule
        nJ = 1e-9*u.J;              % nanojoule
        eV = 1.6021766208e-19*u.J;  % electronvolt
        BTU = 1055.06*u.J;          % British thermal unit (ISO)
        Btu = u.BTU;                % British thermal unit (ISO)
        Btu_IT = 1055.0559*u.J;     % British thermal unit (International Table)
        Btu_th = 1054.3503*u.J;     % British thermal unit (thermochemical)
        kpm = u.kp*u.m;             % kilopond meter
        Ws = u.J;                   % watt-second
        kWh = 3.6e6*u.J;            % kilowatt-hour
        Wh = 3.6e3*u.J;             % watt-hour
        cal = 4.1868*u.J;           % calorie (International Table)
        cal_IT = u.cal;             % calorie (International Table)
        cal_4 = 4.204*u.J;          % calorie (4°C)
        cal_15 = 4.1855*u.J;        % calorie (15°C)
        cal_20 = 4.182*u.J;         % calorie (20°C)
        cal_mean = 4.190*u.J; 	  % calorie (mean)
        cal_th = 4.184*u.J;         % calorie (thermochemical)
        kcal = 1e3*u.cal;           % kilocalorie
        kcal_IT = 1e3*u.cal_IT;     % kilocalorie (International Table)
        Cal = u.kcal;               % large calorie / food calorie
        kcal_4 = 1e3*u.cal_4;       % kilocalorie (4°C)
        kcal_15 = 1e3*u.cal_15;     % kilocalorie (15°C)
        kcal_20 = 1e3*u.cal_20;     % kilocalorie (20°C)
        kcal_mean = 1e3*u.cal_mean; % kilocalorie (mean)
        kcal_th = 1e3*u.cal_th;     % kilocalorie (thermochemical)
        erg = 1e-7*u.J;             % en.wikipedia.org/wiki/Erg
        E_h = 4.359744650e-18*u.J;  % hartree energy
        Ha = u.E_h;                 % hartree
        thm = 1e5*u.BTU;            % therm
        therm = u.thm;              % therm
        quad = 1e15*u.BTU;          % quad
        
        %---- temperature ----
        % For reference: °C = °K-273.15; °F = °R-459.67.
        
        R = u.K*5/9;                % rankine (°F = °R-459.67)
        mK = 1e-3*u.K;              % millikelvin
        uK = 1e-6*u.K;              % microkelvin
        nK = 1e-9*u.K;              % nanokelvin
        deltaK = u.K;               % kelvin (relative temperature)
        deltadegC = u.K;            % celsius (relative, °C = °K-273.15)
        deltadegR = u.R;            % rankine (relative temperature)
        deltadegF = u.R;            % fahrenheit (relative, °F = °R-459.67)
        TP = 1.416808e32*u.K;       % Planck temperature
        
        %---- pressure ----
        
        Pa = u.N/u.m^2;             % pascal
        mPa = 1e-3*u.Pa;            % millipascal
        kPa = 1e3*u.Pa;             % kilopascal
        MPa = 1e6*u.Pa;             % megapascal
        GPa = 1e9*u.Pa;             % gigapascal
        torr = 133.322*u.Pa;        % torr
        Torr = u.torr;              % torr
        mtorr = 1e-3*u.torr;        % millitorr
        bar = 1e5*u.Pa;             % bar
        mbar = 1e-3*u.bar;          % millibar
        kbar = 1e3*u.bar;           % kilobar
        atm = 101325*u.Pa;          % standard atmosphere
        at = u.kgf/u.cm^2;          % technical atmosphere
        psi = u.lbf/u.in^2;         % pound force per square inch
        ksi = 1e3*u.psi;            % kip per square inch
        Msi = 1e6*u.psi;            % million pound force per square inch
        psf = u.lbf/u.ft^2;         % pound force per square foot
        ksf = u.kip/u.ft^2;         % kip per square foot
        Ba = 0.1*u.Pa;              % barye
        pz = u.kPa;                 % pièze
        mmHg = 13.5951*u.kgf/u.m^2; % millimeter of mercury
        cmHg = 10*u.mmHg;           % centimeter of mercury
        mHg = 1e3*u.mmHg;           % meter of mercury
        inHg = 2.54*u.cmHg;         % inch of mercury
        ftHg = 12*u.inHg;           % foot of mercury
        mmH20 = u.kgf/u.m^2;        % millimeter of water (density 1 g/cc)
        mmAq = u.mmH20;             % millimeter of water
        cmH20 = 10*u.mmH20;         % centimeter of water
        cmAq = u.cmH20;             % centimeter of water
        mH20 = 1e3*u.mmH20;         % meter of water
        mAq = u.mH20;               % meter of water
        inH20 = 2.54*u.cmH20;       % inch of water
        inAq = u.inH20;             % inch of water
        wc = u.inH20;               % inch water column
        ftH20 = 12*u.inH20;         % foot of water
        ftAq = u.ftH20;             % foot of water
        
        %---- viscosity ----
        
        St = u.cm^2/u.s;            % stokes
        cSt = u.St/100;             % centistokes
        newt = u.in^2/u.s;          % newt
        P = u.Pa*u.s / 10;          % poise
        cP = u.mPa*u.s;             % centipoise
        reyn = u.lbf*u.s/u.in^2;    % reyn
        
        %---- power ----
        
        W = u.J/u.s;                % watt
        MW = 1e6*u.W;               % megawatt
        GW = 1e9*u.W;               % gigawatt
        kW = 1e3*u.W;               % kilowatt
        mW = 1e-3*u.W;              % milliwatt
        uW = 1e-6*u.W;              % microwatt
        nW = 1e-9*u.W;              % nanowatt
        pW = 1e-12*u.W;             % picowatt
        hp = 550*u.ft*u.lbf/u.s;    % mechanical horsepower (550 ft-lbf/s)
        HP_I = u.hp;                % mechanical horsepower (550 ft-lbf/s)
        hpE = 746*u.W;              % electrical horsepower
        HP_E = u.hpE;               % electrical horsepower
        PS = 75*u.kg*u.g0*u.m/u.s;  % metric horsepower (DIN 66036)
        HP = u.PS;                  % metric horsepower (DIN 66036)
        HP_DIN = u.PS;              % metric horsepower (DIN 66036)
        
        %---- current ----
        
        mA = 1e-3*u.A;              % milliampere
        uA = 1e-6*u.A;              % microampere
        nA = 1e-9*u.A;              % nanoampere
        pA = 1e-12*u.A;             % picoampere
        kA = 1e3*u.A;               % kiloampere
        abA = 10*u.A;               % abampere
        Bi = u.abA;                 % biot
        
        %---- charge ----
        
        C = u.A*u.s;                % coulomb
        e = 1.6021766208e-19*u.C;   % elementary charge
        mC = 1e-3*u.C;              % millicoulomb
        uC = 1e-6*u.C;              % microcoulomb
        nC = 1e-9*u.C;              % nanocoulomb
        pC = 1e-12*u.C;             % picocoulomb
        abC = 10*u.C;               % abcoulomb
        aC = u.abC;                 % abcoulomb
        statC = u.dyn^(1/2)*u.cm;   % statcoulomb
        Fr = u.statC;               % franklin
        esu = u.statC;              % electrostatic unit of charge
        mAh = u.mA*u.hr;            % milliamp-hour
        Ah = u.A*u.hr;              % amp-hour
        
        %---- voltage ----
        
        V = 1*u.J/u.C;              % volt
        kV = 1e3*u.V;               % kilovolt
        MV = 1e6*u.V;               % megavolt
        GV = 1e9*u.V;               % gigavolt
        mV = 1e-3*u.V;              % millivolt
        uV = 1e-6*u.V;              % microvolt
        
        %---- resistance/conductance ----
        
        Ohm = u.V/u.A;              % ohm
        GOhm = 1e9*u.Ohm;           % gigaohm
        MOhm = 1e6*u.Ohm;           % megaohm
        kOhm = 1e3*u.Ohm;           % kiloohm
        mOhm = 1e-3*u.Ohm;          % milliohm
        uOhm = 1e-6*u.Ohm;          % microohm
        nOhm = 1e-9*u.Ohm;          % nanoohm
        abOhm = u.nOhm;             % abohm
        Z0 = 376.730313461*u.Ohm;   % characteristic impedance of vacuum
        R_K = 25812.8074555*u.Ohm;  % von Klitzing constant
        R_K_90 = 25812.807*u.Ohm;   % von Klitzing constant (conventional value)
        S = 1/u.Ohm;                % siemens
        mS = 1e-3*u.S;              % millisiemens
        uS = 1e-6*u.S;              % microsiemens
        nS = 1e-9*u.S;              % nanosiemens
        G0 = 7.7480917310e-5*u.S;   % conductance quantum
        
        
        %---- capacitance ----
        
        F = u.A*u.s/u.V;            % farad
        mF = 1e-3*u.F;              % millifarad
        uF = 1e-6*u.F;              % microfarad
        nF = 1e-9*u.F;              % nanofarad
        pF = 1e-12*u.F;             % picofarad
        
        %---- inductance ----
        
        H = u.Ohm*u.s;              % henry
        mH = 1e-3*u.H;              % millihenry
        uH = 1e-6*u.H;              % microhenry
        nH = 1e-9*u.H;              % nanohenry
        abH = u.nH;                 % abhenry
        kH = 1e3*u.H;               % kilohenry
        MH = 1e6*u.H;               % megahenry
        GH = 1e9*u.H;               % gigahenry
        
        %---- EM ----
        
        T = 1*u.N/(u.A*u.m);        % tesla
        Gs = 1e-4*u.T;              % gauss
        Wb = u.V*u.s;               % weber
        Mx = 1e-8*u.Wb;             % maxwell
        mWb = u.Wb/1000;            % milliweber
        uWb = 1e-6*u.Wb;            % microweber
        nWb = 1e-9*u.Wb;            % nanoweber
        Oe = 250/pi*u.A/u.m;        % oersted
        Gb = 2.5/pi*u.A;            % gilbert
        
        %---- non-dimensionals ----
        
        percent = 0.01;             % %
        pct = u.percent;            % %
        permil = 0.001;             % ‰
        permill = u.permil;         % ‰
        permille = u.permil;        % ‰
        permyriad = 1e-4;           % permyriad
        bp = u.permyriad;           % basis point
        ppm = 1e-6;                 % part per million
        ppb = 1e-9;                 % part per billion
        ppt = 1e-12;                % part per trillion
        ppq = 1e-15;                % part per quadrillion
        
        %---- angles ----
        % Note: angles are dimensionless
        
        rad = 1;                    % radian
        sr = 1;                     % steradian
        turn = 2*pi*u.rad;          % turn
        rev = u.turn;               % revolution = 2*pi radians
        deg = u.turn/360;           % degree
        arcminute = u.deg/60;       % arcminute
        arcmin = u.arcminute;       % arcminute
        arcsecond = u.arcminute/60; % arcsecond
        arcsec = u.arcsecond;       % arcsecond
        grad = u.turn/400;          % gradian
        
        %---- rotational speed ----
        
        rpm = u.rev/u.min;          % revolution per minute
        rps = u.rev/u.s;            % revolution per second
        
        %---- velocity ----
        
        mps = u.m/u.s;              % meter per second
        kyne = u.cm/u.s;            % kyne
        Kyne = u.kyne;              % kyne
        fps = u.ft/u.s;             % foot per second
        fpm = u.ft/u.min;           % foot per minute
        kt = u.nmi/u.hr;            % knot
        kn = u.kt;                  % knot
        kts = u.kt;                 % knot
        knot = u.kt;                % knot
        knot_UK = u.nm_UK/u.hr;     % British imperial knot
        KTAS = u.kt;                % knot
        nmph = u.kt;                % knot
        kph = u.km/u.hr;            % kilometer per hour
        kmh = u.kph;                % kilometer per hour
        mph = u.mi/u.hr;            % mile per hour
        
        %---- volume flow rate ----
        
        cfm = u.ft^3/u.min;         % cubic foot per minute
        cfs = u.ft^3/u.s;           % cubic foot per second
        gpm = u.gal/u.min;          % US customary gallon per minute
        gpm_UK = u.gal_UK/u.min;    % British imperial gallon per minute
        lpm = u.l/u.min;            % liter per minute
        
        %---- fuel economy ----
        
        l_100km = u.l/(100*u.km);   % liter per 100 km
        mpg = u.mi/u.gal;           % mile per gallon
        
        %---- Luminance etc. ----
        
        asb = u.cd/u.m^2;           % apostilb
        sb = u.cd/u.cm^2;           % stilb
        ph = 1e4*u.cd*u.sr/u.m^2;   % phot
        cp = 0.981*u.cd;            % candlepower
        lm = u.cd*u.sr;             % lumen
        lx = u.lm/u.m^2;            % lux
        nx = 1e-3*u.lx;             % nox
        
        %---- other derived SI ----
        
        kat = u.mol/u.s;            % katal
        M = u.mol/u.L;              % molar
        molarity = u.M;             % molarity
        Nms = u.N*u.m*u.s;          % newton-meter-second
        
        %---- radiation ----
        
        Gy = u.J/u.kg;              % gray
        Sv = u.J/u.kg;              % sievert
        Rad =  u.Gy/100;            % absorbed radiation dose
        rem = u.Sv/100;             % roentgen equivalent man
        roentgen = 2.58e-4*u.C/u.kg;% roentgen
        Ly = u.cal_th/u.cm^2;       % langley
        lan = u.Ly;                 % langley
        Bq = 1/u.s;                 % becquerel
        Ci = 3.7e10*u.Bq;           % curie
        
        %---- constants ----
        
        k_B = 1.38064852e-23*u.J/u.K;       % Boltzmann constant
        sigma_SB = 5.670367e-8*u.W/(u.m^2*u.K^4); % Stefan–Boltzmann constant
        h_c = 6.626070040e-34*u.J*u.s;      % Planck constant
        h_bar = u.h/(2*pi);                 % Dirac constant
        mu_B = 9.274009994e-24*u.J/u.T;     % Bohr magneton
        mu_N = 5.050783699e-27*u.J/u.T;     % nuclear magneton
        c = 299792458*u.m/u.s;              % speed of light in vacuum
        c_0 = u.c;                          % speed of light in vacuum
        ly = u.c*u.year;                    % light-year
        lightYear = u.ly;                   % light-year
        mu0 = pi*4e-7*u.N/u.A^2;            % vacuum permeability
        eps0 = u.c^-2/u.mu0;                % vacuum permittivity
        G = 6.67408e-11*u.m^3/u.kg/u.s^2;   % gravitational constant
        N_A = 6.022140857e23/u.mol;         % Avogadro constant
        NA = u.N_A;                         % Avogadro constant
        NAh = u.N_A*u.h_c;                  % molar Planck constant
        M_u = u.g/u.mol;                    % molar mass constant
        K_J = 483597.8525e9*u.Hz/u.V;       % Josephson constant
        K_J_90 = 483597.9*u.Hz/u.V;         % Josephson constant (conv. value)
        F_c = 96485.33289*u.C/u.mol;        % Faraday constant
        alpha = 7.2973525664e-3;            % fine-structure constant
        c1 = 3.741771790e-16*u.W/u.m^2;     % first radiation constant
        c2 = 1.43877736e-2*u.m*u.K;         % second radiation constant
        b_prime = 5.8789238e10*u.Hz/u.K;    % Wien frequency displ. law const.
        b_c = 2.8977729e-3*u.m*u.K;         % Wien wavelength displ. law const.
        R_air = 287.05287*u.J/u.kg/u.K;     % spec. gas const., air (ESDU 77022)
        R_bar = 8.3144598*u.J/u.mol/u.K;    % molar gas constant
        
        %---- digital information ----
        
        nibble = 4*u.bit;                   % nibble
        B = 8*u.bit;                        % byte
        byte = u.B;                         % byte
        octet = u.B;                        % octet
        kB = 1e3*u.B;                       % kilobyte
        MB = 1e6*u.B;                       % megabyte
        GB = 1e9*u.B;                       % gigabyte
        TB = 1e12*u.B;                      % terabyte
        PB = 1e15*u.B;                      % petabyte
        EB = 1e18*u.B;                      % exabyte
        Kibit = 2^10*u.bit;                 % kibibit
        KiB = 2^10*u.B;                     % kibibyte
        Mibit = 2^20*u.bit;                 % mebibit
        MiB = 2^20*u.B;                     % mebibyte
        Gibit = 2^30*u.bit;                 % gibibit
        GiB = 2^30*u.B;                     % gibibyte
        Tibit = 2^40*u.bit;                 % tebibit
        TiB = 2^40*u.B;                     % tebibyte
        Pibit = 2^50*u.bit;                 % pebibit
        PiB = 2^50*u.B;                     % pebibyte
        Eibit = 2^60*u.bit;                 % exbibit
        EiB = 2^60*u.B;                     % exbibyte
        bps = u.bit/u.s;                    % bit per second
        kbps = 1e3*u.bps;                   % kilobit per second
        Mbps = 1e6*u.bps;                   % megabit per second
        Gbps = 1e9*u.bps;                   % gigabit per second
        Tbps = 1e12*u.bps;                  % terabit per second
        
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
        
    end
    
    %% METHODS
    methods
        %%
        
        %% Plotting and display:
        function disp(o)
            try    dispdisp(o);
            catch; builtin('disp',o);  end
        end
    end
end


% 2015-10-06 Created.
