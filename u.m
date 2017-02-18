classdef u < handle
    % u  Physical units. 
    % 
    %   Multiply/divide by u.(unitName) to attach physical units to a variable. For
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
    %   SI base units of kg, m, s, A, K, etc. However, the display of any unit,
    %   whether a base unit or a derived unit, can be customized. For example,
    %   instead of displaying power as [m^2][kg]/s^3, customizing base units could
    %   make it display in terms of [ft^2][lbm]/s^3. Or, if desired for a given
    %   project, all variables with units corresponding to "power" could instead be
    %   displayed in terms of a single derived unit, e.g. W or hp. 
    %   
    %   Display customization is set by whatever myUnits.m file is highest on the
    %   MATLAB search path, so a unique myUnits.m can be placed in a project's
    %   directory to tailor preferred display units for that project. Be sure to
    %   clear the class when changing projects or else the old myUnits.m will remain
    %   in effect.
    % 
    %   Many MATLAB functions won't accept variables with physical units. See u2num.
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
    %   See also myUnits, clear, str2u, u2num,
    %     dispdisp - http://www.mathworks.com/matlabcentral/fileexchange/48637.
    
    %   Copyright Sky Sartorius
    %   www.mathworks.com/matlabcentral/fileexchange/authors/101715
        
    properties (Access = 'protected', Constant = true)
        % Establishes base unit system and display preferences based on whatever
        % myUnits.m file is highest on the MATLAB search path.
        baseUnits = myUnits()
    end
    properties (Constant = true)        
        %% Base units:
        
        m =     u.baseUnits.m       % meter (base unit)
        kg =    u.baseUnits.kg      % kilogram (base unit)
        s =     u.baseUnits.s       % second (base unit)
        A =     u.baseUnits.A       % ampere (base unit)
        K =     u.baseUnits.K       % kelvin (base unit)
        mol =   u.baseUnits.mol     % mole (base unit)
        cd =    u.baseUnits.cd      % candela (base unit)
        USD =   u.baseUnits.USD     % currency (base unit)
        
        %% Derived units:        
        %---- length ----
        
        km = 1e3*u.m;
        dm = 1e-1*u.m;
        cm = 1e-2*u.m;
        mm = 1e-3*u.m;
        um = 1e-6*u.m;
        micron = u.um;
        nm = 1e-9*u.m;
        ang = 1e-10*u.m;            % ångström
        in = 2.54*u.cm; 
        mil = 1e-3*u.in;
        ft = 12*u.in;
        kft = 1e3*u.ft;
        yd = 3*u.ft;
        mi = 5280*u.ft;
        nmi = 1852*u.m;             % nautical mile
        NM = u.nmi;                 % nautical mile
        a0 = .529e-10*u.m;          % Bohr radius
        
        %---- area ---- 
        
        ha = 10000*u.m^2;
        hectare = u.ha;
        ac = 43560*u.ft^2;
        acre = u.ac;
        
        %---- volume ----
        
        cc = (u.cm)^3;
        L = 1000*u.cc;
        mL = u.cc;
        cuin = 16.387064*u.mL;      % cubic inch
        gal = 231*u.cuin;           % US gallon
        quart = u.gal/4;
        pint = u.quart/2;
        cup = u.pint/2;
        floz = u.cup/8;             % fluid ounce
        Tbls = u.floz/2;
        tsp = u.Tbls/3;
        
        %---- acceleration ----
        
        g0 = 9.80665*u.m/u.s^2;     % standard gravity
        gn = u.g0;                  % standard gravity
        Gal = u.cm/u.s^2;           % en.wikipedia.org/wiki/Gal_(unit)
        
        %---- force ----
        
        N = u.kg*u.m/u.s^2;
        kN = 1000*u.N;
        dyn = 1e-5*u.N; 
        lbf = 4.4482216152605*u.N;  % pound force
        kgf = u.kg*u.g0;            % kilogram force
        kp = u.kgf;                 % kilopond
        p = u.kp/1000;              % pond
        sn = u.kN;                  % Sthene
        
        %---- mass ----
        
        gram = 1e-3*u.kg;
        g = u.gram;                 % gram
        mg = 1e-3*u.gram;
        lbm = 0.45359237*u.kg;      % pound mass
        lb = u.lbm;                 % pound mass
        st = 14*u.lbm;              % stone
        stone = u.st;
        slug = u.lbf/(u.ft/u.s^2);
        oz = (1/16)*u.lbm;
        amu = 1.660539040e-27*u.kg; % atomic mass unit
        Da = u.amu;                 % atomic mass unit
        t = 1000*u.kg;
        tonne = u.t;
        mug = u.kgf/(u.m/u.s^2);    % metric slug
        hyl = u.mug;
        TMU = u.mug;                % technische Masseneinheit
        
        %---- more force ----
        
        pdl = u.lbm*u.ft/u.s^2;     % poundal
        gramForce = u.gram*u.g0;
        gf = u.gramForce;           % gram force
        ozf = u.oz*u.g0;            % ounce force
        
        %---- time ----
        
        ms = 1e-3*u.s;
        us = 1e-6*u.s;
        ns = 1e-9*u.s;
        ps = 1e-12*u.s;
        min = 60*u.s;
        hr = 60*u.min;
        day = 24*u.hr;
        week = 7*u.day;
        fortnight = 2*u.week;
        year = 365.25*u.day;        % Julian year; defines light-year
        month = u.year/12;          % 1/12th Julian year
        
        %---- frequency ----
        
        % hertz; note: Incompatible with angle and angular velocity units.
        Hz = 1/u.s;                 
        kHz = 1e3 *u.Hz;
        MHz = 1e6 *u.Hz;
        GHz = 1e9 *u.Hz;
        
        %---- energy ----
        
        J = u.N*u.m;
        MJ = 1e6*u.J;
        kJ = 1e3*u.J;
        mJ = 1e-3*u.J;
        uJ = 1e-6*u.J;
        nJ = 1e-9*u.J;
        eV = 1.6022e-19*u.J;
        BTU = 1.0550559e3*u.J;
        kWh = 3.6e6*u.J;            % kilowatt-hour
        Wh = 3.6e3*u.J;             % watt-hour
        cal = 4.1868*u.J;
        kCal = 1e3*u.cal;
        erg = 1e-7*u.J;             % en.wikipedia.org/wiki/Erg 
        quad = 1e15*u.BTU;          % en.wikipedia.org/wiki/Quad_(unit)
        
        %---- temperature ----  
        % For reference: °C = °K-273.15; °F = °R-459.67.
        
        R = u.K*5/9;
        mK = 1e-3*u.K;
        uK = 1e-6*u.K;
        nK = 1e-9*u.K;
        
        %---- pressure ----
        
        Pa = u.N/u.m^2;
        mPa = u.Pa/1000; 
        kPa = u.Pa*1000;
        MPa = u.kPa*1000;
        torr = 133.322*u.Pa;
        mtorr = 1e-3*u.torr;
        bar = 1e5*u.Pa;
        mbar = 1e-3*u.bar;
        atm = 101325*u.Pa;
        psi = u.lbf/u.in^2; 
        ksi = 1000*u.psi; 
        Msi = 1000*u.ksi;
        psf = u.lbf/u.ft^2;
        Ba = 0.1*u.Pa;              % Barye
        pz = u.kPa;                 % pièze 
        mmHg = 133.322387415*u.Pa;  
        inHg = 25.4*u.mmHg;
        
        %---- viscosity ----
        
        St = u.cm^2/u.s;            % stokes (kinematic viscosity)
        cSt = u.St/100;
        P = u.Pa * u.s / 10;        % poise (dynamic viscosity)
        cP = u.mPa * u.s;
        
        %---- power ----
        
        W = u.J/u.s;
        MW = 1e6*u.W;
        kW = 1e3*u.W;
        mW = 1e-3*u.W;
        uW = 1e-6*u.W;
        nW = 1e-9*u.W;
        pW = 1e-12*u.W;
        hp = 550*u.ft*u.lbf/u.s;    % mechanical horsepower
        hpE = 746*u.W;              % electrical horsepower
        PS = 75*u.kg*u.g0*u.m/u.s;  % metric horsepower (DIN 66036 definition)
        
        %---- current ----

        mA = 1e-3*u.A;
        uA = 1e-6*u.A;
        nA = 1e-9*u.A;
        
        %---- charge ----
        
        C = u.A*u.s;                % coulomb
        e = 1.6022e-19*u.C;         % elementary charge
        mC = 1e-3*u.C;
        uC = 1e-6*u.C;
        nC = 1e-9*u.C;
        pC = 1e-12*u.C;
        
        mAh = u.mA*u.hr;            % milliamp-hour
        Ah = u.A*u.hr;              % amp-hour
        
        %---- voltage ----
        
        V = 1*u.J/u.C;              % volt
        kV = 1e3*u.V;
        mV = 1e-3*u.V;
        uV = 1e-6*u.V;
        
        %---- resistance ----
        
        Ohm = u.V/u.A; 
        MOhm = 1e6*u.Ohm;
        kOhm = 1e3*u.Ohm;
        mOhm = 1e-3*u.Ohm;
        S = 1/u.Ohm;                % siemens
        
        %---- capacitance ----
        
        F = u.A*u.s/u.V;            % farad 
        mF = 1e-3*u.F;
        uF = 1e-6*u.F;
        nF = 1e-9*u.F;
        pF = 1e-12*u.F;
        
        %---- inductance ----
        
        H = u.Ohm*u.s;              % henry 
        mH = 1e-3*u.H;
        
        %---- EM ----
        
        T = 1*u.N/(u.A*u.m);        % tesla
        gauss = 1e-4*u.T;
        Wb = u.V*u.s;               % weber 
        mWb = u.Wb/1000;
        uWb = 1e-6*u.Wb;
        nWb = 1e-9*u.Wb;
        
        %---- fundamental constants ----
        % See http://www.efunda.com/units/show_constants.cfm for more.
        
        kB = 1.38e-23*u.J/u.K;              % Boltzmann constant
        sigma_SB = 5.670e-8 * u.W/(u.m^2 * u.K^4); % Stefan–Boltzmann constant
        h = 6.62607004e-34 * u.J*u.s;       % Planck constant
        hbar = u.h/(2*pi);                  % Dirac constant
        mu_B = 9.27400999e-24 * u.J/u.T;    % Bohr magneton
        mu_N = 5.050783699e-27 * u.J/u.T;   % nuclear magneton
        c = 299792458 * u.m/u.s;            % speed of light in vacuum
        eps0 = 8.8541878176204e-12* u.C/(u.V*u.m); % vacuum permittivity
        mu0 = 1.2566370614359e-6 * u.J/(u.m*u.A^2); % vacuum permeability
        
        % specific gas constant for air (ESDU 77022 definition)
        Rair = 287.05287*u.J/u.kg/u.K;
        
        ly = u.c*u.year;            % light-year
        lightYear = u.ly;
        
        %---- non-dimensionals ---- 
        
        percent = 0.01;             % %
        pct = u.percent;
        permil = 0.001;             % ‰
        permill = u.permil;
        permille = u.permil;
        permyriad = 1e-4;           % ?
        bp = u.permyriad;           % basis point
        ppm = 1e-6;                 % part per million
        ppb = 1e-9;                 % part per billion
        ppt = 1e-12;                % part per trillion
        ppq = 1e-15;                % part per quadrillion
        
        %---- angles ----
        % Note: angles are dimensionless
        
        rad = 1;                    % radian
        sr = 1;                     % steradian 
        turn = 2*pi*u.rad;
        rev = u.turn;               % revolution = 2*pi radians
        deg = u.turn/360;           % degree
        arcminute = u.deg/60;
        arcsecond = u.arcminute/60;
        grad = u.turn/400;          % gradian
        
        %---- rotational speeds ----
        
        rpm = u.rev/u.min; 
        
        %---- speeds ----
        
        mps = u.m/u.s;
        fps = u.ft/u.s;
        kt = u.nmi/u.hr;
        kn = u.kt;
        kts = u.kt;
        knot = u.kt;
        KTAS = u.kt; 
        nmph = u.kt;
        kph = u.km/u.hr;
        mph = u.mi/u.hr;
        fpm = u.ft/u.min;
        
        %---- volume flow rates ---- 
        
        cfm = u.ft^3/u.min;         % cubic feet per minute
        cfs = u.ft^3/u.s;           % cubic feet per second
        
        %---- other derived SI ----
        
        kat = u.mol/u.s;            % katal 
        lm = u.cd*u.sr;             % lumen 
        lx = u.lm/u.m^2;            % lux
        
        %---- currency ----
        % See also mathworks.com/matlabcentral/fileexchange/47255
        
        dollar = u.USD;
        cent = u.USD/100;
        
    end
    
    %% METHODS
    methods
        %% Plotting and display methods:
        function disp(o)
            try    dispdisp(o);
            catch; builtin('disp',o);  end
        end
    end
end


% 2015-10-06 Created.
