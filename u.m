classdef u < handle
    % u  Class with DimVars (dimensioned variables) as constant properties.
    %
    %   A dimensioned variable contains both a value (any valid numeric type) and
    %   units (any combination of mass, length, time, etc.). Math operations
    %   performed on dimensioned variables will automatically perform dimensional
    %   analysis to ensure that units are consistent.
    %
    %   To create a variable with physical units information, MULTIPLY by the
    %   appropriate DimVar in u. For example, to define a length of 5 inches: length
    %   = 5 * u.in. You can also use the function str2u.
    %
    %   Math operations on DimVars can create new units. For example:
    %       speed  = 2000 * u.rpm;
    %       torque = 3 * str2u('ft-lbf');
    %       power  = speed * torque; % Returns units of power (e.g. Watts).
    %
    %   If you perform operations in which the units cancel, a normal variable will
    %   be returned. For example, 100 * u.m^2 / u.acre returns a double (and is also
    %   a way to do a simple unit conversion from m² to acres).
    % 
    %   Many MATLAB functions won't accept dimensioned variables (though the list is
    %   growing). See U2NUM.
    % 
    %   The myUnits function allows customization of e.g. preferred display units
    %   for an individual project through the use of the units function.
    %
    %   Calling u by itself will display all available units in u.
    %
    %   IMPORTANT: Both the m-files (u.m, units.m, unitsOf.m, etc.) and the
    %   directory ..\@DimVar must be placed on the current path (or in the
    %   current working directory).
    % 
    %   See also myUnits, units, clear, str2u, u2num,
    %     dispdisp - http://www.mathworks.com/matlabcentral/fileexchange/48637.
    
    %   Units on FEX http://www.mathworks.com/matlabcentral/fileexchange/38977.
    %   Author:     Sky Sartorius
    %               www.mathworks.com/matlabcentral/fileexchange/authors/101715

    properties (Access = 'protected', Constant = true)
        %% Establish base unit system and display preferences:
        U = myUnits()
    end
    properties (Constant = true)
        %% Set base units:
        m = u.U.m
        kg = u.U.kg
        s = u.U.s
        A = u.U.A
        K = u.U.K
        mol = u.U.mol
        cd = u.U.cd
        USD = u.U.USD
        
        %% Other units:
        %----currency----
        cent = u.USD/100;
        
        %----fundamental constants ----
        g0 = 9.80665*u.m/u.s^2; % %exact
        gn = u.g0;
        
        %----units that need to come early----
        lbm = 0.45359237*u.kg; %exact
        
        %------- length ----
        km = 1e3*u.m;
        dm = 1e-1*u.m;
        cm = 1e-2*u.m;
        mm = 1e-3*u.m;
        um = 1e-6*u.m;
        micron = u.um;
        nm = 1e-9*u.m;
        ang = 1e-10*u.m;
        in = 2.54*u.cm; %exact
        mil = 1e-3*u.in;
        ft = 12*u.in;
        kft = 1e3*u.ft;
        yd = 3*u.ft;
        mi = 5280*u.ft;
        nmi = 1852*u.m;
        NM = u.nmi;
        a0 = .529e-10*u.m;
        
        %------- area ------- % 
        ha = 10000*u.m^2;
        hectare = u.ha;
        ac = 43560*u.ft^2;
        acre = u.ac;
        
        %------- volume -------
        cc = (u.cm)^3;
        L = 1000*u.cc;
        mL = u.cc;
        cuin = 16.387064*u.mL; %exact
        gal = 231*u.cuin; %updated (US gallon)
        quart = u.gal/4;
        pint = u.quart/2;
        cup = u.pint/2;
        floz = u.cup/8;
        Tbls = u.floz/2;
        tsp = u.Tbls/3;
        
        %------- acceleration -------
        Gal = u.cm/u.s^2; % en.wikipedia.org/wiki/Gal_(unit)
        
        %---- force -------
        N = u.kg*u.m/u.s^2;
        kN = 1000*u.N;
        dyn = 1e-5*u.N; 
        lbf = 4.4482216152605*u.N; %updated, exact
        kgf = u.kg*u.g0; % kilogram force
        kp = u.kgf; % kilopond
        p = u.kp/1000; % pond
        sn = u.kN; % https://en.wikipedia.org/wiki/Sthene
        pdl = u.lbm*u.ft/u.s^2; % poundal
        
        %----- mass ---------
        gram = 1e-3*u.kg;
        g = u.gram; % Don't confuse u.g with u.g0
        mg = 1e-3*u.gram;
        lb = u.lbm;
        st = 14*u.lbm;
        stone = u.st;
        slug = u.lbf/(u.ft/u.s^2); % updated, exact
        oz = (1/16)*u.lbm;
        amu = 1.66e-27*u.kg;
        t = 1000*u.kg;
        tonne = u.t;
        mug = u.kgf/(u.m/u.s^2);% metric slug
        hyl = u.mug;
        TMU = u.mug; % technische Masseneinheit
        
        %---- more force ------- 2014-11-13
        gramForce = u.gram*u.g0;
        gf = u.gramForce;
        ozf = u.oz*u.g0;
        
        %---- time -------
        ms = 1e-3*u.s;
        us = 1e-6*u.s;
        ns = 1e-9*u.s;
        ps = 1e-12*u.s;
        min = 60*u.s;
        hr = 60*u.min;
        day = 24*u.hr;
        week = 7*u.day;
        fortnight = 2*u.week;
        year = 365.25*u.day; %exact Julian year - defines light-year
        month = u.year/12; %inexact/poorly defined
        
        %---- frequency ----
        Hz = 1/u.s; % INCOMPATIBLE with angle and angular velocity units.
        % Should NOT be u.turn/u.s (open to debate). 2013-10-30
        kHz = 1e3 *u.Hz;
        MHz = 1e6 *u.Hz;
        GHz = 1e9 *u.Hz;
        
        %----- energy -----
        J = u.N*u.m;
        MJ = 1e6*u.J;
        kJ = 1e3*u.J;
        mJ = 1e-3*u.J;
        uJ = 1e-6*u.J;
        nJ = 1e-9*u.J;
        eV = 1.6022e-19*u.J;
        BTU = 1.0550559e3*u.J;
        kWh = 3.6e6*u.J;
        Wh = 3.6e3*u.J;
        cal = 4.1868*u.J;
        kCal = 1e3*u.cal;
        erg = 1e-7*u.J; % https://en.wikipedia.org/wiki/Erg 
        quad = 1e15*u.BTU; %https://en.wikipedia.org/wiki/Quad_(unit)
        
        %---- temperature ---
        R = u.K*5/9;
        % C = K-273.15; just FYI - don't uncomment
        % F = R-459.67; just FYI - don't uncomment
        mK = 1e-3*u.K;
        uK = 1e-6*u.K;
        nK = 1e-9*u.K;
        
        %---- pressure -----
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
        Ba = 0.1*u.Pa;
        pz = u.kPa; % pièze 
        mmHg = 133.322387415*u.Pa;
        inHg = 25.4*u.mmHg; %3386.389*u.Pa 
        
        %---- viscosity ----
        St = u.cm^2/u.s; % Stokes (kinematic viscosity)
        cSt = u.St/100;
        P = u.Pa * u.s / 10; % Poise (dynamic) en.wikipedia.org/wiki/Poise
        cP = u.mPa * u.s;
        
        %----- power --- ---
        W = u.J/u.s;
        MW = 1e6*u.W;
        kW = 1e3*u.W;
        mW = 1e-3*u.W;
        uW = 1e-6*u.W;
        nW = 1e-9*u.W;
        pW = 1e-12*u.W;
        hp = 550*u.ft*u.lbf/u.s; %mechanical horsepower, exact;
        hpE = 746*u.W; %electrical horsepower, exact;
        PS = 75*u.kg*u.g0*u.m/u.s; %DIN 66036 definition of metric horsepower.
        
        %----- Current ------
        % u.A = 1*u.C/u.s;
        mA = 1e-3*u.A;
        uA = 1e-6*u.A;
        nA = 1e-9*u.A;
        
        %------ Charge ------
        C = u.A*u.s;  
        e = 1.6022e-19*u.C;
        mC = 1e-3*u.C;
        uC = 1e-6*u.C;
        nC = 1e-9*u.C;
        pC = 1e-12*u.C;
        
        %------ Voltage -----
        V = 1*u.J/u.C;
        kV = 1e3*u.V;
        mV = 1e-3*u.V;
        uV = 1e-6*u.V;
        
        %----- Resistance/capacitance/inductance ------ 
        Ohm = u.V/u.A; 
        MOhm = 1e6*u.Ohm;
        kOhm = 1e3*u.Ohm;
        mOhm = 1e-3*u.Ohm;
        S = 1/u.Ohm; % Siemens
        
        % Capacitance
        F = u.A*u.s/u.V;
        mF = 1e-3*u.F;
        uF = 1e-6*u.F;
        nF = 1e-9*u.F;
        pF = 1e-12*u.F;
        
        % Inductance
        H = u.Ohm*u.s;
        mH = 1e-3*u.H;
        
        % Capacity 
        mAh = u.mA*u.hr;
        Ah = u.A*u.hr;
        
        %---- EM -----
        T = 1*u.N/(u.A*u.m); %Tesla
        gauss = 1e-4*u.T;
        Wb = u.V*u.s; % Weber 
        mWb = u.Wb/1000;
        uWb = 1e-6*u.Wb;
        nWb = 1e-9*u.Wb;
        
        %----Fundamental constants ----
        kB = 1.38e-23*u.J/u.K;
        sigma_SB = 5.670e-8 * u.W/(u.m^2 * u.K^4);
        h = 6.626e-34 * u.J*u.s; % Planck constant
        hbar = u.h/(2*pi);
        mu_B = 9.274e-24 * u.J/u.T;
        mu_N = 5.0507866e-27 * u.J/u.T;
        c = 299792458*u.m/u.s; %exact speed of light
        eps0 = 8.8541878176204e-12* u.C/(u.V*u.m);
        mu0 = 1.2566370614359e-6 * u.J/(u.m*u.A^2);
        
        % Ideal specific gas constant for dry air; value from ESDU 77022
        Rair = 287.05287*u.J/u.kg/u.K;
        
        %
        ly = u.c*u.year; % 1 light-year
        lightYear = u.ly;
        
        %----non-dimensionals---- 
        percent = 0.01; %
        pct = u.percent;
        permil = 0.001; % ‰
        permill = u.permil;
        permille = u.permil;
        permyriad = 1e-4; % ?
        bp = u.permyriad; % basis point
        ppm = 1e-6; % part per million
        ppb = 1e-9; % part per billion
        ppt = 1e-12; % part per trillion
        ppq = 1e-15; % part per quadrillion % caution: approaching eps
        
        %----angles----
        % Note: angles are dimensionless.
        rad = 1; 
        sr = 1;
        turn = 2*pi*u.rad;
        rev = u.turn;
        deg = u.turn/360;
        arcminute = u.deg/60;
        arcsecond = u.arcminute/60;
        grad = u.turn/400;
        
        %----rotational speeds----
        rpm = u.rev/u.min; 
        
        %----speeds ----
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
        
        %----volume flow rates ---- 
        cfm = u.ft^3/u.min;
        cfs = u.ft^3/u.s;
        
        %----
        kat = u.mol/u.s;
        lm = u.cd*u.sr;
        lx = u.lm/u.m^2;
        
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
%   Copyright 2015 Sky Sartorius
%   www.mathworks.com/matlabcentral/fileexchange/authors/101715