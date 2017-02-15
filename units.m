function varargout = units(unitSystem,varargin)
% UNITS returns a struct. Each field of the struct contains a variable of
% the DimVar (Dimensioned Variable) class. 
% 
%   HIGHLY RECOMMENDED: Use the classdef file u.m (which contains constant
%   properties for all the units) and myUnits.m. In this case it is not
%   necessary to set the struct u using the UNITS function, nor is it 
%   necessary to pass the units struct between functions.
% 
%   u = UNITS creates the UNITS struct, u.
% 
%   UNITS with no output assigns the units struct to the variable 'u' in
%   the caller workspace. No warning is given if u already exists.
% 
%   UNITS baseUnitStr or UNITS(baseUnitStr) creates a units struct based on
%   the base units as indicated by baseUnitStr. Some common possibilities
%   are:
%     MKS       - Standard SI/metric (default).
%     CGS       - cm-g-s.
%     MMGS      - mm-g-s.
%     FPS/EE/AE - ft-lb-s / English Engineering / Absolute English.
%     FSS/BG    - ft-slug-s / Gravitational FPS / Technical FPS / 
%                 British Gravitational.
%     Verbose   - Spelled-out base units (meter, kilogram, second, etc.).
% 
%   UNITS(baseUnitSystem) creates a units struct based on the base units
%   provided in the n-by-2 cell array baseUnitSystem. The first element in
%   each row of the cell array is the string for the base unit, and the
%   second element is the number that the unit must be multiplied by to get
%   to the corresponding SI base unit. The rows must be in order of length,
%   mass, time, current, absolute temperature. Not all rows must be
%   provided if only, e.g., non-SI length and mass are desired. Any number
%   of additional base units may be provided (with a multiplier of 1), e.g.
%   'USD' or 'mol'.
% 
%   UNITS(..., dispUnits) allows for displaying derived units (or base
%   units) using preferred display strings instead of combinations of base
%   units, e.g. displaying in watts or horsepower instead of
%   [m^2][kg]/s^3. See documentation for details on use of dispUnits.
% 
%   See also u, myUnits, u2num, unitsOf, DimVar.unitslist, DimVar.strconv, 
%     DimVar.displayunitsstruct, duration, u2duration, 
%     DISPDISP - http://www.mathworks.com/matlabcentral/fileexchange/48637.

% Original seed for this class by Rob deCarvalho - many more units are now
% supported, definitions have been updated to reflect official exact
% conversions, and many useful features and methods have been added and
% modified.
%       http://www.mathworks.com/matlabcentral/fileexchange/authors/22148
%       http://www.mathworks.com/matlabcentral/fileexchange/10070


%% ----- Set up the fundamental dimensions over which to assign units -----
%  Note: Add any dimension you want (e.g. 'dollar') to any list, but never
%  delete or reorder the fundamental base units.


defaultUnitSystem = 'SI'; % Put your favorite base unit system here.

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

%% Shorthand unit systems
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

%%

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
    % 2013-07-15/Sartorius: got rid of eval scheme
end

%------ Get back into SI units for defining all subsequent units.------
for i = 1:nSIunits
    u.(siUnitSystem{i,1}) = u.(dimensionNames{i})*multipliers(i);
end 

%% -------- Define useful units  --------
%----currency----
u.cent = u.USD/100;
        
%----fundamental constants ----
u.g0 = 9.80665*u.m/u.s^2; % %exact
u.gn = u.g0;

%----units that need to come early----
u.lbm = 0.45359237*u.kg; %exact

%------- length ----
u.km = 1e3*u.m; 
u.dm = 1e-1*u.m;
u.cm = 1e-2*u.m;
u.mm = 1e-3*u.m;
u.um = 1e-6*u.m;
u.micron = u.um;
u.nm = 1e-9*u.m;
u.ang = 1e-10*u.m;
u.in = 2.54*u.cm; %exact
u.mil = 1e-3*u.in;
u.ft = 12*u.in;
u.kft = 1e3*u.ft;
u.yd = 3*u.ft;
u.mi = 5280*u.ft;
u.nmi = 1852*u.m;
u.NM = u.nmi; % added 2013-04-16
u.a0 = .529e-10*u.m;

%------- area ------- % added 2014-08-18
u.ha = 10000*u.m^2; 
u.hectare = u.ha;
u.ac = 43560*u.ft^2;
u.acre = u.ac;

%------- volume -------
u.cc = (u.cm)^3;
u.L = 1000*u.cc;
u.mL = u.cc;
u.cuin = 16.387064*u.mL; %exact
u.gal = 231*u.cuin; %updated (US gallon)
u.quart = u.gal/4;
u.pint = u.quart/2;
u.cup = u.pint/2;
u.floz = u.cup/8;
u.Tbls = u.floz/2;
u.tsp = u.Tbls/3;

%------- acceleration -------
u.Gal = u.cm/u.s^2; % added 2014-02-16 en.wikipedia.org/wiki/Gal_(unit)

%---- force -------
u.N = u.kg*u.m/u.s^2;
u.kN = 1000*u.N;
u.dyn = 1e-5*u.N; % changed from full-spelled dyne 2014-02-16
u.lbf = 4.4482216152605*u.N; %updated, exact
u.kgf = u.kg*u.g0; % kilogram force
u.kp = u.kgf; % kilopond, added 2013-09-05
u.p = u.kp/1000; % pond, added 2014-02-16
u.sn = u.kN; % added 2014-02-16 https://en.wikipedia.org/wiki/Sthene
u.pdl = u.lbm*u.ft/u.s^2; % poundal added 2014-02-16

%----- mass ---------
u.gram = 1e-3*u.kg;
u.g = u.gram; % Don't confuse u.g with u.g0
u.mg = 1e-3*u.gram;
u.lb = u.lbm;
u.st = 14*u.lbm;
u.stone = u.st;
u.slug = u.lbf/(u.ft/u.s^2); % updated, exact
u.oz = (1/16)*u.lbm;
u.amu = 1.66e-27*u.kg;
u.t = 1000*u.kg;
u.tonne = u.t;
u.mug = u.kgf/(u.m/u.s^2);% metric slug  added 2014-02-16
u.hyl = u.mug;
u.TMU = u.mug; % technische Masseneinheit

%---- more force ------- 2014-11-13
u.gramForce = u.gram*u.g0;
u.gf = u.gramForce;
u.ozf = u.oz*u.g0;

%---- time -------
u.ms = 1e-3*u.s;
u.us = 1e-6*u.s;
u.ns = 1e-9*u.s;
u.ps = 1e-12*u.s;
u.min = 60*u.s;
u.hr = 60*u.min;
u.day = 24*u.hr;
u.week = 7*u.day;
u.fortnight = 2*u.week;
u.year = 365.25*u.day; %exact Julian year - defines light-year
u.month = u.year/12; %inexact/poorly defined

%---- frequency ---- 
u.Hz = 1/u.s; % INCOMPATIBLE with angle and angular velocity units.
% Should NOT be u.turn/u.s (open to debate). 2013-10-30
u.kHz = 1e3 *u.Hz;
u.MHz = 1e6 *u.Hz;
u.GHz = 1e9 *u.Hz;

%----- energy -----
u.J = u.N*u.m;
u.MJ = 1e6*u.J;
u.kJ = 1e3*u.J;
u.mJ = 1e-3*u.J;
u.uJ = 1e-6*u.J;
u.nJ = 1e-9*u.J;
u.eV = 1.6022e-19*u.J;
u.BTU = 1.0550559e3*u.J;
u.kWh = 3.6e6*u.J;
u.Wh = 3.6e3*u.J;
u.cal = 4.1868*u.J;
u.kCal = 1e3*u.cal;
u.erg = 1e-7*u.J; % https://en.wikipedia.org/wiki/Erg added 2014-02-16
u.quad = 1e15*u.BTU; %https://en.wikipedia.org/wiki/Quad_(unit)

%---- temperature ---
u.R = u.K*5/9;
% C = K-273.15; just FYI - don't uncomment
% F = R-459.67; just FYI - don't uncomment
u.mK = 1e-3*u.K;
u.uK = 1e-6*u.K;
u.nK = 1e-9*u.K;

%---- pressure -----
u.Pa = u.N/u.m^2;
u.mPa = u.Pa/1000; % added 2013-09-29
u.kPa = u.Pa*1000;
u.MPa = u.kPa*1000;
u.torr = 133.322*u.Pa;
u.mtorr = 1e-3*u.torr;
u.bar = 1e5*u.Pa;
u.mbar = 1e-3*u.bar;
u.atm = 101325*u.Pa;
u.psi = u.lbf/u.in^2; %updated, exact
u.ksi = 1000*u.psi; % added 2013-10-29/Sartorius
u.Msi = 1000*u.ksi; % added 2013-10-29/Sartorius
u.psf = u.lbf/u.ft^2;
u.Ba = 0.1*u.Pa;
u.pz = u.kPa; % pièze added 2014-02-16
u.mmHg = 133.322387415*u.Pa; % added 2014-02-16
u.inHg = 25.4*u.mmHg; %3386.389*u.Pa added 2014-02-16

%---- viscosity ----
% added 2013-09-29/Sartorius
u.St = u.cm^2/u.s; % Stokes (kinematic viscosity)
u.cSt = u.St/100;
u.P = u.Pa * u.s / 10; % Poise (dynamic) en.wikipedia.org/wiki/Poise
u.cP = u.mPa * u.s;

%----- power --- ---
u.W = u.J/u.s;
u.MW = 1e6*u.W;
u.kW = 1e3*u.W;
u.mW = 1e-3*u.W;
u.uW = 1e-6*u.W;
u.nW = 1e-9*u.W;
u.pW = 1e-12*u.W;
u.hp = 550*u.ft*u.lbf/u.s; %mechanical horsepower, exact;
u.hpE = 746*u.W; %electrical horsepower, exact;
u.PS = 75*u.kg*u.g0*u.m/u.s; %DIN 66036 definition of metric horsepower.

%----- Current ------
% u.A = 1*u.C/u.s;
u.mA = 1e-3*u.A;
u.uA = 1e-6*u.A;
u.nA = 1e-9*u.A;

%------ Charge ------
u.C = u.A*u.s;  % Switched from coul to C 2014-02-09
u.e = 1.6022e-19*u.C;
u.mC = 1e-3*u.C; %added 2014-02-09
u.uC = 1e-6*u.C;
u.nC = 1e-9*u.C;
u.pC = 1e-12*u.C;

%------ Voltage -----
u.V = 1*u.J/u.C;
u.kV = 1e3*u.V;
u.mV = 1e-3*u.V;
u.uV = 1e-6*u.V;

%----- Resistance/capacitance/inductance ------ added 2014-02-09
u.Ohm = u.V/u.A; %?
u.MOhm = 1e6*u.Ohm;
u.kOhm = 1e3*u.Ohm;
u.mOhm = 1e-3*u.Ohm;
u.S = 1/u.Ohm; % Siemens

% Capacitance
u.F = u.A*u.s/u.V;
u.mF = 1e-3*u.F;
u.uF = 1e-6*u.F;
u.nF = 1e-9*u.F;
u.pF = 1e-12*u.F;

% Inductance
u.H = u.Ohm*u.s;
u.mH = 1e-3*u.H;

% Capacity % added 2014-09-11
u.mAh = u.mA*u.hr;
u.Ah = u.A*u.hr;

%---- EM -----
u.T = 1*u.N/(u.A*u.m); %Tesla
u.gauss = 1e-4*u.T;
u.Wb = u.V*u.s; % Weber added 2014-02-16
u.mWb = u.Wb/1000;
u.uWb = 1e-6*u.Wb;
u.nWb = 1e-9*u.Wb;

%----Fundamental constants ----
u.kB = 1.38e-23*u.J/u.K;
u.sigma_SB = 5.670e-8 * u.W/(u.m^2 * u.K^4);
u.h = 6.626e-34 * u.J*u.s; % Planck constant
u.hbar = u.h/(2*pi);
u.mu_B = 9.274e-24 * u.J/u.T;
u.mu_N = 5.0507866e-27 * u.J/u.T;
u.c = 299792458*u.m/u.s; %exact speed of light
u.eps0 = 8.8541878176204e-12* u.C/(u.V*u.m);
u.mu0 = 1.2566370614359e-6 * u.J/(u.m*u.A^2);

% Ideal specific gas constant for dry air; value from ESDU 77022
u.Rair = 287.05287*u.J/u.kg/u.K; 

%
u.ly = u.c*u.year; % 1 light-year; added 2015-02-16
u.lightYear = u.ly;

%----non-dimensionals---- 2014-11-19
u.percent = 0.01; %
u.pct = u.percent;
u.permil = 0.001; % ‰
u.permill = u.permil;
u.permille = u.permil;
u.permyriad = 1e-4; % ?
u.bp = u.permyriad; % basis point
u.ppm = 1e-6; % part per million
u.ppb = 1e-9; % part per billion
u.ppt = 1e-12; % part per trillion
u.ppq = 1e-15; % part per quadrillion % caution: approaching eps

%----angles---- 
% Note: angles are dimensionless.
% Updated 2013-10-30/Sartorius
u.rad = 1;
u.sr = 1;
u.turn = 2*pi*u.rad;
u.rev = u.turn;
u.deg = u.turn/360;
u.arcminute = u.deg/60;
u.arcsecond = u.arcminute/60;
u.grad = u.turn/400;

%----rotational speeds----
u.rpm = u.rev/u.min; % Added 2013-10-30

%----speeds ----
u.mps = u.m/u.s; % added 2016-02-15
u.fps = u.ft/u.s;
u.kt = u.nmi/u.hr;
u.kn = u.kt;
u.kts = u.kt;
u.knot = u.kt;
u.KTAS = u.kt; % added 2013-09-10
u.nmph = u.kt;
u.kph = u.km/u.hr;
u.mph = u.mi/u.hr;
u.fpm = u.ft/u.min;

%----volume flow rates ---- % added 2014-06-12
u.cfm = u.ft^3/u.min;
u.cfs = u.ft^3/u.s;

%----
u.kat = u.mol/u.s;
u.lm = u.cd*u.sr;
u.lx = u.lm/u.m^2;

%%
if nargout
    varargout = {u};
else
    % Use assignin to have e.g. "units fps" syntax (no "u =...").
    assignin('caller','u',u);
%     disp('Units struct assigned to variable ''u''.');
end
    
end


    

% Revision history/Sartorius
%{
added useful units, especially ones for for aerospace: slug, psf, lbm, 
    fps, kt, deg, etc.
    uploaded to file exchange
 
    added tonne, t, fpm; updated lbf definition and therefore its
    dependents (slug, etc.); changed g0; made year exact julian year
    (defines light-year); added horsepower units; lots of other stuff
    addded kilogram force
changed documentation, refined examples, u2num, unitOf, and added meshgrid
function
2013-04-12
    all new subsasgn method file
    reworked whole help block and published showdemo html 
    small tweaks to u2num and unitsOf help blocks
    uploaded new version to FEX
2013-04-22
    added nnz method
2013-07-15/Sartorius
    Changed DimensionedVariable.m into a more familiar-form classdef file
    highly modified display method with three display alternatives
    Uploaded to FEX
2013-07-17/Sartorius - removed unwanted method help block behavior
2013-07-19/Sartorius: fixed issue with no-numerator display; uploaded
2013-08-28/Sartorius Static method unitslist created; uploaded to FEX
2013-09-02/Sartorius Added u.kp; improved unitslist method
2013-09-09 uploaded to FEX
2013-10-29 added ksi and Msi; added static method strconv
2013-10-30 added angles and rpm; updated strconv; uploaded to FEX
2014-02-10 added electrical units; added option input to select base system
2014-03-09 packaged up with lots of new stuff; uploaded to FEX
2014-05-14 updated several methods to get fuller functionality
2014-05-14 simplified several methods, mostly eliminating extra isa(...)
2014-05-15 added compatible method, greatly simplifying a lot of stuff
calls
2014-05-16 added clearcanceledunits method, also simplifying.
2014-05-16 made subsref method way simpler and faster (no more eval).
2014-08-25 began changing name from Dimensioned Variable to DimVar
2014-11-19 no more auto display of available units (static method replaces
this functinality. Help block and showdemo updated.
2015-06-04 fixed error with Rankine. Uploaded to FEX.
2015-08-20 dispUnits functionality input checking; documentation
2016-01-14 upload to FEX.
2016-01-27 added mm-g-s.
2016-02-15 added mps for m/s.
%}
