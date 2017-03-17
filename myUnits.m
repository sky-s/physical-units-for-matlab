function u = myUnits()
% myUnits  Sets the base set of fundamental units and display preferences for
% variables with physical units. Copy this file to your project directory and
% modify to customize display of physical units for a particular project.
% 
%   Add as many fields as desired (DimVar or otherwise) to the output struct of
%   this function if project-specific units are desired, accessible via
%   u.U.(field).
%   
%   <a href="matlab:edit myUnits">Open myUnits.m in editor</a>.
% 
%   See also u, units.
% 
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Returned struct must have fields: m, kg, s, A, K, mol, cd, and currency.
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Basic SI, no special display for derived units.
% u = units;

%% Basic FPS, no special display for derived units.
% u = units('FPS');

%% SI with display of common derived units.
% u = units('SI','-SI');

%% SI with display of common derived units, using special symbols.
% Requires 2 calls to units.m. Important: use same base units for both calls!

baseUnits = 'SI';

u = units(baseUnits); % Use for following definition of display units.

dispUnits = {
    'N'     u.N
    'Pa'    u.Pa
    'J'     u.J
    'W'     u.W
    'C'     u.C
    'V'     u.V
    'F'     u.F
    char(937) u.Ohm
    'S'     u.S
    'Wb'    u.Wb
    'T'     u.T
    'H'     u.H
    '¤'     u.currency % Symbol for generic currency.
    }; 

u = units(baseUnits,dispUnits);


%% Highly customized display with lots of custom display derived units:
% Requires 2 calls to units.m. Important: use same base units for both calls!
% 
% baseUnits = 'FPS'; 
% 
% u = units(baseUnits); % Use for following definition of display units.
% 
% dispUnits = {
%     'RPM/V'	u.rpm/u.V
%     'N-m/A' u.N*u.m/u.A
%     'W-hr'	u.W*u.hr
%     'RPM'   u.rpm
%     'W'     u.W
%     'V'     u.V
%     ['N-m/' char(8730) 'W'] u.N*u.m/sqrt(u.W)
%     'W-hr/kg'	u.Wh/u.kg
%     'W/kg'  u.W/u.kg
%     'mAh'	u.mAh
%     char(937)	u.Ohm
%     '°K'	u.K
%     'µF'	u.uF
%     'lbf'   u.lbf
%     'hr'    u.hr
%     'psf'   u.psf
%     };
% 
% u = units(baseUnits,dispUnits);

%% Bespoke base units, including using preferred, non-generic currency:
% Requires 2 calls to units.m. Important: use same base units for both calls!
% 
% baseUnits = {
%     'cm'    100     % 100   cm  / m
%     'g'     1000    % 1000  g   / kg
%     's'     1
%     'A'     1
%     'K'     1
%     'mol'   1
%     'cd'    1
%     'bit'   1
%     'USD'   1       % 1     USD / generic currency
%     };
% 
% u = units(baseUnits); % Use for following definition of display units.
% 
% dispUnits = {
%     'N'     u.N
%     'Pa'    u.Pa
%     'J'     u.J
%     'W'     u.W
%     'C'     u.C
%     'V'     u.V
%     'F'     u.F
%     char(937)   u.Ohm
%     'S'     u.S
%     'Wb'    u.Wb
%     'T'     u.T
%     'H'     u.H
%     '$'     u.USD
%     };
% 
% u = units(baseUnits,dispUnits);

%% Do not use phsical units at all and just define all base units as 1. 
% (Use when calls to DimVar methods are significanly slowing execution.)

% [u.m, u.kg, u.s, u.A, u.K, u.mol, u.cd, u.bit, u.currency] = deal(1);

% While writing code using physical units makes code faster to develop, easier
% to read, and easier to debug, it may slow the code down due to all the calls
% to the DimVar methods. One good solution to working with the class is to
% develop the code using the class, but _after_ it's all tested and debugged,
% use this line in myUnits.m to do all math with normal variables.
end

