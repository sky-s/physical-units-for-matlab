function u = myUnits()
% myUnits  Sets the base set of fundamental units and display preferences for
% variables with physical units. Copy this file to your project directory and
% modify to customize display of physical units for a particular project.
%   
%   <a href="matlab:edit myUnits">Open myUnits.m in editor</a>.
% 
%   See also u, units.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The returned struct must have the fields m, kg, s, A, K, mol, cd, and USD.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Basic SI, no derived units.
% u = units;

%% Basic FPS, no derived units.
% u = units('FPS');

%% SI with common derived units.
% u = units('SI','-SI');

%% SI with common derived units and using special symbols.
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
    char(937)   u.Ohm
    'S'     u.S
    'Wb'    u.Wb
    'T'     u.T
    'H'     u.H
    '$'     u.USD % Alternative for more generic currency: ¤
    }; 

u = units(baseUnits,dispUnits);


%% Highly customized display with lots of custom display derived units:
% 
% baseUnits = 'FPS'; 
% 
% u = units(baseUnits);
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

%% Do not use phsical units at all and just define all base units as 1. 
% % (Use when calls to DimVar methods are significanly slowing execution.)
% [u.m, u.kg, u.s, u.A, u.K, u.mol, u.cd, u.USD] = deal(1);
end

