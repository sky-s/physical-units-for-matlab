function u = myUnits()
% myUnits  Sets the base set of fundamental units. Copy this to your working
% directory to customize for a particular project.
% 
%   By default, the returned struct must have the fields m, kg, s, A, and K.
% 
%   See also u, units.

%% Basic SI:
% u = units ('SI','-SI');

%% Verbose
% u = units('verbose fps');

%% Fancy SI:

baseUnits = 'SI';

u = units(baseUnits);

dispUnits = {'N'	u.N
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
    '$'     u.USD}; %¤ alternate

u = units(baseUnits,dispUnits);


%% Non-SI; more elaborate display:
% baseUnits = 'fps';
% % dispUnits = '-SI';
% 
% 
% u = units(baseUnits);
% dispUnits = {'RPM/V'	u.rpm/u.V
%     'N-m/A'     u.N*u.m/u.A
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
%     'lbf' u.lbf
%     'g' u.gram
%     'hr' u.hr
%     'psf' u.psf};
% 
% u = units(baseUnits,dispUnits);

%% No units for the sake of speed:
% baseUnits = {'m' 'kg' 's' 'A' 'K' 'mol' 'cd' 'USD'};
% u = cell2struct(num2cell(ones(size(baseUnits))),baseUnits,2);
end

%   Copyright 2015 Sky Sartorius
%   www.mathworks.com/matlabcentral/fileexchange/authors/101715