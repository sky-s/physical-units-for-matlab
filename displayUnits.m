function dispUnits = displayUnits()
% displayUnits  List of preferred units for displaying quantities with physical
% units. Copy this file to your project directory and modify to customize
% physical units for a particular project.
% 
%   Variables without a preferred display unit (all variables if displayUnits
%   returns an empty cell) will be displayed using a combination of the base
%   units defined in baseUnitSystem.
% 
%   The list returned by displayUnits may be either a cellstr with a list of
%   fields of u (or valid input to str2u) or a 2-column cell array with the
%   first column for display strings and the second containing DimVars with the
%   corresponding unit.
% 
%   Remember to clear classes whenever changing displayUnits.
%   
%   <a href="matlab:edit displayUnits">Open displayUnits.m in editor</a>.
% 
%   See also u, baseUnitSystem, str2u.

%% Standard SI.
dispUnits = {'N' 'Pa' 'J' 'W' 'C' 'V' 'F' 'Ohm' 'S' 'Wb' 'T' 'H'};

%% More elaborate SI.
% dispUnits = {
%     'N'         u.N
%     'Pa'        u.Pa
%     'J'         u.J
%     'W'         u.W
%     'C'         u.C
%     'V'         u.V
%     'F'         u.F
%     char(937)   u.Ohm
%     'S'         u.S
%     'Wb'        u.Wb
%     'T'         u.T
%     'H'         u.H
%     '¤'         u.currency
%     }; 

%% FPS.
% dispUnits = {'lbf' 'psf' 'hp' 'V' 'Ohm'};

%% Use only base units.
% dispUnits = {};