% A test script for labeling axes with units.
%   <a href="matlab:runtests('testScript_plotAxesUnits')">run tests</a>

% Set up by clearing class
clear all
close all
a = 1:5;
b = a.^2;

%% overriding with new plots
fig('hold off','-reset')
plot(a,b)

% Replace normal with units
plot(a*u.lb,b);

% Replace units with normal
plot(b,a)

% Replace units with different units
plot(a*u.hp,b*u.kW)

plot(a*u.K,b)


%% polar, polarplot (TODO)
% fig polars -reset
% clf
% t = 0 : .01 : 2*pi;
% polarplot(t, sin(2*t).*cos(2*t)*u.smoot, '--r');

%% Multiplots
fig multiplotOneCall -reset

% Incompatible units.
ME = shoulderror('DimVar:incompatibleUnits','plot',a*u.lb,b,a*u.kgf,sqrt(b));
clf

h = plot(a*u.lb,b,a*u.kg,sqrt(b));
xd = {h.XData};
assert(isequal(h(1).XData,a)) % lb should dominate
assert(isequal(h(2).XData,a*u.kg/u.lb)) % should be converted to lb
assert(~isequal(xd{:})) % Should have different scales.

%% Plotting into existing held axes
fig holding -reset

clf
plot(a,b)
hold on

% Replace normal with units - I wish this would error, and maybe it's possible
% to know that there exists data in the axes and not just a blank canvas
plot(a*u.lb,b);

% Replace units with normal - I wish this would error, but likely impossible
% since it won't call overloaded method
plot(b,a)

% Replace units with different units
shouldalert('plot',a*u.hp,b*u.kW)
shouldalert('plot',a*u.K,b)

%% yyaxis
fig yyaxis -reset
yyaxis right
plot(a*u.lb,b*u.K);
yyaxis left
plot(a*u.lb,b*u.N);
hold on
shouldalert('plot',a*u.kg,b*u.N/exp(1)) % Changing x units
yyaxis right
hold on
shouldalert('plot',a*u.kg,b*u.W) % Changing y units
