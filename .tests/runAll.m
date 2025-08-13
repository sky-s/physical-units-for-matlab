% Run all Physical Units test scripts.
addpath('./.tests/alertChecking/');
addpath('./.tests/fig/');

R{1} = runtests('testScript_subsasgn');
R{2} = runtests('testScript_offsetUnits');
R{3} = runtests('testScript_plotAxesUnits'); close all
R{4} = runtests('testScript_units');
R{5} = runtests('testScript_plotting'); close all
R{6} = runtests('testScript_noBaseUnits');
R{7} = runtests('testScript_SIprefixes');

rmpath('./.tests/alertChecking/');
rmpath('./.tests/fig/');

r = [R{:}];
