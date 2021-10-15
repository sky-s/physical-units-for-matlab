% Run all Physical Units test scripts.
addpath('./.tests/alertChecking/');

R{1} = runtests('testScript_subsasgn');
R{2} = runtests('testScript_offsetUnits');
R{3} = runtests('testScript_plotAxesUnits'); close all
R{4} = runtests('testScript_units');
R{5} = runtests('testScript_plotting'); close all
R{6} = runtests('testScript_noBaseUnits');

rmpath('./.tests/alertChecking/');

%% Examine failed tests.
r = [R{:}];
expectedToFail = contains(string({r.Name}),"expectedToFail");
failuresOfInterest = r([r.Failed] & ~expectedToFail);

if isempty(failuresOfInterest)
    disp('ALL PASSED EXPECTED TO PASS.')
else
    disp(failuresOfInterest)
    error('UNEXPECTED FAILURES.')
end
