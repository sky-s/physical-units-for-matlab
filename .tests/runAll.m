%% test
feature('locale')
feature('DefaultCharacterSet')
exp = {
    '°'                 % 0 Degrees
    '²'                 % 1 Squared character.
    '³'                 % 2 Cubed character.
    '(^per |^per-|^/)'  % 3 Leading 'per' special case.
    '( per |-per-)'     % 4 Replace per with /
    ')('                % 6 Multiply back-2-back parens.
    ']['                % 7 Multiply back-2-back brackets.
    '-u\.'              % 9 - leading unit is *.
    }

%% Run all Physical Units test scripts.
addpath('./.tests/alertChecking/');
addpath('./.tests/fig/');

feature('DefaultCharacterSet', 'UTF-8');

R{1} = runtests('testScript_subsasgn');
R{2} = runtests('testScript_offsetUnits');
R{3} = runtests('testScript_plotAxesUnits'); close all
R{4} = runtests('testScript_units');
R{5} = runtests('testScript_plotting'); close all
R{6} = runtests('testScript_noBaseUnits');

rmpath('./.tests/alertChecking/');
rmpath('./.tests/fig/');


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
