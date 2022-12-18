% A test script for overloaded DimVar subsasgn method.
%   <a href="matlab:runtests('testScript_subsasgn')">run tests</a>

% Set up by clearing class
clear all

%% assign dim to dim
v1d = u.m*(1:4);
v1d([2,6]) = 9*u.m; 
assert(isequal(v1d,u.m*[1     9     3     4     0     9]));

%% new variable
clear v
v(3) = u.furlong; 
assert(isequal(v,u.furlong*[0 0 1]))

%% assign normal to dim
v1d = u.m*(1:4); %#ok<NASGU>
shoulderror('DimVar:subsasgn:invalidAssignment','v1d(3) = 8');

%% remove element, using a variable (errors with doubles)
v1d = u.m*(1:4);
a = [];
shoulderror('v1d(3) = a;');

%% remove element, special [] syntax
v1d = u.m*(1:4);
v1d(1) = []; 
assert(isequal(v1d,u.m*[2 3 4]));

%% assign NaN to DimVar allowed case
v1d = u.m*(1:4);
v1d(3) = NaN; 
assert(isequaln(v1d,u.m*[1 2 NaN 4]));

%% assign NaNdim to DimVar
v1d = u.m*(1:4);
v1d(3) = u.m*NaN; 
assert(isequaln(v1d,u.m*[1 2 NaN 4]));

%% assign to empty DimVar
ve = []*u.m;
ve(3) = u.m;
assert(isequal(ve,u.m*[0 0 1]))

%% don't allow shrinking array without consistent units
ve = (1:4)*u.m;
shoulderror('DimVar:incompatibleUnits','ve(3) = []*u.kg;');

%% don't allow changing units of all-nan array
v = nan(3,5)*u.m;
shoulderror('DimVar:incompatibleUnits','v(5) = u.kg');

%% assign to empty, wrong unit
ve = []*u.m;
shoulderror('DimVar:incompatibleUnits','ve(3) = u.kg;');

%% assign dim to empty norm**
ve = [];
ve(3) = u.m;

% test for ideal / aspirational toolbox behavior (reason doesn't work: doesn't
% call overloaded method):
% assert(isequal(ve,u.m*[0 0 1]))

% test that passes in current state of toolbox: 
shoulderror("assert(isequal(ve,u.m*[0 0 1]))");

%% assign dim to empty normal, subsasgn call
ve = [];
S.type = '()';
S.subs = {3};
ve = subsasgn(ve,S,9*u.m);
assert(isequal(ve,u.m*[0 0 9]))

%% assign dim to normal**
v1 = 1:4;
% test for ideal / aspirational toolbox behavior (reason doesn't work: doesn't
% call overloaded method):
% shoulderror("v1([2,6]) = 9*u.lb;")

% test that passes in current state of toolbox: 
v1([2,6]) = 9*u.lb;

%% assign dim to normal, subsasgn call
v1 = 1:4;
S.type = '()';
S.subs = {[2 6]};
shoulderror("subsasgn(v1,S,9*u.m)");

%% assign dim to empty normal
v = [];
v = subsasgn(v,substruct('()',{3}),u.kg);
assert(isequal(v,u.kg*[0 0 1]))

%% assign DimVar to NaN**
v1n = NaN(1,4);
v1n(3) = u.m;
% test for ideal / aspirational toolbox behavior (reason doesn't work: doesn't
% call overloaded method):
% assert(isequal(v1n,u.m*[NaN NaN 1 NaN]))

% test that passes in current state of toolbox: 
shoulderror("assert(isequal(v1n,u.m*[NaN NaN 1 NaN]))")

