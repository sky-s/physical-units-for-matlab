% A general test script for physical units plotting.
%   <a href="matlab:runtests('testScript_plotting')">run tests</a>

% Set up by clearing class
close all
clear classes
clear all

displayUnits = {'ft','lb','rpm'};
sqrt(u.acre)
a = 1:50;
b = u.rpm*(5:5:600);
c = u.kg*(1:50);
d = u.Pa*(5:3:20);

%% line
figure
line(a,a*u.ft) 
ylabel('length');
assert(axishasunits(0,'ft'))

%% line2
figure 
line(gca,'XData',a*u.acre,'YData',sqrt(a*u.acre),'LineWidth',7)
xlabel('area')
ylabel('length')
assert(axishasunits('acre','ft'))%FIXME: sqrt(a*u.acre) results in m instead of ft

%% histogram
r = u.R*randp([1,19,22],[1e4,1]);
figure
histogram(r);
xlabel temp
assert(axishasunits('R'))

% historgram (2)
figure
histogram(r,'BinWidth',u.R*.5,'BinLimits',u.R*[15 50])
xlabel temperature
assert(axishasunits('R'));

% histogram (3)
figure
histogram(r/u.m,u.R/u.m*(10:20))
xlabel temp/length
assert(axishasunits(1))

%% histogram2
x = randn(1000,1)*u.kg;
y = randn(1000,1)*u.m;
Xedges = [-Inf -2:0.4:2 Inf]*u.kg;
Yedges = [-Inf -2:0.4:2 Inf]*u.m;
figure
histogram2(x,y,Xedges,Yedges);
% nbins = 5;
% h = histogram2(x,y,nbins);
xlabel mass
ylabel length
assert(axishasunits('kg','m'))

%% plot (1)
% note: variable a has no units.
figure
plot(a,c);
ylabel mass
assert(axishasunits(0,1));

%% plot (2)
figure
plot(c,a);
xlabel mass
assert(axishasunits(1,0));

%% plot (3)
figure
plot(c)
assert(axishasunits(0,1))

%% plot (4)
shoulderror('DimVar:incompatibleUnits',@() plot(a,c,c,a));

%% plot (5)
figure
plot(b',(b').^2,'-ok');
xlabel rpm
ylabel rpm²

%% plot (6)
shoulderror('DimVar:incompatibleUnits','plot',...
    a,a.^2*u.km,'-k',a/2,a,'or','LineStyle',':','Color','m','LineWidth',2);

%% plot (7)
figure
h = plot(b,b.^2,'-k',b/2,b.^2/.4,'or',...
    'LineStyle',':','Color','m','LineWidth',2);
xlabel rpm
ylabel rpm2
assert(all(h(1).XData == h(2).XData*2))

%% plot (8)
figure
plot(b,b.^2,'LineWidth',2);
xlabel rpm
ylabel rpm2
assert(axishasunits('rpm',1))

%% plot (9)
figure
plot(b,b.^2,'LineStyle','--');
xlabel rpm
ylabel rpm2
assert(axishasunits('rpm',1))

%% plot (10)
figure
plot(b,b.^2,'LineWidth',2,'LineStyle',':'); 
xlabel rpm
ylabel rpm2
assert(axishasunits('rpm',1))

%% plot3 (1)
figure
plot3(a,a*u.m,a.^2,'o');
ylabel length
assert(axishasunits(0,'m',0))

%% plot3 (2)
figure
plot3(scd(c,'lb_m'),a*u.m,c.*a,'o');
xlabel mass; ylabel length; zlabel mass
assert(axishasunits('lb_m','m','kg'))

%% plot3 (3)
figure
plot3(b,b.^2,sqrt(b),'-k',...
    b/2,4*b.^2,sqrt(b),'or','LineStyle',':','Color','m')
xlabel rpm
ylabel rpm2
zlabel sqrt(rpm)

%% plot3 (4)
figure
plot3(c,a*u.m,c.*a,'LineStyle',':') 
xlabel mass; ylabel length; zlabel mass
assert(axishasunits(1,1,1))

%% plot3 (5)
figure
plot3(gca,c,a*u.m,c.*a,'LineWidth',2)
xlabel mass; ylabel length; zlabel mass
assert(axishasunits('kg','m','kg'))

%% plot3 (6)
shoulderror('DimVar:incompatibleUnits','plot3',a,a.^2,sqrt(a),'-k',...
    a/2*u.m,4*a.^2,sqrt(a),'or','LineStyle',':','Color','m');

shoulderror('DimVar:incompatibleUnits','plot3',a*u.m,a.^2,sqrt(a),'-k',...
    a/2,4*a.^2,sqrt(a),'or','LineStyle',':','Color','m');

%% plotting with wrong scd
% make sure plot doesn't look at wrong scd and instead gets from displayunits
plot3(scd(c,'hp'),scd(a*u.m,'K'),scd(c.*a,'lb_m'),'LineStyle',':') 
xlabel mass; ylabel length; zlabel mass
assert(axishasunits('kg','m','lb_m'))

%% fill (1)
figure
t = (1/16:1/8:1)'*2*pi;
x = cos(t)*u.cm;
y = sin(t)*u.kg;
fill(x,y,'r')
xlabel length
ylabel mass
assert(axishasunits('cm','kg',0));

%% fill (2)
X = [0 1 1 2; 1 1 2 2; 0 0 1 1];
Y = [1 1 1 1; 1 0 1 0; 0 0 0 0];
Z = [1 1 1 1; 1 0 1 0; 0 0 0 0]*u.kW;
C = [0.5000 1.0000 1.0000 0.5000;
     1.0000 0.5000 0.5000 0.1667;
     0.3330 0.3330 0.5000 0.5000];

figure
fill3(X,Y,Z,C)
zlabel power
assert(axishasunits(0,0,'kW'))

%% patch
x = [1 2 1]; y = [3 4 5]; z = [9 8 2];
figure
patch(x*u.m,y,'b')
xlabel length
assert(axishasunits('m'))

% patch (2)
figure
patch(x,y*u.kgf,'b')
ylabel force
assert(axishasunits(0,'kgf'))

% patch (3)
figure
patch(x*u.kg,y*u.m,scd(z*u.kW/u.kg,'kW/kg'),[.1 .4 .9])
xlabel mass
ylabel length
zlabel powerpermass
view(3)
assert(axishasunits('kg','m','kW/kg'))

% patch (4)
figure
patch('XData',x,'YData',y*u.kW,'ZData',z*u.kgf)
ylabel power
zlabel force
view(3)
assert(axishasunits(0,'kW','kgf'))

% patch (5)
figure
v2 = [2 4; 2 8; 8 4; 5 0; 5 2; 8 0]*u.hp;
f2 = [1 2 3; 
    4 5 6];
patch('Faces',f2,'Vertices',v2,'FaceColor','green')
xlabel power
ylabel power
zlabel power
assert(axishasunits('hp','hp'))

%% expected to fail: patch (6) (struct doesn't call overloaded method)
figure
v2 = [2 4; 2 8; 8 4; 5 0; 5 2; 8 0]*u.hp;
f2 = [1 2 3; 
    4 5 6];
s = struct('Faces',f2,'Vertices',v2,'FaceColor','c');
shoulderror('patch(s)');
xlabel power
ylabel power
zlabel power

%% surf
[x,y]=meshgrid(2:5,12:14);
z = x.^2./y.^2;

% surf (1)
figure
surf(magic(3)*u.m)
zlabel length
assert(axishasunits(0,0,'m'))

% surf (2)
figure
surf(magic(3)*u.m, rand(3))
zlabel length
assert(axishasunits(0,0,'m'))

% surf (3)
figure
surf('ZData',magic(3)*u.kgf)
zlabel force
assert(axishasunits(0,0,'kgf'))

% surf (4)
figure
surf(x*u.m,y,z)
xlabel length
assert(axishasunits('m',0,0))

% surf (5)
figure
surf(x,y*u.kgf,z,rand(size(z)))
ylabel force
assert(axishasunits(0,'kgf',0))

% surf (6)
figure
surf(x*u.kg,y*u.m,z*u.kW/u.kg)
xlabel mass
ylabel length
zlabel powerPerMass
assert(axishasunits('kg','m',1))

% surf (7)
figure
surf('XData',x,'YData',y*u.kW,'ZData',z*u.kgf)
ylabel power
zlabel force
assert(axishasunits(0,'kW','kgf'))

%% mesh, surfc, meshc, surfl, etc.
% TODO

%% text
figure
text(u.ft,u.lb,u.rpm,'text label','Color','red','FontSize',14)
view(3)
xlabel length
ylabel mass
zlabel rpm
assert(axishasunits('ft','lb','rpm'))

%% contourX
x = linspace(-2*pi,2*pi);
y = linspace(0,4*pi);
[X,Y] = meshgrid(x,y);
Z = sin(X)+cos(Y);

% contourc
C = contourc(x*u.ft,y*u.lb,Z*u.rpm);

% contour (1)
figure
contour(X*u.ft,Y*u.lb,Z*u.rpm)
xlabel length
ylabel mass
assert(axishasunits('ft','lb',0))

% contour3
figure
contour3(X*u.ft,Y*u.lb,Z*u.rpm)
xlabel length
ylabel mass
zlabel rpm
view(3)
assert(axishasunits('ft','lb','rpm'))

% contourf
figure
contourf(X*u.ft,Y*u.lb,Z*u.rpm)
xlabel length
ylabel mass
assert(axishasunits('ft','lb',0))

%% bar
figure
subplot(3,1,1), bar(rand(10,5)*u.kWh,'stacked'), colormap(cool)
assert(axishasunits(0,'kWh'))
subplot(3,1,2), bar((0:.25:1)*u.mi,rand(5)*u.lbf,1)
assert(axishasunits('mi','lbf',0))
subplot(3,1,3), bar(rand(2,3)*u.K,.75,'grouped')
assert(axishasunits(0,'K',0))

%% barh
figure
subplot(3,1,1), barh(rand(10,5)*u.kWh,'stacked'), colormap(cool)
assert(axishasunits('kWh'))
subplot(3,1,2), barh((0:.25:1)*u.mi,rand(5)*u.lbf,1)
assert(axishasunits('lbf','mi'))
subplot(3,1,3), barh(rand(2,3)*u.K,.75,'grouped')
assert(axishasunits('K',0))

%% barp3
figure
subplot(1,2,1)
bar3(peaks(5)*u.m)
assert(axishasunits(0,0,'m'))
subplot(1,2,2)
bar3((1:5)*u.kg,rand(5)*u.mi,'stacked')
assert(axishasunits(0,'kg','mi'))

%% bar3h
figure
subplot(1,2,1), bar3h(peaks(5)*u.m,.5)
assert(axishasunits(0,'m',0))
subplot(1,2,2), bar3h((1:5)*u.kg,rand(5)*u.m,'stacked')
assert(axishasunits(0,'m','kg'))

%% multiPlot with different per-variable display units
figure
h1= plot(d,double(d)*u.mm);
h1.LineWidth=8;
hold on
h2 = plot(d,scd(double(d)*u.cm,'mm')/10);
h2.LineWidth=3;
title('2 lines should overlap')
assert(all(abs(h1.XData-h2.XData)<sqrt(eps)))
assert(all(abs(h1.YData-h2.YData)<sqrt(eps)))
assert(axishasunits('Pa','mm'))

%% multiPlot with different per-variable display units, same call
figure
h= plot(d,double(d)*u.mm,d,double(d)*u.cm/10);
h(1).LineWidth=8;
h(2).LineWidth=3;
assert(all(abs(h(1).XData-h(2).XData)<sqrt(eps)))
assert(all(abs(h(1).YData-h(2).YData)<sqrt(eps)))
assert(axishasunits('Pa','mm'))
title('2 lines should overlap')

%% multiPlot with struct style input
[x,y]=meshgrid(2:5,12:14);
z = x.^2./y.^2;

figure
surf('XData',x,'YData',y*u.kW,'ZData',z*u.kgf)
ylabel power
zlabel force


hold on
% patch
x = [1 2 1]; y = [3 4 5]; z = [9 8 2];
patch('XData',x,'YData',y*u.kW,'ZData',z*u.kgf)
ylabel power
zlabel force
view(3)
title('a surf and a patch, about the same scale')
assert(axishasunits(0,'kW','kgf'))

%% _lim functions
figure
linem(u.smoot*rand(3,10))
view(3)
shouldalert('xlim',gca,u.ft*[-6 6]);
ylim(u.smoot*[0 2]);
zlim(u.smoot*[0 2.5])
assert(isequal(zlim,[0 2.5]))

%% helper
function allCorrect = axishasunits(x,y,z)
correct = true(3,1);
ax = gca;
if ~isempty(x)
    correct(1) = checkfmt(ax.XAxis.TickLabelFormat,x);
end
if nargin > 1 && ~isempty(y)
    correct(2) = checkfmt(ax.YAxis.TickLabelFormat,y);
end
if nargin > 2 && ~isempty(z)
    correct(3) = checkfmt(ax.ZAxis.TickLabelFormat,z);
end
allCorrect = all(correct);
end

function tf = checkfmt(s,x)
if ischar(x)
    tf = strcmp(s(4:end),x);
else
    if x
        tf = numel(s) > 3; % '%g xxx'
    else
        tf = strcmp(s,'%g');
    end
end
end

%% external dependencies
function R = randp(P, varargin)
% RANDP pseudorandom integers from a specified discrete distribution
%    R = RANDP(P, N) returns an N-by-N matrix containing pseudorandom
%    integers drawn from a specified discrete distribution on 1:numel(P).
%    The distribution is specified by the relative values of P so that a
%    value K is present approximately "P(K)/sum(P) times in the matrix R. 
%    All values of P should => 0, NaNs are set to 0.
%
%    The other arguments to RANDP specify the size of R in the same way as
%    matlab's own functions do: RANDP(P, N) returns an N-by-N matrix,
%    RANDP(P,M,N) and RANDP(P, [M N]) return M-by-N arrays, etc.
%
%    Examples:
%       % random values from [1 2 4] and a bias for 2
%       R = randp([1 2 0 1], 1, 100) ;  % 100 values
%       histc(R, 1:4)            % -> ~25 ~50 0 ~25
%       
%       % create a random, but biased DNA sequence
%       C ='AGCT', P = [4 1 1 2]
%       DNA = C(randp(P, 1, 50))
%       
%    Also see RAND, RANDPERM
%             RANDPERMBREAK, RANDINTERVAL, RANDSWAP (MatLab File Exchange)
%             RANDSAMPLE (Stats Toolbox)
% Created for Matlab R13+, last tested in 2018a
% version 3.0 (mar 2019)
% (c) Jos van der Geest
% http://www.mathworks.com/matlabcentral/fileexchange/authors/10584
% email: samelinoa@gmail.com
%
% File history:
% 1.0 (nov 2005) - created
% 1.1 (nov 2005) - modified slightly to check input arguments to RAND first
% 1.2 (aug 2006) - fixed bug when called with scalar argument P
% 2.0 (feb 2009) - use HISTC for creating the integers (faster and simplier
%                  than previous algorithm)
% 2.1 (dec 2017) - updated for newer releases
% 2.2 (nov 2018) - updated comments, use cumsum more efficiently
% 3.0 (mar 2019) - updated to avoid cumulative round-off errors in weights,
%                  added check on P, updated code, improved help section
narginchk(2,Inf) ;
P = P(:) ;
if ~isnumeric(P) || isempty(P) 
    error('RANDP:InvalidProbabilityArgument', ...
        'First argument should be a non-empty numerical array.') ;
elseif any(P < 0) 
    error('RANDP:InvalidProbabilitiesNegative', ...
        'All probabilities should be 0 or larger.') ;
elseif any(isinf(P))
   error('RANDP:InvalidProbabilitiesInf', ...
        'Probabilities should be finite.') ;
end
% let rand do all the argument checking for the other arguments
try
    R = rand(varargin{:}) ;    
catch ME
    rethrow(ME) ;
end
P(isnan(P)) = 0 ;
if all(P == 0)
    warning('RANDP:ZeroProbabilities', ...
        'All zero probabilities -> returning zeros') ;
    R(:) = 0 ;
elseif numel(P) == 1
    % a single probability -> all 1
    R(:) = 1 ;
else
    % use histc with numel(P) bins with specific sizes to bin the elements
    P = cumsum(P) ;
    P = [0 ; P] ./ P(end) ;
    P(end) = 1 ; % exact to eliminate cumulative round-off errors
    [~, R] = histc(R, P) ;
end
% Note that RANDP is older than RANDSAMPLE which adopted this technique :-)
end
function h = linem(varargin)
% LINEM  Line using single 2 x n (or 3 x n) matrix as input instead of x, y,
% (and z) vectors. LINEM(M) is equivalent to line(M(1,:),M(2,:)) or
% line(M(1,:),M(2,:),M(3,:)).
% 
%   Scalars and row vectors (including handles, Name/Value pairs, etc.) are
%   preserved and passed to the line function.
% 
%   Example.
%     linem(randn(3,100),'LineStyle','none','Marker','x')
% 
%   See also plot, plot3, line, plotm, transpose.
%   Copyright 2017 Sky Sartorius
%   Contact: www.mathworks.com/matlabcentral/fileexchange/authors/101715
varargin = parsematrixplotting(varargin{:});
h_ = line(varargin{:});
if nargout
    h = h_;
end
end
function varargin = parsematrixplotting(varargin)
% PARSEMATRIXPLOTTING  Expand 2 x n (or 3 x n) matrix in a set of inputs into x,
% y, (and z) vector inputs for plotting functions. 
% 
%   Scalars and row vectors (including handles, Name/Value pairs, etc.) are
%   preserved and passed to the outup.
% 
%   Example.
%     v = parsematrixplotting(randn(3,100),'LineStyle','none','Marker','x');
%     line(v{:})
% 
%   See also plotm, linem, transpose.
%   Copyright 2017 Sky Sartorius
%   Contact: www.mathworks.com/matlabcentral/fileexchange/authors/101715
for i = numel(varargin):-1:1
    M = varargin{i};
    if size(M,1) >= 2
        nDims = size(M,1);
        if nDims > 3
            error('Input must have no more than 3 rows.')
        end
        M = mat2cell(M,ones(nDims,1));
        varargin = [varargin(1:i-1) {M{:}} varargin(i+1:end)]; %#ok<CCAT1>
    end
end
end