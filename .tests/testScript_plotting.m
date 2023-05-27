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

T = table;
T.a = a(:);
T.b = b(1:50)';
T.c = c';
T.d = u.Pa*rand(50,1);

%% line
figure
line(a,a*u.ft) 
ylabel('length');
assert(axishasunits(0,'ft'))

%% line into ax
figure
ax = gca;
line(ax,a,a*u.ft);
assert(axishasunits(0,'ft'))

%% line2
figure 
line(gca,'XData',a*u.acre,'YData',sqrt(a*u.acre),'LineWidth',7)
xlabel('area')
ylabel('length')
assert(axishasunits('acre','ft'))

%% line w/ param= syntax
figure
line(a,a*u.ft,Marker='o',LineWidth=5);
assert(axishasunits(0,'ft'))

%% histogram
r = u.R*(25 + 5*randn([1e4,1]));
figure
histogram(gca,r);
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
assert(axishasunits('lb','ft','lb_m'))

%% Area
fig area
area(magic(5)*u.nmi)
assert(axishasunits(0,'nmi')); clf
area(magic(5)*u.nmi,-6*u.nmi)
assert(axishasunits(0,'nmi')); clf
area(6:10,magic(5)*u.nmi,-6*u.nmi)
assert(axishasunits(0,'nmi')); clf
area((6:10)*u.K,magic(5)*u.nmi,-6*u.nmi)
assert(axishasunits('K','nmi'))
area((6:10)*u.K,magic(5)*u.nmi)
assert(axishasunits('K','nmi'))
area((6:10),magic(5)*u.nmi)
assert(axishasunits(0,'nmi'))
area((6:10)*u.K,magic(5))
assert(axishasunits('K'))


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
assert(axishasunits('hp','hp'))

%% patch (6) TODO
figure
v2 = [2 4; 2 8; 8 4; 5 0; 5 2; 8 0]*u.hp;
f2 = [1 2 3; 
    4 5 6];
s = struct('Faces',f2,'Vertices',v2,'FaceColor','c');

% TODO: ideally this would succeed, but with the struct input, the overloaded
% method isn't called.
% Ideal test:
% patch(s);
% assert(axishasunits('hp','hp'))

% Current test:
shoulderror("patch(s)");
title('Failed patch');

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

%% bar3
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

%% scatter
figure
ax = gca;
scatter(ax,a*u.lb,c,a*u.kph)
assert(axishasunits('lb','kg'))

%% scatter on table
figure
scatter(T,"a","b")

%% swarmchart3
figure
x=u.kg*randi(4,1000,1);
y=u.rpm*randi(4,1000,1);
z=u.J*randn(1000,1);
swarmchart3(x,y,z)
assert(axishasunits('kg','rpm','J'))

%% bubblechart3
figure
bubblechart3(a,c.^2,b(1:50),sqrt(b(1:50)),a)
assert(axishasunits(0,'kg²','rpm'))

%% mesh
figure
ax = gca;
[x,y]=meshgrid(2:5,12:14);
z = x.^2./y.^2;
meshz(ax,z*u.nmi)
assert(axishasunits(0,0,'nmi'))

figure
meshc(x*u.rpm,y,z*u.nmi,sqrt(z))
assert(axishasunits('rpm',0,'nmi'))


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
R = u.smoot*rand(3,10);
line(R(1,:),R(2,:),R(3,:));
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
