% A general test script for physical units.
%   <a href="matlab:runtests('testScript_units')">run tests</a>

% Set up by clearing class
close all
clear all

% baseUnitSystem = 'none'
displayUnits = {'ft','lb','rpm'};
a = 1:50;
b = u.rpm*(5:5:600);
c = u.kg*(1:50);
d = u.Pa*(5:3:20);

%% u.m example 1
rotationSpeed = 2500 * u.rpm;
torque = 95 * str2u('ft-lbf');  % Use alternate string-based definition.
power = rotationSpeed * torque; % Returns variable with units of power.
horsePower = power / u.hp;      % Convert/cancel units.

%% u.m example 2
100 * u.acre/u.ha;  % Convert 100 acres to hectares.
u.st/u.kg;          % Return conversion factor for stone to kilos.

%% u.m example 3
fieldSize = 3*u.sqkm
% Muliplies and divides remove any per-variable custom display units.
% It's nice to display in the units that make sense for that variable.
rate = 3.7*u.acre/u.day
rate = scd(rate,'sqm/hr')
timeNeeded = fieldSize/rate
timeNeeded = scd(timeNeeded,'month')

%% str2u examples
str2u('kg-m²/s^3')
str2u('-5km/s')
str2u('-5 km / s')

%% unitsOf
V = 8*u.g0;
assert(unitsOf(V) == V/u2num(V));
assert(unitsOf(V) == V/double(V));
assert(unitsOf(9) == 9/double(9));
assert(isequal(1,unitsOf(magic(4))))
assert(isequal(u.ft,unitsOf(magic(4))*u.ft))

%% colon
shoulderror('DimVar:colon:incrementRequired','u.cm:u.km');
g = u.ft:u.in:u.yd;

%% validateattributes
validateattributes(u.m,{'numeric'},{'positive'})
validateattributes(u.m,{'DimVar'},{'positive'})
validateattributes(u.m,{'numeric','DimVar'},{'finite'})

validateattributes(u.m,{'double','DimVar'},{'finite'})

shoulderror("MATLAB:invalidType",...
    @validateattributes,u.m,{'double'},{'nonnegative'});

%% compatible1
shoulderror('DimVar:incompatibleUnits',@()compatible(-5,u.mi));
%% compatible2
compatible(u.Pa,u.psi);
%% iscompatible3a
assert(isequal(false,iscompatible(u.mi,-5)));
%% iscompatible3b
assert(isequal(true,iscompatible(u.mi,u.mm)));
%% compatible4
compatible(u.mi,u.m,u.mm,u.km)
%% compatible5
compatible(u.mi)

%% iscompatible1
assert(isequal(true,iscompatible(u.kg,u.lb,u.mg)))
%% iscompatible2
assert(isequal(true,iscompatible(u.m)))
%% iscompatible3
assert(isequal(true,iscompatible(4,5,6)))
%% iscompatible4
assert(isequal(false,iscompatible(u.kg,u.K,u.mg)))

%% isempty
assert(isempty([]*u.slug));
%% isfinite
assert(isfinite(-u.g0));
%% isnumeric
assert(isnumeric(u.cent));
%% isnan
assert(isnan(nan*u.GB));
%% isinf
assert(isinf(-inf*u.cal));
%% isreal
assert(~isreal(1i*u.hr));

%% istype
assert(istype(u.m,'Length'))
assert(istype(u.m,'DimVar'))
assert(~istype(u.kg,'Length'))
assert(istype(u.kg,'Mass'))
assert(istype(u.day,'Time'))
assert(istype(u.R,'Temperature'))
assert(istype(u.USD,'Currency'))
assert(istype(u.acre,'Area'))
assert(istype(u.gal,'Volume'))
assert(istype(u.g0,'Acceleration'))
assert(istype(u.kN,'Force'))
assert(istype(u.kJ,'Energy'))
assert(istype(u.kPa,'Pressure'))
assert(istype(u.W,'Power'))
assert(istype(u.kph,'Velocity'))

%% duration
shoulderror(@()u2duration(-5));
assert(isduration(u2duration(u.day)))
assert(isequal(duration2u(minutes(5)),5*u.min))
assert(duration2u(hours(3.3))==3.3*u.hour)
assert(isa(u2duration(u.year),'duration'));
assert(isequal(u2duration(u.year_Gregorian),years(1)))

%% meshgrid
A = meshgrid(a,b);
meshgrid(a,b,c);
[B, C] = meshgrid(b,c,d);
assert(ndims(B)==3);
assert(isa(C,'DimVar'));
assert(~isa(A,'DimVar'));

%% ndgrid
[D,C,B,A] = ndgrid(d,c,b,a,5:9);
assert(ndims(B)==5);
assert(isa(C,'DimVar'));
assert(~isa(A,'DimVar'));

%% unitConversionFactor1
assert(isequal(unitconversionfactor('m²',u.cm^2),10000));
%% unitConversionFactor2
assert(isequal(unitconversionfactor(u.ft,u.in),12));
%% unitConversionFactor3
assert(isequal(unitconversionfactor('kW',u.W),1000));
%% unitConversionFactor4
shoulderror(@() unitconversionfactor(u.lbf,u.gram));
%% unitConversionFactor5
unitconversionfactor(u.lbf,u.gram,'Force',true);
%% unitConversionFactor6
assert(isequal(unitconversionfactor(u.m,u.kg,'Force',1),u.m/u.kg));
%% unitConversionFactor7
assert(isequal(unitconversionfactor('m',u.km),u.m/u.km));
%% unitConversionFactor8
assert(isequal(unitconversionfactor(u.kW,'W'),1000));
%% unitConversionFactor9
assert(isequal(unitconversionfactor('ft','in'),12));
%% unitConversionFactor10
assert(isequal(unitconversionfactor(u.radian,u.arcminute),u.rad/u.arcminute));
assert(isequal(unitconversionfactor('MPa/min','(lbf/cm^2)/hr'),...
    u.MPa/u.min/(u.lbf/u.cm^2/u.hr)));
assert(isequal(unitconversionfactor('deg/min',1/u.hour),u.deg/u.min*u.hr));
assert(isequal(unitconversionfactor(u.meter,u.kilogram,'Force',true),u.m/u.kg));
if exist('symunit','file')==2
    assert(isequal(unitconversionfactor(symunit('HP'),u.W),u.HP/u.W));
end

%% str2u1
str2u('kg-m')
%% str2u1.2
str2u('kg-m^2')
%% str2u1.3
str2u('s/kg-m')
%% str2u1.1
assert(u.W==str2u('kg-m²/s^3'));
%% str2u1.1
assert(u.W==str2u('m²-kg/s^3'));
%% str2u1.2
assert(str2u('s^3/m²-kg')==1/u.W);
%% str2u1.3
assert(str2u('m²-s^-3-kg')==u.W);
%% str2u1.4
assert(str2u('per m²-s^-3-kg')==1/u.W);
%% str2u2
assert(u.W==str2u('m^2.0*kg*s^-3'));
%% str2u3
assert(u.rpm==str2u('rev/min'));
assert(isequal(u.rpm,str2u('rpm')));
%% str2u4
assert(str2u('kg*K^(-2/3)-mol')==u.kg*u.K^-(2/3)*u.mol)
%% str2u4.1
assert(str2u('kg-m²/s^-4*K^(1/2)-mol')==u.kg*u.m^2*u.s^4*u.K^.5*u.mol)
%% str2u4.2
assert(str2u('kg-m²/s^-4*K^-(1/2)-mol')==u.kg*u.m^2*u.s^4*u.K^-.5*u.mol)
%% str2u4.4
assert(str2u('kg-m²/s^-4*K^(-1/2)-mol')==u.kg*u.m^2*u.s^4*u.K^-.5*u.mol)
%% str2u4.3
assert(str2u('kg-m²/s^-4*K^.5-mol')==u.kg*u.m^2*u.s^4*u.K^.5*u.mol)
%% str2u4.5
% assert(str2u('km/h-(K/g0)')==u.km/(u.h*(u.K/u.g0))) % not supported
%% str2u5
assert(str2u('ft+in')==u.ft+u.in);
%% str2u6
assert(str2u('hr-a_0-K')==u.hr*u.a_0*u.K);
%% str2u plain (nonDimVar) input
assert(isequal(str2u('5'),str2num('5'))); %#ok<*ST2NM>
assert(isequal(str2u('-4'),str2num('-4')));
%%
assert(isequal(str2u('+0.4+3j'),str2num('+0.4+3j')));
%%
assert(isequal(str2u('pi'),str2num('pi')));
assert(isequal(str2u(''),str2num('')));

%% str2u protect against internal variable 'eval'ing
% exp, rep, (also nargin, e.g.)
shoulderror(@() str2u('exp'));
shoulderror(@() str2u('rep'));
shoulderror(@() str2u('nargin')); % builtin str2num can't do it.
shoulderror(@() str2u('5 rep'));

%% str2u leading numeric or value
assert(str2u('-.5a0/s')==-.5*u.a0/u.s)
assert(str2u('+05 a0/s')==5*u.a0/u.s)
assert(isequal(str2u('-.4 ft'),-.4*u.ft))
assert(isequal(str2u('.4ft'),.4*u.ft))
assert(isequal(str2u('.4 ft'),.4*u.ft))
assert(str2u('-5 1/kt')==-5/u.kt)


% assert(isequal(str2u('pi J'),.4*u.ft))

%% str2u denominator stuff
assert(str2u('1/kt')==1/u.kt)


%% str2u strings
assert(isequal(str2u(["5 ft" "horsepower"]),{5*u.ft u.horsepower}));

%% str2uWithDashDenominator
assert(str2u('lbf/a_0-hr-K')==str2u('lbf/(a_0-hr-K)'));
% assert(str2u('km/h-(K/s)')==str2u('km/(h-(K/s)'));

%% more str2u
assert(str2u(' / hour')==1/u.hr) 
assert(str2u('per-hr')==1/u.hr)
assert(str2u('5/hr')==5/u.hr)
assert(str2u('5 per hour')==5/u.hr)
assert(str2u('+.1 mph')==.1*u.mph)
assert(str2u('5 mile per hour')==5*u.mph)
assert(str2u('per mile')==1/u.mi)
assert(isequal(str2u('1'),1))
assert(isequal(str2u('5 percent'),0.05))
assert(str2u('53 lb/hp-hr-K-mol')==53*u.lb/(u.hp*u.hr*u.K*u.mol))
assert(str2u(' per lb-mi')==1/(u.lb*u.mi))

%% histcounts
R = randn([100,1])*u.nmi;
[n,edge,bin] = histcounts(double(R));
[n2,edge2,bin2] = histcounts(R);
assert(isequal(n,n2))
assert(isequal(bin,bin2))
compatible(R,edge2)
assert(isequal(double(edge2),edge))

%% histcounts2

%% interp1q
assert(isequal(5,interp1q(u.slug*(1:5)',(2:2:10)',2.5*u.slug)))
assert(5*u.m==interp1q(u.slug*(1:5)',u.m*(2:2:10)',2.5*u.slug))
shoulderror('DimVar:incompatibleUnits','interp1q',(1:5)',(2:2:10)',2.5*u.slug);

%% interp1
assert(5==interp1(u.slug*(1:5)',(2:2:10)',2.5*u.slug))
assert(5*u.m==interp1(u.slug*(1:5)',u.m*(2:2:10)',2.5*u.slug))
shoulderror('DimVar:incompatibleUnits','interp1',(1:5)',(2:2:10)',2.5*u.slug);
assert(5*u.m==interp1(u.slug*(1:5)',u.m*(2:2:10)',2.5*u.slug,'pchip'))
assert(5*u.m==interp1(u.slug*(1:5)',u.m*(2:2:10)',2.5*u.slug,'pchip',u.km))
assert(5*u.m==interp1(u.slug*(1:5)',u.m*(2:2:10)',2.5*u.slug,'pchip','extrap'))

%% interp2

[X,Y,V] = peaks(10); [Xq,Yq] = meshgrid(-3:.1:3,-3:.1:3);
Vq = interp2(X,Y,V*u.K,Xq,Yq); 
% surf(Xq,Yq,Vq);
%% interp2 - 2
[X,Y,V] = peaks(10); 
[Xq,Yq] = meshgrid(-3:.1:3,u.ft*(-5:.1:3));
V = u.K*V;
Y = u.ft*Y;
Vq = interp2(X,Y,V,Xq,Yq,'nearest', -20*u.K);
Vq = interp2(V,Xq,Yq,'nearest', -20*u.K);
Xq = Xq + 3;
Yq = Yq + 6*u.ft;
Vq = interp2(V,Xq,Yq);

% surf(Xq,Yq,Vq)

%% interp3
[X,Y,Z,V] = flow(10);
[Xq,Yq,Zq] = meshgrid(.1:.25:10,-3:.25:3,-3:.25:3);
Vq = interp3(X,Y,Z,V*u.K,Xq,Yq,Zq); 
Vq = interp3(u.m*X,Y,Z,V*u.K,u.m*Xq,Yq,Zq,'spline'); 
Vq = interp3(u.m*X,Y,Z,V*u.K,u.m*Xq,Yq,Zq,'spline',8*u.K); 
Vq = interp3(u.m*X,Y*u.kg,Z,V*u.K,u.m*Xq,Yq*u.kg,Zq);
% slice(Xq,Yq,Zq,Vq,[6 9.5],2,[-2 .2]), shading flat

%% interpn
f = @(x,y,z,t) t.*exp(-x.^2 - y.^2 - z.^2);

[x,y,z,t] = ndgrid(-1:0.2:1,-1:0.2:1,-1:0.2:1,0:2:10);
v = f(x,y,z,t);

%     Construct a finer grid:
[xq,yq,zq,tq] = ndgrid(-1:0.05:1,-1:0.08:1,-1:0.05:1,0:0.5:10);

%     Interpolate f on the finer grid by using splines:
vq = interpn(x,y,z,t,v*u.K,xq,yq,zq,tq,'spline');
vq = interpn(x,y,z,t,v*u.K,xq,yq,zq,tq,'spline',nan*u.K);
vq = interpn(x,y,z,t,v*u.K,xq,yq,zq,tq);
vq = interpn(x,u.kg*y,z,t,v,xq,yq*u.kg,zq,tq,'spline');

%     And finally, visualize the function:
% nframes = size(tq, 4);
% for j = 1:nframes
%     slice(yq(:,:,:,j), xq(:,:,:,j), zq(:,:,:,j), vq(:,:,:,j),0,0,0);
%     caxis([0 10]);
%     M(j) = getframe;
% end
% movie(M);

%% norm
assert(isequal(u.m*norm(a),norm(u.m*a)));
assert(isequal(u.W*norm(magic(3)),norm(u.W*magic(3))));

%% mustBe___
mustBeFinite(u.rpm);
mustBeGreaterThan(u.km,u.m);
mustBeGreaterThanOrEqual(u.hp,u.W);
mustBeLessThan(u.lbm,u.kg);
mustBeLessThanOrEqual(0*u.K,273*u.K);
mustBeNegative(-u.mi);
mustBeNonnegative(u.Wb);
mustBeNonpositive(0*u.km);
mustBeNonzero(-u.s);
mustBePositive(u.GB);

shoulderror(@()mustBeFinite(nan*u.rpm))
shoulderror(@()mustBeLessThan(u.lbm,u.lb))

%% string
t =  {u.hp; sqrt(u.K)};
assert(abs(1-str2u(string(u.m*u.hp))/str2u('(m)(hp)'))<1e-5)
assert(abs(1-str2u(string(u.m*u.hp))/str2u('[m][hp]'))<1e-5)
assert(isequal(string(t),[string(u.hp); string(sqrt(u.K))]))
assert(isequal(string(t),string(str2u(string(t)))))

assert(isequal(num2cell(u.A*magic(4)),str2u(string(u.A*magic(4)))))
 
%% displayParser
% clear all
u.Ohm/u.K
u.AvogadroConstant*u.kg/u.kg
u.kg
u.N
u.N*u.m/u.m
