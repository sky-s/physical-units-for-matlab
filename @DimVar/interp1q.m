function yi = interp1q(x,y,xi)
% See also interp1q.

compatible(x,xi);
outUnits = unitsOf(y);
yi = outUnits*interp1q(double(x),y/outUnits,double(xi));