function tf = istype(v,name)
% ISTYPE  Determine if input DimVar is of a specified category of units. 
% 
%   ISTYPE(dv,'Unit Category') returns true if dv is a DimVar of the category.
%   Categories are case sensitive, and valid options are: Length, Mass, Time,
%   Temperature, Currency, Area, Volume, Acceleration, Force, Energy, Pressure,
%   Power, Velocity. Otherwise, built-in isa is used.
% 
%   See also isa, u.

tf = false;

switch name
    case 'Length'
        if nnz(v.exponents) == 1 && v.exponents(1) == 1
            tf = true;
        end
    case 'Mass'
        if nnz(v.exponents) == 1 && v.exponents(2) == 1
            tf = true;
        end
    case 'Time'
        if nnz(v.exponents) == 1 && v.exponents(3) == 1
            tf = true;
        end
    case 'Temperature'
        if nnz(v.exponents) == 1 && v.exponents(5) == 1
            tf = true;
        end
    case 'Currency'
        if nnz(v.exponents) == 1 && v.exponents(9) == 1
            tf = true;
        end
    case 'Area'
        if nnz(v.exponents) == 1 && v.exponents(1) == 2
            tf = true;
        end
    case 'Volume'
        if nnz(v.exponents) == 1 && v.exponents(1) == 3
            tf = true;
        end
    case 'Acceleration'
        if nnz(v.exponents) == 2 && isequal(v.exponents(1:3),[1 0 -2])
            tf = true;
        end
    case 'Force'
        if nnz(v.exponents) == 3 && isequal(v.exponents(1:3),[1 1 -2])
            tf = true;
        end
    case 'Energy'
        if nnz(v.exponents) == 3 && isequal(v.exponents(1:3),[2 1 -2])
            tf = true;
        end
    case 'Pressure'
        if nnz(v.exponents) == 3 && isequal(v.exponents(1:3),[-1 1 -2])
            tf = true;
        end
    case 'Power'
        if nnz(v.exponents) == 3 && isequal(v.exponents(1:3),[2 1 -3])
            tf = true;
        end
    case 'Velocity'
        if nnz(v.exponents) == 2 && isequal(v.exponents(1:3),[1 0 -1])
            tf = true;
        end
    otherwise
        tf = isa(v,name);
end