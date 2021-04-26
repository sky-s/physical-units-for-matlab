# Physical Units Toolbox
[![View Physical Units Toolbox on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/38977-physical-units-toolbox)

Enables operations using hundreds of supported physical units of measurement and physical constants.

If the Physical Units Toolbox is on your MATLAB path, there is nothing to initialize, import, add to your workspace, or pass to functions - simply multiply/divide by `u.(unitName)` to attach physical units to a variable. 
For example, to define a speed using a supported unit: `carSpeed = 100 * u.kph`. 
Or, define a speed with an unsupported unit as a combination of supported units: `snailSpeed = 20 * u.m/u.week`.

Variables with physical units attached are of the class DimVar ("dimenensioned variable"). Math operations performed on dimensioned variables will automatically perform dimensional analysis and can create new units or cancel units and return a normal variable.

On a global or per-project basis, you can customize to use the base unit system of your choice (e.g. ft-lb-s instead of the SI m-kg-s), as well as customize preferred display units. Display units for any given variable can also be customized. Variables will display in the command window or plot, etc. in terms of those units.

Most common Matlab functions will work with physical units, including many types of plots (with added axis labels).
