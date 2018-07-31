function str = mat2str(dv,varargin)
% str = mat2str(dv)  Appends units to normal mat2str such that eval(str) returns
% a DimVar. Additional arguments are passed to built-in mat2str.
% 
%   See also mat2str.

%   Sky Sartorius 
%   www.mathworks.com/matlabcentral/fileexchange/authors/101715


[dispVal,~,unitStr] = displayparser(dv);


% Interpret everything prior to the first alphabetic character (incl. case of
% leading - or .) as the value.

exp = {'^[-+.0-9]+' ')('  ']['  '([A-Za-z]+\w*)' '-(?=[A-Za-z]+)'  '²'  '³' };
rep = {'$0*'        ')*(' ']*[' 'u.$0'           '*'               '^2' '^3'};

str = [mat2str(dispVal,varargin{:}) ' * ' regexprep(strtrim(unitStr),exp,rep)];