function [varargout] = meshgrid(varargin)

for ii = 1:nargin
    I{ii} = [',varargin{' int2str(ii) '}']; %#ok<*AGROW>
    O{ii} = [',varargout{' int2str(ii) '}'];
    if isa(varargin{ii},'DimVar')
        varargout{ii} = varargin{ii};
        I{ii} = [I{ii} '.value'];
        O{ii} = [O{ii} '.value'];
    end
end

% get rid of extra comma
I{1}(1) = '';
O{1}(1) = '';

eval(['[' strcat(O{:}) '] = meshgrid(' strcat(I{:}) ');']);