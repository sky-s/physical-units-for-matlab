function nist(varargin)
% nist  Open National Institute of Standards and Technology Constants, Units,
% and Uncertainty website.
% 
%   nist(s), where s is a string representing the constant, will attempt to open
%   the page corresponding to that constant. A cell string or character array
%   for s or nist(s1,s2,s3,...) opens multiple pages.
% 
%   If any argument is '-search', the NIST seach results are used.
% 
%   Example: Open physics.nist.gov/cgi-bin/cuu/Value?z0 in a browser.
%       nist('Z0')
% 
%   Example: Search for Planck and Boltzmann.
%       nist Planck Boltzmann -search
%       nist(["Planck" "Boltzmann"],'-search')
% 
%   See also u.

%   Copyright Sky Sartorius
%   www.mathworks.com/matlabcentral/fileexchange/authors/101715

if ~nargin
    web('https://physics.nist.gov/cuu/Constants/index.html','-browser')
    return
end

searchFlagInd = strcmpi(varargin,'-search');

if any(searchFlagInd)
    varargin = varargin(~searchFlagInd);
    base = 'https://physics.nist.gov/cgi-bin/cuu/Results?search_for=%s';
else
    base = 'https://physics.nist.gov/cgi-bin/cuu/Value?%s';
end


for i = 1:numel(varargin)
    arg = cellstr(varargin{i});
    for j = 1:numel(arg)
        web(sprintf(base,lower(arg{j})),'-browser')
    end
end
