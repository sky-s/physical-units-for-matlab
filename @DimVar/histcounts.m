function [n,v,bin] = histcounts(v,varargin)
% Same as histcounts where returned edges have same units as X. 
% 
%   See also histcounts.

% [N,edges] = histcounts(X)
% [N,edges] = histcounts(X,edges)
% histcounts(...,'BinWidth',BW)
% histcounts(...,'BinLimits',[BMIN,BMAX]))

numArgInd = cellfun(@isnumeric,varargin);

% All numeric inputs must be compatible EXCEPT if 2nd argument is a scalar.
if nargin > 1 && isnumeric(varargin{1}) && isscalar(varargin{1})
    % [N,edges] = histcounts(X,nbins)
    compatible(v,varargin{numArgInd(2:end)});
else
    compatible(v,varargin{numArgInd});
end

varargin(numArgInd) = cellfun(@double,varargin(numArgInd),...
    'UniformOutput',false); % double returns v.value for DimVars.
[n,v.value,bin] = histcounts(v.value,varargin{:});
end
