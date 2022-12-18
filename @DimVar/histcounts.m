function [n,v,bin] = histcounts(v,varargin)
% Same as histcounts where returned edges have same units as X. 
% 
%   See also histcounts.

% [N,edges] = histcounts(X)
% [N,edges] = histcounts(X,nbins)
% [N,edges] = histcounts(X,edges)
% histcounts(...,'BinWidth',BW)
% histcounts(...,'BinLimits',[BMIN,BMAX]))
% All numeric inputs must be compatible EXCEPT if 2nd argument is a scalar.

numArgInd = cellfun(@isnumeric,varargin);
compatible(v,varargin{numArgInd});
varargin(numArgInd) = cellfun(@double,varargin(numArgInd),...
    'UniformOutput',false);
[n,v.value,bin] = histcounts(v.value,varargin{:});
end
