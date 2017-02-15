function [v1,I] = sort(v1,varargin)

[v1.value, I] = sort(v1.value,varargin{:});
