function [s] = num2str(v, varargin)

[dispVal,~,unitStr] = displayparser(v);
s = num2str(dispVal, varargin{:});

if isempty(s)
    s = '[]';
end

s = strcat(s,[' ' unitStr]);