function v = vertcat(v,varargin)

for i = 1:length(varargin)
    vi = varargin{i};
    if compatible(v,vi)
        v.value = [v.value; vi.value];
    end
end

% 2014-05-16/Sartorius: simplified for speed; expanded to handle unlimited
% number on inputs.

% TODO: get rid of inelegant loop.