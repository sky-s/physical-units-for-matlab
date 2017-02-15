function v = horzcat(v,varargin)

for i = 1:length(varargin)
    vi = varargin{i};
    if compatible(v,vi)
        v.value = [v.value vi.value];
    end
end

% 2014-05-14/Sartorius: simplified for speed; got rid of error code typo;
% expanded to handle unlimited number on inputs.

% TODO: get rid of inelegant loop.