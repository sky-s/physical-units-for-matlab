function dispdisp(o)
% DISPDISP  Verbose display of structures and objects.
%   DISPDISP(Z), where Z is a struct or object, is the same as calling disp(Z),
%   except for each struct or object field or property of Z, the output of that
%   field or property's disp method will be displayed if it normally only
%   requires one line to do so.
% 
%   Use DISPDISP to overload the disp method for your classes whose properties
%   might contain other objects.
%
%   See also struct, classdef, matlab.mixin.CustomDisplay.

% DISPDISP is a standalone utility available on the file exchange at
% https://www.mathworks.com/matlabcentral/fileexchange/48637.

% Input checking.
narginchk(1,1);

validType = isobject(o) | isstruct(o);
if ~validType || ~isscalar(o)
    try
        builtin('disp',o);
        return
    catch
        error('Input must be scalar MATLAB object or struct.')
    end
end

S = evalc('builtin(''disp'',o)');

names = fieldnames(o);
a = cellfun('length',names);
fieldWidth = 4 + max(a);
for i = 1:length(names)
    n = names{i};
    try
        val = o.(n);
    catch
        % Value retrieval error.
    end
    
    if isobject(val) || isstruct(val)
        s = evalc('disp(val)');
        
        if isscalar(regexp(s,'[\r\n]+')) % If it's a single (populated) line.
            % Clean up extra line breaks and white space:
            if isstruct(val)
                s = strtrim(regexprep(s,'\s+',' '));
            else
                s = strtrim(s);
            end
            
            leader = sprintf(['%' num2str(fieldWidth) 's: '],n);

            % Replace leader and rest of line up to line break.
            S = regexprep(S,[leader '.*\n??'],[leader s],'dotexceptnewline');
        end
    end
end
disp(S(1:end-1)) % Remove last \n.
end

% Copyright (c) 2016, Sky Sartorius
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
