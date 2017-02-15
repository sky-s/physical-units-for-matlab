function v1 = subsasgn(v1,S,v2)
% v1 = subsasgn(v1,ind,v2) is the same as v1(ind) = v2. v1 is a DimVar. v2
% must match the units of v1 -OR- v2 can be undimensioned NaN or []. If v1
% is all NaN, v1's units will be changed to match those of v2 (with a
% warning), even if v2 is all NaN.

% There are 21 possible scenarios for which this will get called for the 5
% possible inputs for each of v1 and v2: 
%   [] or new
%   [] DimVar
%   DimVar
%   NaN
%   NaN DimVar
% There are 21 instead of 25 because of 4 cases that will go to normal
% subsasgn. I've done what testing I can, and the 6 most common cases all
% have the desired behavior. However, please report unexpected behavior.

if isempty(v1)
    v1 = []*unitsOf(v2);
end

if isempty(v2)
    makesub([]);
    return
end

if isa(v2,'DimVar')
    if isequal(v2.exponents,v1.exponents)
        if isempty(v2.value)
            %don't allow shrinking array without consistent units
            makesub([]);
        else
            makesub(v2.value);
        end
    elseif all(isnan(v1.value(:)))
        % initialized NaN array - throw out exponents
        v1.exponents = v2.exponents;
        warning(['Units of all-NaN assignee array "' inputname(1)...
            '" changed to match units of assignment "' inputname(3) '".'])
        makesub(v2.value)
    else
        error('Assignement of DimVar to DimVar must have consistent units.')
    end
    
elseif isempty(v2)
    makesub([]);
elseif isnan(v2) 
    % isnan returns TRUE for NaN*u.xx, hence catching ALL DimensionedVars
    % in the first nested if above
    makesub(nan);
else
    error('Assignment must be NaN, [], or consistent DimVar.')
end


    function makesub(v2)
        v1.value(S.subs{:}) = v2;
    end

end

% Revision history
%{
2013-04-12/Sartorius
   -In response to FEX user feedback ("Very cool, I use it extensively and
    it really helps. One thing that I encountered is that assignements such
    as vect(index)=[] which remove elements in a vector don't work when
    using a Dimensioned Variable vector. Also preallocating arrays e.g.
    before a loop is not possible when we don't know the unit of the
    content at the moment of allocation. Maybe NaN and 0 should be treated
    specially, in a sense these have no unit and any unit at the same
    time.")
   -I threw out the entire old function that had lots of EVAL calls and
    implemented this more elegant version. I hope it works for everything,
    but I'm not sure.
   -Changing the units of v1 for all-nan v1 might have problems, but I hope
    not.
%}