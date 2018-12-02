function v1 = subsasgn(v1,S,v2)
% v1 = subsasgn(v1,ind,v2) is the same as v1(ind) = v2. v1 is a DimVar. v2 must
% match the units of v1 -OR- v2 can be undimensioned NaN or []. If v1 is not a
% DimVar, this method is only called when using the subsasgn(...) syntax.
% Otherwise there is no error checking, and assignment will use the builtin
% subsasgn and, for example, convert v2 to double, so be careful.

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

if isempty(v1) && ~isa(v1,'DimVar')
    % Covers new variable case. I do not know any way to detect the difference
    % between a new variable and simply v1 = [].
    v1 = []*scd(DimVar(v2.exponents,1),v2.customDisplay);
end

if isa(v2,'DimVar')
    if isequal(v2.exponents,v1.exponents)
            v1.value(S.subs{:}) = v2.value;
    else
        error('Assignement of DimVar to DimVar must have consistent units.')
    end
    
elseif isempty(v2) || all(isnan(v2))
    % Allow deletion or assigning non-DimVar NaN without unit matching.
    v1.value(S.subs{:}) = v2;
else
    error('Assignment must be NaN, [], or consistent DimVar.')
end

end