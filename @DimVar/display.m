function display(v) %#ok<DISPLAY>
% This method is purely for more aesthetic display of DimVars in the command
% window. It is not necessary, as disp is already overloaded.

% Note: display used to contain important functionality with the syntax:
% [numeratorString, denominatorString, v] = display(v).

[dispVal,~,~,numString,denString] = displayparser(v);

showName = inputname(1);
if isempty(showName)
    showName = 'ans';
end

if strcmp(matlab.internal.display.formatSpacing,'loose')
    looseLine = '\n';
else
    looseLine = '';
end

% Display.
fprintf([looseLine '%s =\n' looseLine],showName)
disp(dispVal)
if isempty(denString)
    fprintf(['\t\t %s\n' looseLine], numString);
else
    numLength = length(numString);
    denLength = length(denString);
    barLength = max(numLength,denLength);
    hbar = repmat('-',1,barLength);
    
    numBuffer = repmat(' ',1,max(0,floor((denLength-numLength)/2)));
    
    denBuffer = repmat(' ',1,max(0,floor((-denLength+numLength)/2)));
    
    fprintf(['\t\t%s%s\n\t\t%s\n\t\t%s%s\n' looseLine],...
        numBuffer,numString,hbar,denBuffer,denString);
    
    % % Alternate using cprintf to underline:
    % if barLength <= numLength
    %     fprintf('\t\t%s',numBuffer);
    %     cprintf('_black','%s\n',numString);
    % else
    %     gapLength = barLength - numLength;
    %     frontPad = floor(gapLength/2);
    %     backPad = gapLength-frontPad;
    %     fprintf('\t\t');
    %     cprintf('_black',...
    %         [repmat(' ',1,frontPad) '%s' repmat(' ',1,backPad) '\n'],numString);
    % end
    % fprintf(['\t\t%s%s\n' looseLine],denBuffer,denString);
end