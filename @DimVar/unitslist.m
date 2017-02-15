function [listString,list] = unitslist(varargin)
% UNITSLIST goes through an m-file and makes a list of units used and their
% values that can be pasted into the code (replacing the u = units; call),
% thereby significantly speeding up the code by eliminating the use of the
% DimVar class.
%
%   UNITSLIST prompts to select one or more m-files that contain u.xxx
%   calls and displays in the command window a list of units used in the
%   file(s) and their values in the form u.xxx = val; with a comment
%   describing the unit.
% 
%   UNITSLIST(fName1,fName2,fName3,...) searches the files indicated by the
%   file names provided.
% 
%   UNITSLIST(...,'-COPY') additionally copies the list to the system
%   clipboard.
% 
%   [listString] = UNITSLIST returns the list as a single string.
%
%   [listString, listCell] = UNITSLIST returns the list as a cell array of
%   strings.
%
%   UNITSLIST will only find calls if the variable name 'u' has been used
%   for the units structure.
%
%   Note that UNITSLIST can easily be used to create an alternative UNITS
%   function that returns the struct u with float values instead of
%   DimVars.
%
%   See also UNITS.

%   Copyright 2013 Sky Sartorius
%   www.mathworks.com/matlabcentral/fileexchange/authors/101715

ind = strcmpi('-copy',varargin);
if any(ind)
    copyToClipboard = true;
    varargin = varargin(~ind);
else
    copyToClipboard = false;
end

ind = strcmpi('-eval',varargin);
if any(ind)
    evalMe = true;
    varargin = varargin(~ind);
else
    evalMe = false;
end
%   Providing the '-EVAL' argument will clear the variable u and evaluate
%   the list in the caller workspace, creating a new u struct.
%
% This option isn't documented in the help block because I can't think of a
% use case for it.

if isempty(varargin)
    [fnames, pth] = uigetfile('*.m','Choose file(s)','MultiSelect','on');
    fnames = fullfile(pth,fnames);
else
    fnames = varargin;
end
if ~iscell(fnames)
    fnames = {fnames};
end

u = units; %#ok<NASGU>
list = {};
for i = 1:length(fnames)
    str = stripcomments(fileread(fnames{i}));
    newItems = regexp(str,'u\.\w+','match');
    list = [list; newItems']; %#ok<*AGROW>
end

% Get rid of duplicates in the list.
list = unique(list);

% Get values.

listStr = '';
for i = 1:length(list)
    try
        val = eval(list{i});
        if isa(val,'DimVar')
            [~,descStr] = unitsOf(val);
            val = u2num(val);
        else
            descStr = '';
        end
        list{i} = sprintf('%-11s = %-22.14g; %% %s',list{i},val,descStr);
    catch ME
        if strcmp(ME.identifier,'MATLAB:nonExistentField')
            list{i} = ... %sprintf('\b');
                sprintf('%% "%s" not defined in UNITS function.',list{i});
        else
            throw(ME)
        end
    end
    listStr = sprintf('%s\n%s',listStr,list{i});
end


if evalMe
    evalin('caller','clear u');
    evalin('caller',listStr);
end

if nargout
    listString = listStr;
else
    disp(listStr);
end

if copyToClipboard
    clipboard('copy',listStr);
    fprintf('\n***List text copied to system clipboard.***\n\n');
end
end

function ostr = stripcomments(istr)
% From http://www.mathworks.com/matlabcentral/fileexchange/4645
% MLSTRIPCOMMENTSSTR Strip comments from a string with MATLAB code.
%
%   OSTR = MLSTRIPCOMMENTSSTR(ISTR) takes the input string ISTR, assuming
%   it contains MATLAB code, strips all MATLAB comments, and returns the
%   result as the output string OSTR.  The input string must be a row
%   vector or a cell array of row vectors.  Any input string may be a
%   multi-line string, that is a string containing newline characters, for
%   example a string containing the entire contents of an m-file.
%
%   Example:  Read all MATLAB code from an m-file opened for reading on
%   file identifier "fid" and strip all comments:
%
%       str = fread(fid);               % data as numerical column vector
%       str = char(str)';               % data as character row vector
%       str = mlstripcommentsstr(str);  % strip all MATLAB comments
%
%   See also MLSTRIPCOMMENTSFILE, MLSTRIPCOMMENTSFID.

%   Author:      Peter J. Acklam
%   Time-stamp:  2004-03-18 08:23:19 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

% Regex for removing a comment on a line possibly containing MATLAB code.
%
% Use a persistent variable so that the variable is set up only once.
%
% Find the part of the line which is not a comment and then replace the
% whole line with the part that was not a comment.
%
% A comment is started by a percent sign, but not all percent signs start a
% comment.  A percent sign might be embedded in a string, in which case the
% percent sign obviously does not start a comment.  So we must look for
% percent signs which are not embedded in strings.  To do this we must find
% all strings.  Strings are delimited by single quotes, but not all single
% quotes start strings.  A single quote might also be a tranpose operator.
% So we must look for single quotes which are not tranpose operators.
%
% The MATLAB grammar (syntax and semantics) is not available to the public,
% but from my knowledge, a single quote is a transpose operator if and only
% if one of the following criteria is satisfied
%
%   1) it follows a right closing delimiter
%   2) it follows right after a variable name or function name
%   3) it follows right after a period (as part of the ".'" operator)
%   4) it follows right after another tranpose operator
%
% This can be immediately after a right closing bracket ("]"), right
% closing parenthesis (")"), right closing brace ("}"), an English letter
% ("a" to "z" and "A" to "Z", a digit ("0" to "9"), an underline ("_"), a
% period ("."), or another single quote ("'").

persistent mainregex

if isempty(mainregex)
    
    mainregex = [ ...
' (                   ' ... % Grouping parenthesis (content goes to $1).
'   ( ^ | \n )        ' ... % Beginning of string or beginning of line.
'   (                 ' ... % Non-capturing grouping parenthesis.
'                     ' ...
'' ... % Match anything that is neither a comment nor a string...
'       (             ' ... % Non-capturing grouping parenthesis.
'           [\]\)}\w.]' ... % Either a character followed by
'           ''+       ' ... %    one or more transpose operators
'         |           ' ... % or else
'           [^''%]    ' ... %   any character except single quote (which
'                     ' ... %   starts a string) or a percent sign (which
'                     ' ... %   starts a comment).
'       )+            ' ... % Match one or more times.
'                     ' ...
'' ...  % ...or...
'     |               ' ...
'                     ' ...
'' ...  % ...match a string.
'       ''            ' ... % Opening single quote that starts the string.
'         [^''\n]*    ' ... % Zero or more chars that are neither single
'                     ' ... %   quotes (special) nor newlines (illegal).
'         (           ' ... % Non-capturing grouping parenthesis.
'           ''''      ' ... % An embedded (literal) single quote character.
'           [^''\n]*  ' ... % Again, zero or more chars that are neither
'                     ' ... %   single quotes nor newlines.
'         )*          ' ... % Match zero or more times.
'       ''            ' ... % Closing single quote that ends the string.
'                     ' ...
'   )*                ' ... % Match zero or more times.
' )                   ' ...
' [^\n]*              ' ... % What remains must be a comment.
];
    
    % Remove all the blanks from the regex.
    mainregex = mainregex(~isspace(mainregex));
    
end

% Initialize output.
ostr = istr;

% If input is a character array, put it into a cell array.  We'll later make
% sure output is a character array if input is a character array.
if ischar(ostr)
    ostr = {ostr};
end

% Iterate over each element in the cell array.
for i = 1 : numel(ostr)
    
    % Get the ith input string.
    str = ostr{i};
    
    % All character arrays must be row vectors.
    if (ndims(str) > 2) | (any(size(str) > 0) & (size(str, 1) ~= 1))
        error('All character arrays must be row vectors.');
    end
    
    % Note:  With MATLAB 6.5.1 it seems that the '^' and '$' anchors in the
    % regexes '(^|\n)' and '($|\n)' are not working (bug?), so use a
    % workaround:  prepend and append a newline and remove them later.
    
    use_workaround = 1;
    
    if use_workaround
        lf   = sprintf('\n');     % LF = Line feed character
        str = [lf, str, lf];      % prepend and append an LF
    end
    
    % 2013-09-02/Sartorius 'tokenize' lines removed
    % Remove comment lines.
    str = regexprep(str, '(^|\n)[^\S\n]*%[^\n]*', '$1');
    
    % Remove a comment on a line possibly containing MATLAB code.
    str = regexprep(str, mainregex, '$1');
    
    % Remove trailing whitespace from each line.
    str = regexprep(str, '[^\S\n]+($|\n)', '$1');
    
    % Compress multiple blank lines.
    %str = regexprep(str, '\n\n\n+', sprintf('\n\n'));
    
    if use_workaround
        str = str(2:end-1);       % remove leading and trailing LF
    end
    
    % Insert ith string into output cell array.
    ostr{i} = str;
    
end

% If the input was a character array, make sure output is so too.
if ischar(istr)
    ostr = ostr{1};
end
end

% 2013-08-28/Sartorius - created
% 2013-09-02/Sartorius
%   Got rid of lots of code by using regexp
%   Integrated fileexchange/4645 to handle comments better
%   Added error catching for file-defined u fields etc.
% 2013-09-24/Sartorius added eval and copy options
