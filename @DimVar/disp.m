function disp(v)
[val,~,appendStr] = displayparser(v); %#ok<ASGLU>

t = evalc('disp(val)');
if isempty(t)
    t = sprintf('\t[]');
end

switch get(0,'FormatSpacing')
    case 'compact'
        fprintf('%s\t%s\n',deblank(t),appendStr);
    case 'loose'
        fprintf('%s\t%s\n\n',deblank(t),appendStr);
end