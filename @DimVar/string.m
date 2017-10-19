function s = string(v)

s = string(cellfun(@num2str,num2cell(v),'UniformOutput',false));