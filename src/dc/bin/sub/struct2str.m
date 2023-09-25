function str = struct2str(s)
% STRUCT2STR  Convert structure array to string
f = fieldnames(s);
n = length(f);

for i=1:n,
    temp = getfield(s,f{i});
    a = length(temp) > 1;
    b = isstruct(temp);
    c = iscell(temp);
    d = isempty(temp);
    if a | b | c | d
        temp = nan;
    end
%     temp
    val(i,1) = temp;
end

str1 = char(f);
str1(:,end+1) = ' ';
str1(:,end+1) = '=';
str1(:,end+1) = ' ';
str1(:,end+1) = ' ';
str2 = num2str(val);
for i=1:size(str2,1)
    str2(i,:) = strrep(str2(i,:),'NaN','...');
end

str = [str1 str2];