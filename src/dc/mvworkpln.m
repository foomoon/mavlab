function mvworkpln(hg,x)
% MVWORKPLN  Move Work Plane
%
%   mvworkpln(hg,x)
h = get(hg,'children');
% x0 = get(hh(1),'XData');
for i=1:length(h)
    oldx = get(h(i),'xdata');
    newx = ones(size(oldx))*x;
    set(h(i),'XData',newx);
    
    tag = get(h(i),'tag');
    if strcmp(tag,'dcworkplane')
       xdata = get(h(i),'xdata');
       set(h(i),'cdata',xdata); 
    end
end