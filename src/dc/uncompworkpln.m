function uncompworkpln(h,az,el,xn)

m = length(xn);
for i=1:m
    set(h(i),'xdata',xn{i});
    if strcmp(get(h(i),'type'),'patch')
        set(h(i),'facealpha',.4,'edgealpha',1,'hittest','on',...
            'visible','on','edgecolor','k');
    elseif strcmp(get(h(i),'type'),'line')
        set(h(i),'color','k','hittest','off','linestyle','-');
    end
end
view(az,el);