function [hh,az,el,xn]=compworkpln(h,n)
% COMPWORKPLN  Complile Work Planes (Animated)
%
%   [hh,az,el,xn]=compworkpln(h,n)

fig = gcf;
ud = get(fig,'userdata');
set(ud.hsurf,'visible','off');

hline = findobj(h,'type','line');
set(hline,'linestyle',':');

% get(h(n),'children')
hcurline = findall(h(n),'type','line');
set(hcurline,'linestyle','-','hittest','on');

c0 = findall(h(n),'tag','dcworkplane');
set(c0,'edgecolor','r');
x0 = mean(get(c0,'xdata'));

steps = 15; 
freq = 4*steps;  % hz

m = length(h);
k = 1;
for i=1:m
    children = get(h(i),'children'); %findobj(h(i),'tag','dcworkplane');
    for j=1:length(children)
        hh(k) = children(j);
        xn{k} = get(hh(k),'xdata');
        d(k) = mean(xn{k})-x0;
        k = k + 1;
    end
end

m = length(hh);
for i=1:steps   
    for j=1:m
        x = xn{j} - (i/steps)*d(j);
        children = hh(j); %get(h(j),'children');
        set(children,'xdata',x);        
    end
    drawnow
    pause(1/freq);
end

% Set patch transparancy to zero
hpatch = findobj(h,'type','patch');
set(hpatch,'facealpha',0,'hittest','off');


set(c0,'edgecolor','k');

% Rotate to Front view
[az,el]=view;

daz = az - (-90);
del = el - 0;

for i=1:steps
    azn = az - (i/steps)*daz;
    eln = el - (i/steps)*del;
    view(azn,eln);
    drawnow
    pause(2/freq);
end


set(hpatch,'visible','off');
ud.hh = hh; 
ud.az = az;
ud.el = el;
ud.xn = xn;
set(gca,'userdata',ud);
