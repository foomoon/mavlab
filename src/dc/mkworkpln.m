function hg = mkworkpln(pos,w,h)
% MKWORKPLN  Make Work Plane
%
%  hg = mkworkpln(pos,w,h)
x = [0 0 0 0];
y = [-1 1 1 -1];
z = [-1 -1 1 1];

x = x+pos;
y = y*w/2;
z = z*h/2;

hg = hggroup('parent',gca,'handlevisibility','on');
h = patch(x,y,z,'facecolor','interp','facealpha',.5,'handlevisibility','on');
set(h,'zdata',z);
bdown = ['h=findall(gca,''tag'',''dcworkplane''); m=find(h==gcbo);'...
         'p=cell2mat(get(h,''parent'')); pp = get(h(m),''parent'');'...
         'n = find(p==pp);  [hh,az,el,xn]=compworkpln(p,n);'...
         'ud.hh=hh; ud.az=az; ud.el=el; ud.xn=xn;'...
         'set(gca,''userdata'',ud);'];
set(h,'buttondownfcn',bdown,'parent',hg,'tag','dcworkplane');

XData = get(h,'XData');
set(h,'CData',XData)

set(h,'facelighting','none');

% set(h,'AmbientStrength',0);
% set(h,'DiffuseStrength',0);
% set(h,'SpecularStrength',0);
% set(h,'SpecularExponent',1);
% set(h,'SpecularColorReflectance',1);