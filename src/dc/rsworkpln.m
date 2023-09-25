function rsworkpln(hg,w,h)
% RSWORKPLN  Resize Work Plane
%
%   rsworkpln(hg,w,h)
hwork = findall(hg,'type','patch','tag','dcworkplane');

y = get(hwork,'ydata');
z = get(hwork,'zdata');
w0 = max(y)-min(y);
h0 = max(z)-min(z);

y = y*w/w0;
z = z*h/h0;

set(hwork,'ydata',y);
set(hwork,'zdata',z);