function triSurfWing(p,t),


tc = t(find(t(:,4) == 0),1:3);
tf = t(find(t(:,4) == 1),1:3);
tb = t(find(t(:,4) == 2),1:3);


trisurf(tc,p(:,1),p(:,2),p(:,3),'facecolor','k','edgecolor','none');
hold on
trisurf(tf,p(:,1),p(:,2),p(:,3),'facecolor',[.9 .9 .9],'edgecolor','none','facealpha',.8);
trisurf(tb,p(:,1),p(:,2),p(:,3),'facecolor','k','edgecolor','none');

axis off
axis equal


% [in on] = inpolygon(p(:,1),p(:,2),p(tc,1),p(tc,2));
% % in = abs(in - on);
% % in = find(~in);
% % p(in,3) = nan;
% 
% maxZ = max(p(:,3));
% shiftZ = 1*maxZ; 
% img = imread('cf090.jpg');
% h = textureSurf(p(:,1),p(:,2),p(:,3)+shiftZ,img);
% 
% 
% x = get(h,'Xdata');
% y = get(h,'Ydata');
% z = get(h,'Zdata');
% 
% theta = linspace(0,2*pi,100);
% xv = cos(theta) + 2;
% yv = sin(theta) - 3;
% 
% [in on] = inpolygon(x,y,xv,yv);%p(tf,1),p(tf,2));
% pts = in & ~on;
% z(pts) = nan;
% set(h,'Zdata',z);
% 
% % l = light('Position',[-50 -15 29]);