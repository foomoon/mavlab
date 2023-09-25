
fig = figure;
ud.hsurf(1) = surface(zeros(2),zeros(2),zeros(2),'facecolor','w',...
                   'facealpha',1,'edgealpha',.2,'edgecolor',0*[1 1 1],...
                   'facelighting','gouraud');
ud.hsurf(2) = surface(zeros(2),zeros(2),zeros(2),'facecolor','w',...
                   'facealpha',1,'edgealpha',.2,'edgecolor',0*[1 1 1],...
                   'facelighting','gouraud');
set(fig,'userdata',ud,'color','w');
light

n = 4;
xp = linspace(0,1,n);
crv(1)=nrbmak([[0 0 0 0 0 0]; -.75*[0 .1 .4 .4 .1 0]; .75*[.5 .5 .45 -.25 -.5 -.5]],[0 0 0 0 .25 .75 1 1 1 1]);

for i=1:n, 
    crv(i) = crv(1);
    crv(i).coefs(1,:)=xp(i);
    crv(i).coefs(2,:)=crv(i).coefs(2,:)/i; 
    crv(i).coefs(3,:)=crv(i).coefs(3,:)/i;
end

hold on; 

for i=1:n, 
    p = nrbeval(crv(i),linspace(0,1,100));
    x = p(1,:); y = p(2,:); z = p(3,:);
    hline(i) = line(x,y,z,'color','k');
%     hline(i) = gnurbs(crv(i));
end


for i=1:n
    hwork(i) = mkworkpln(xp(i),1,1);
end


view(-45,30); 
drawax = gca;
set(drawax,'projection','perspective');
axis equal; 
grid off; 

% hline=flipud(findall(gca,'type','line')); 

for i=1:n
    set(hline(i),'userdata',crv(i),'buttondownfcn',... 
        'nrb=get(gcbo,''userdata''); gnurbs(nrb);',...
        'hittest','off','parent',hwork(i))
end
        
% Add origin
line(0,0,0,'marker','+','markeredgecolor',[1 1 1]*.5,'markersize',10)

bdown = 'ud = get(gca,''userdata''); uncompworkpln(ud.hh,ud.az,ud.el,ud.xn); workpln2srf';
hbut = uicontrol('style','pushbutton','callback',bdown,'string','reset',...
                 'units','normalized');

set(gcf,'toolbar','figure');
axis off

ax = axes;
set(ax,'units','normalized','position',[.1,.95,.8,.05],'zlim',[-1/100,0]);

for i=1:n
    mkcontrolpt(hwork(i));
end

line([0,1],[0,0],[-1 -1]/100,'color','k','linewidth',1); 
axis off;

hold off
workpln2srf;