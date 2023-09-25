function workpln2srf(fig)

if ~nargin
    fig = gcf; %findall(gcf,'type','axes');
end

% while isempty
h = findall(fig,'tag','dcworkplane');
hp = get(h,'parent');
hp = cell2mat(hp);
% i = i+1;

for i=1:length(hp)
    x(i) = getworkpln(hp(i),'position');
    crv(i) = getworkpln(hp(i),'nurb');
    xx = crv(i).coefs(1,:);
    crv(i).coefs(1,:) = 0*xx + x(i);
end

[x,ind]=sort(x);
crv = crv(ind);

srf = nrbloft(crv);
axes(get(hp(1),'parent'))
hold on;
[X,Y,Z] = nrb2xyz(srf,25*[1,1]);
ud = get(fig,'userdata');
set(ud.hsurf(1),'Xdata',X,'Ydata',Y,'Zdata',Z,'visible','on');
set(ud.hsurf(2),'Xdata',X,'Ydata',-Y,'Zdata',Z,'visible','on');
% surface(X,Y,Z,'facecolor','b','edgealpha',.1);



function [x,y,z] = nrb2xyz(nrb,res)
% NRB2XYZ  Evaluate NURBS with a given resolution

n = length(nrb.order);
if n > 1
    nu = res(1);
    nv = res(2);
    p = nrbeval(nrb,{linspace(0,1,nu),linspace(0,1,nv)});
else
    p = nrbeval(nrb,linspace(0,1,res(1)));
end
x = squeeze(p(1,:,:));
y = squeeze(p(2,:,:));
z = squeeze(p(3,:,:));
