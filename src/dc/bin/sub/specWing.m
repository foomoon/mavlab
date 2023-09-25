function stats = specWing(wing),

wing.mirror = 0;
[X,Y] = dcBuild(wing,[]);
Y = Y-max(Y(:));
% Get Border Data 
x1 = fliplr(X(1,:)); x2 =  X(:,1)'; x3 = X(end,:); x4 = flipud(X(:,end))';
y1 = fliplr(Y(1,:)); y2 =  Y(:,1)'; y3 = Y(end,:); y4 = flipud(Y(:,end))';
x = [x1 x2 x3 x4];
y = [y1 y2 y3 y4];
xl = X(:,1);
yl = Y(:,1);
xt = X(:,end);
yt = Y(:,end);

Area = 2*polyarea(x,y);
[mAC mChord] = mac(xl,yl,xt,yt,Area);
AC = abs(Y(1,1) - mAC)/wing.chord;
AR = wing.span/mChord; %(wing.span)^2/Area;

stats.PA = Area;
stats.MC = mChord;
stats.AC = AC;
stats.AR = AR;

    
        
function [ac,mean_ac] = mac(xl,yl,xt,yt,A)

% figure; plot(xl,yl,'r',xt,yt,'g')
x = linspace(xl(1),xl(end),5000)';
yl = interp1(xl,yl,x);
yt = interp1(xt,yt,x);
xl = x;

c = yl-yt;
dx = [0; diff(xl)];
% mean_ac = sum(c.^2)/sum(c);
mean_ac = 2*sum(dx.*c.^2)/A;
xb = sum(yl.*c)/sum(c);
ac = xb - (mean_ac/4);