function nrbToolpath(srf,r,step,res)

hold on
% Parametric coordinates s,t
% s = linspace(0,1,res(1));
% t = linspace(0,1,res(2));
% t = 0;

% Evaluate position of actual surface at s,t
% p = nrbeval(srf,{s t});
% x = squeeze(p(1,:,:));
% y = squeeze(p(2,:,:));
% z = squeeze(p(3,:,:));

% NURBS Derivative
dsrf = nrbderiv(srf);

L = 12;
W = 7;

% Number of loops
niter = W/2/step-1;


nptsu = L/res;

% ==============
% Start For loop
for i=0:niter-1
    
    start = 0*(0 + i*step/L);
    stop  = (1 - i*step/L);
    npts = round(nptsu*(stop-start));
    s = linspace(start,stop,npts);
    p = nrbeval(srf,{s [0 1]});
    y = squeeze(p(2,:,:));
    h = abs(diff(y,[],2));
    ya = (sum(y,2)/2)';
    
    % Need to check for div/0
    h(h<eps) = 1/eps;
    t = (i*step./h)';
    p = nrbeval(srf,{s t});
    y = squeeze(p(2,:,:));
    ind = y(:,end)'<ya;
    s(ind) = [];
    t(ind) = [];

    for j=1:length(t)
        % Evaluate derivative at s,t
        [p1, dp] = nrbdeval(srf, dsrf, {s(j), t(j)});

        % Normalize Tangent/Binormal
        T = vecnorm(dp{1});
        B = vecnorm(dp{2});
        % Compute Normal via cross product and normalize it
        N = veccross(T,B);
        N = vecnorm(N);
        % Scale Normal to r
        N = N*r;

%         plot3(p1(1,:),p1(2,:),p1(3,:),'ro');
        h=quiver3(p1(1,:),p1(2,:),p1(3,:),N(1,:),N(2,:),N(3,:),0);
        set(h,'color','k');
    end
end
% End loop
% ==============
axis equal
hold off
% % figure; 
% if any(size(x))==1,
%     line(x,y,z)
% else
%     surf(x,y,z)
% end
% hold on
% plot3(p1(1,:),p1(2,:),p1(3,:),'ro');
% quiver3(p1(1,:),p1(2,:),p1(3,:),N(1,:),N(2,:),N(3,:),0);
% axis equal
% hold off