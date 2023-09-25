function [srf,crvs] = nrbWing(wing)

wing.mirror = 0;
x = wing.foils.x(:,1);
y = wing.foils.y(:,1);
maxy = max(y);
xx = cosdist(12);
y = spline(x,y,xx);
x = xx;
z = zeros(size(x));

% crv = nrbinterp(x,y,z);
crv = pp2nrb(cscvn([x(:)'; y(:)']));
newx = crv.coefs(1,:)';
newy = crv.coefs(2,:)';
wing.foils.x = newx;
wing.foils.y = newy;

% nrbplot(crv,100);
% hold on; plot(newx,newy,'r.')

err = (max(newy)-maxy)/max(newy);
camber_mod = wing.camber/(1-err);
wing.camber = camber_mod;
[X,Y,Z] = dcBuild(wing,10);

n = size(X,1);
% curves = struct;
for i=1:n 
    curves(i) = crv;
    curves(i).coefs(1:3,:) = [X(i,:); Y(i,:); Z(i,:)];
%     curves(i)=nrbinterp([X(i,:); Y(i,:); Z(i,:)]); 
    crvs{i} = curves(i);
end

srf = nrbloft(curves);



% 
% pts1 = [X(:,1)'; Y(:,1)'; Z(:,1)'];
% pts2 = [X(:,end)'; Y(:,end)'; Z(:,end)'];
% pts3 = [X(1,:); Y(1,:); Z(1,:)];
% pts4 = [X(end,:); Y(end,:); Z(end,:)];
% 
% order = 3;
% k     = [size(pts1,2) size(pts2,2) size(pts3,2) size(pts4,2)];
% n = order + k - 2*(order-1);
% 
% knt1 = [zeros(1,order-1) linspace(0,1,n(1)) ones(1,order-1)];
% knt2 = [zeros(1,order-1) linspace(0,1,n(2)) ones(1,order-1)];
% knt3 = [zeros(1,order-1) linspace(0,1,n(3)) ones(1,order-1)];
% knt4 = [zeros(1,order-1) linspace(0,1,n(4)) ones(1,order-1)];
% % knt1 = aptknt(pts1(1,:),n);
% % knt2 = aptknt(pts2(1,:),n);
% % knt3 = aptknt(fliplr(pts3(2,:)),n);
% % knt4 = aptknt(fliplr(pts4(2,:)),n);
% 
% crv1 = nrbmak(pts1,knt1);
% crv2 = nrbmak(pts2,knt2);
% crv3 = nrbmak(pts3,knt3);
% crv4 = nrbmak(pts4,knt4);
% 
% % figure; hold on
% % nrbplot(crv1,30)
% % nrbplot(crv2,30)
% % nrbplot(crv3,30)
% % nrbplot(crv4,30)
% % 
% srf = nrbloft(crv1,crv2,crv3,crv4);
% % hold on;
% % nrbplot(srf,[30,30])
% 
% % surf(X,Y,Z,'edgecolor','none'); axis equal
% 
% 
% 
% function srf = nrbloft(u1,u2,v1,v2)
% 
% if nargin ~= 4
%   error('Incorrect number of input arguments');
% end
% 
% r1 = nrbruled(u1, u2);
% r2 = nrbtransp(nrbruled(v1, v2));
% t  = nrb4surf(u1.coefs(:,1), u1.coefs(:,end), u2.coefs(:,1), u2.coefs(:,end));
% % This is a fix (Daniel Claxton)
% t = nrbtransp(t);
% 
% % raise all surfaces to a common degree
% du = max([r1.order(1), r2.order(1), t.order(1)]);
% dv = max([r1.order(2), r2.order(2), t.order(2)]);
% r1 = nrbdegelev(r1, [du - r1.order(1), dv - r1.order(2)]);
% r2 = nrbdegelev(r2, [du - r2.order(1), dv - r2.order(2)]);
% t  = nrbdegelev(t,  [du - t.order(1),  dv - t.order(2)]);
% 
% % merge the knot vectors, to obtain a common knot vector
% 
% % U knots
% k1 = r1.knots{1};
% k2 = r2.knots{1};
% k3 = t.knots{1};
% k = unique([k1 k2 k3]);
% n = length(k);
% kua = [];
% kub = [];
% kuc = [];
% for i = 1:n
%   i1 = length(find(k1 == k(i)));
%   i2 = length(find(k2 == k(i)));
%   i3 = length(find(k3 == k(i)));
%   m = max([i1, i2, i3]);
%   kua = [kua k(i)*ones(1,m-i1)];  
%   kub = [kub k(i)*ones(1,m-i2)];
%   kuc = [kuc k(i)*ones(1,m-i3)];
% end  
% 
% % V knots
% k1 = r1.knots{2};
% k2 = r2.knots{2};
% k3 = t.knots{2};
% k = unique([k1 k2 k3]);
% n = length(k);
% kva = [];
% kvb = [];
% kvc = [];
% for i = 1:n
%   i1 = length(find(k1 == k(i)));
%   i2 = length(find(k2 == k(i)));
%   i3 = length(find(k3 == k(i)));
%   m = max([i1, i2, i3]);
%   kva = [kva k(i)*ones(1,m-i1)];  
%   kvb = [kvb k(i)*ones(1,m-i2)];
%   kvc = [kvc k(i)*ones(1,m-i3)];
% end  
% 
% r1 = nrbkntins(r1, {kua, kva});
% r2 = nrbkntins(r2, {kub, kvb});
% t  = nrbkntins(t,  {kuc, kvc});
% 
% C=(r1.coefs(:,:,1)-r1.coefs(:,:,end))./(t.coefs(:,:,1)-t.coefs(:,:,end));
% % C([1 3 4],:) = [];
% C(isnan(C))=0;
% C = sqrt(C(1,:).^2 + C(2,:).^2 + C(3,:).^2);
% for i=1:size(r1.coefs,3)
% %     r2.coefs(3,:,i) = C.*(r2.coefs(3,:,i)-t.coefs(3,:,i))+t.coefs(3,:,i);
% % r2.coefs(3,:,i) = C.*(r2.coefs(3,:,i)-r1.coefs(3,:,i))+r1.coefs(3,:,i);
% end
% 
% % combine coefficient to construct Coons surface
% for i=1:4
%     coefs(i,:,:) = r1.coefs(i,:,:) + r2.coefs(i,:,:) - t.coefs(i,:,:);
% end
% 
% 
% srf = nrbmak(coefs, r1.knots);
