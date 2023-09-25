function P = getborder(X,Y,Z,npts)
X1 = X(1,:)';
X2 = X(:,end);
X3 = fliplr(X(end,:))';
X4 = flipud(X(:,1));

Y1 = Y(1,:)';
Y2 = Y(:,end);
Y3 = fliplr(Y(end,:))';
Y4 = flipud(Y(:,1));

Z1 = Z(1,:)';
Z2 = Z(:,end);
Z3 = fliplr(Z(end,:))';
Z4 = flipud(Z(:,1));


X = [X1; X2; X3; X4];
Y = [Y1; Y2; Y3; Y4];
Z = [Z1; Z2; Z3; Z4];

% Length of border
Ltot = curvelength(X,Y,Z);
L(1) = curvelength(X1,Y1,Z1);
L(2) = curvelength(X2,Y2,Z2);
L(3) = curvelength(X3,Y3,Z3);
L(4) = curvelength(X4,Y4,Z4);

n = ceil((L/Ltot)*npts);

P1 = curvspace([X1 Y1 Z1],n(1));
P2 = curvspace([X2 Y2 Z2],n(2));
P3 = curvspace([X3 Y3 Z3],n(3));
P4 = curvspace([X4 Y4 Z4],n(4));

P = [P1(1:n(1)-1,:); P2(1:n(2)-1,:); P3(1:n(3)-1,:); P4(1:n(4),:)];

% figure;
% line(P(:,1),P(:,2),P(:,3),'marker','.'); axis equal
% title(['npts: ' num2str(size(P,1))])
