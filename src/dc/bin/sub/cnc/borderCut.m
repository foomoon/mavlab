function [X,Y,Z] = borderCut(X1,Y1,Z1,D,b),
        
bRat = b*D/2;

Z = Z1 - (D/2)*(1-cos(asin(b)));
N = curvnormal([X1 Y1]);
P = [-bRat*N] + [X1 Y1];

X = P(:,1);
Y = P(:,2);
 
X([end]) = [];
Y([end]) = [];
Z([end]) = [];

X = [X(1); X; X(end)];
Y = [Y(1); Y; Y(end)];
Z = [0; Z; 0];