function K1 = plateElement(D,X1,X2,X3,Y1,Y2,Y3)
% Discrete Kirchhoff Triangular (DKT) element

% Find the area of each element    
A = triarea([X1 X2 X3],[Y1 Y2 Y3]);


K1 = zeros(9,9);

e_gauss = [0.5;0.5;0.0];
n_gauss = [0.0;0.5;0.5];

x23 = X2 - X3;  y23 = Y2 - Y3;
x31 = X3 - X1;  y31 = Y3 - Y1;
x12 = X1 - X2;  y12 = Y1 - Y2;

l23 = sqrt((x23^2 + y23^2));
l31 = sqrt((x31^2 + y31^2));
l12 = sqrt((x12^2 + y12^2));

P4 = -6*x23/(l23^2);  P5 = -6*x31/(l31^2);  P6 = -6*x12/(l12^2);
q4 = 3*x23*y23/(l23^2);  q5 = 3*x31*y31/(l31^2);  q6 = 3*x12*y12/(l12^2);
r4 = 3*y23^2/(l23^2);  r5 = 3*y31^2/(l31^2);  r6 = 3*y12^2/(l12^2);
t4 = -6*y23/(l23^2);  t5 = -6*y31/(l31^2);  t6 = -6*y12/(l12^2);

for j = 1:3

    e = e_gauss(j);
    n = n_gauss(j);

    Hxe = [P6*(1-2*e)+(P5-P6)*n, q6*(1-2*e)-(q5+q6)*n, -4+6*(e+n)+r6*(1-2*e)-n*(r5+r6), -P6*(1-2*e)+n*(P4+P6),...
        q6*(1-2*e)-n*(q6-q4), -2+6*e+r6*(1-2*e)+n*(r4-r6), -n*(P5+P4), n*(q4-q5), -n*(r5-r4)];
    Hye = [t6*(1-2*e)+(t5-t6)*n, 1+r6*(1-2*e)-(r5+r6)*n, -q6*(1-2*e)+n*(q5+q6), -t6*(1-2*e)+n*(t4+t6),...
        -1+r6*(1-2*e)+n*(r4-r6), -q6*(1-2*e)-n*(q4-q6), -n*(t4+t5), n*(r4-r5), -n*(q4-q5)];
    Hxn = [-P5*(1-2*n)-e*(P6-P5), q5*(1-2*n)-e*(q5+q6), -4+6*(e+n)+r5*(1-2*n)-e*(r5+r6), e*(P4+P6), e*(q4-q6),...
        -e*(r6-r4), P5*(1-2*n)-e*(P4+P5), q5*(1-2*n)+e*(q4-q5), -2+6*n+r5*(1-2*n)+e*(r4-r5)];
    Hyn = [-t5*(1-2*n)-e*(t6-t5), 1+r5*(1-2*n)-e*(r5+r6), -q5*(1-2*n)+e*(q5+q6), e*(t4+t6), e*(r4-r6),...
        -e*(q4-q6), t5*(1-2*n)-e*(t4+t5), -1+r5*(1-2*n)+e*(r4-r5),-q5*(1-2*n)-e*(q4-q5)];

    Bp = (1/(2*A))*[y31*Hxe+y12*Hxn;-x31*Hye-x12*Hyn;-x31*Hxe-x12*Hxn+y31*Hye+y12*Hyn];

    K1 = 2*A*Bp'*D*Bp + K1;

end

K1 = K1/6;