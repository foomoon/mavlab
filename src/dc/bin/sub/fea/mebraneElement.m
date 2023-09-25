function K3 = mebraneElement(T,X1,X2,X3,Y1,Y2,Y3),

% Find the area of each element    
A = polyarea([X1;X2;X3],[Y1;Y2;Y3]);
    
Bw = [Y2-Y3,Y3-Y1,Y1-Y2;...
      X3-X2,X1-X3,X2-X1];

Kw = (T/(4*A))*Bw'*Bw;

K3 = zeros(9,9);

K3(1,:) = [Kw(1,1),0,0,Kw(1,2),0,0,Kw(1,3),0,0];
K3(4,:) = [Kw(2,1),0,0,Kw(2,2),0,0,Kw(2,3),0,0];
K3(7,:) = [Kw(3,1),0,0,Kw(3,2),0,0,Kw(3,3),0,0];

