function Fn = p2force(nodes,tri,dP)

M = eye(3)+1;

X = nodes(:,1);
Y = nodes(:,2);
tx = X(tri(:,1:3));
ty = Y(tri(:,1:3));
dP = dP(tri(:,1:3));

N_elements = size(tri,1);
Fn = zeros(size(X));

for i = 1:N_elements
    
    A = triarea(tx(i,:),ty(i,:));    
    Fm = (A/12)*dP(i,1:3)*M;     
    Fn(tri(i,1)) = Fn(tri(i,1)) + Fm(1);   %#ok<AGROW>
    Fn(tri(i,2)) = Fn(tri(i,2)) + Fm(2);   %#ok<AGROW>
    Fn(tri(i,3)) = Fn(tri(i,3)) + Fm(3);   %#ok<AGROW>
end

% F(DOF(:,1),1) = Fnorm(:);