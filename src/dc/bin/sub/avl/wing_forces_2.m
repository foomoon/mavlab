function [F] = wing_forces_2(trinodes,trielements,dCP,velocity)


%% Nodal information:
X = trinodes(:,1);
Y = trinodes(:,2);
Z = trinodes(:,3);
tri = trielements(:,1:3);

%% Number of Nodes/Elements
N_nodes = length(X);
N_elements = length(tri);

%% Find Global DOF
GlobalDOF = N_nodes*3;

%% Assign 3 DOF for each node
DOF = 1:GlobalDOF;
DOF = reshape(DOF,3,N_nodes)';


%% intitialize the force vector
F = zeros(GlobalDOF,1);
Fnorm = zeros(size(X));


%% calculate the pressure difference over the wing, using dCP

dP = .5 * 1.229 * (velocity^2) * dCP;
        
for i = 1:N_elements
      
    X1 = X(trielements(i,2))*2.54/100;
    Y1 = Y(trielements(i,2))*2.54/100;
  
    X2 = X(trielements(i,3))*2.54/100;
    Y2 = Y(trielements(i,3))*2.54/100;
  
    X3 = X(trielements(i,4))*2.54/100;
    Y3 = Y(trielements(i,4))*2.54/100;
    
    A = polyarea([X1;X2;X3],[Y1;Y2;Y3]);
    
    Fm = A*(dP(trielements(i,2)) + dP(trielements(i,3)) + dP(trielements(i,4))) / 3;
    
    Fnorm(trielements(i,2)) = Fnorm(trielements(i,2)) + Fm/3;
    Fnorm(trielements(i,3)) = Fnorm(trielements(i,3)) + Fm/3;
    Fnorm(trielements(i,4)) = Fnorm(trielements(i,4)) + Fm/3;
    
end

F(DOF(:,1),1) = Fnorm(:);

    