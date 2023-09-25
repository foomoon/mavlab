function [EDOF,F,BC] = cnect(tri,p,d,f,bc)
% CNECT  Element Connectivity and Force/boundary condition assembly
N_nodes = size(p,1);
N_elements = size(tri,1);

GlobalDOF = N_nodes*d;

F = sparse(GlobalDOF,1);
BC = zeros(numel(bc),2);

%% Assign d DOF for each node 
DOF = 1:GlobalDOF;
DOF = reshape(DOF,d,N_nodes)';

EDOF = [(1:N_elements)',DOF(tri(:,1),:),DOF(tri(:,2),:),DOF(tri(:,3),:)];

if nargin > 3
    %% Initialize force vector
    F = sparse(GlobalDOF,1);

    %% Set Default Force vector
    for i=1:size(f,2)
        F(DOF(:,i),1) = f(:,i);
    end
end

if nargin > 4
    BC = DOF(bc(:,1),:)';
    BC = BC(:);
    tmp = bc(:,2:d+1)';
    BC(:,2) = tmp(:);
    rm = isnan(BC(:,2));
    BC(rm,:) = [];
end