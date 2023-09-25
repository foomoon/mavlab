function [K] = assemble(trinodes,trielements,varargin),
% ASSEMBLE - Compute Stiffness Matrix for Mav
%   
%   [K] = ASSEMBLE(nodes,tri,properties) returns the
%   stiffness matrix (K) for the wing geometry
%   
%   nodes      - nodal values (n-by-3)
%   tri        - element array (n-by_4) *4th column is element type

%   properties - structure array containing following fields:
%       PlyAngles - vector of layer orientation angles (1-by-#Layers) 
%                   ie [45 90 45] is three layers alternating between 45
%                   and 90 degrees
%       Tension   - Latex Pretension (0,1, or 2)
%       BatLayers - Number of layers for battons 
%


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

%% Initialize stiffness matrix
K = sparse(GlobalDOF,GlobalDOF);

%% Assign 3 DOF for each node 
DOF = 1:GlobalDOF;
DOF = reshape(DOF,3,N_nodes)';

%% Assign Nodal DOF to each triangular element
Edof = [[1:N_elements]',DOF(tri(:,1),:),DOF(tri(:,2),:),DOF(tri(:,3),:)];


% Input Check
switch nargin
    case 2,
        %% load carbon fiber layup schedule:
        ply_angles = [45 90 45];
        %% load the number uni-directional carbon fiber layers for each batten
        N_batten_layers = 3;
        %% load the latex pretension (0 - low, 1 - medium, 2 - high)
        Pretension = 0;
    case 3,
        prop = varargin{1};
        ply_angles = prop.PlyAngles;
        Pretension = prop.Tension;
        N_batten_layers = prop.BatLayers;
    otherwise,
        if nargin > 2,
            error('Incorrect Number of inputs. Please enter nodes, elements and properties')
        end
end

%% MATERIAL PROPERTIES

% Calculate the latex stress resultant
T = (1/8)*4^(Pretension + 1);

% Calculate the [0-90] bi-directional carbon fiber D matrix
Dp = cfHooke(ply_angles);

% Calculate the uni-directional D matrix
Db = battenHooke(N_batten_layers);


%% Calculate the coordinates of the three nodes of each element
X1 = X(trielements(:,1));
Y1 = Y(trielements(:,1));
% Z1 = Z(trielements(:,1));

X2 = X(trielements(:,2));
Y2 = Y(trielements(:,2));
% Z2 = Z(trielements(:,2));

X3 = X(trielements(:,3));
Y3 = Y(trielements(:,3));
% Z3 = Z(trielements(:,3));

%% Label each element as carbon fiber (0), latex (1), or batten (2)    
element_type = trielements(:,4);


%% Assemble Global Stiffness matrix from element Stiffness matrices
for i = 1:N_elements,      
    
    switch element_type(i),        
        
        case 0,
        %% Calculate stiffness matrix of each carbon fiber element
        K0 = plateElement(Dp,X1(i),X2(i),X3(i),Y1(i),Y2(i),Y3(i));
        K = assem(Edof(i,:),K,K0);    
        
        case 1,
        %% Calculate stiffness matrix of each latex element        
        K1 = membraneElement(T,X1(i),X2(i),X3(i),Y1(i),Y2(i),Y3(i));
        K = assem(Edof(i,:),K,K1);
        
        case 2,
        %% Calculate stiffness matrix of each batten element
        K2 = plateElement(Db,X1(i),X2(i),X3(i),Y1(i),Y2(i),Y3(i));
        K = assem(Edof(i,:),K,K2);
         
    end        
    
end
