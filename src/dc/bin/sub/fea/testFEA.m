% function testFEA
clear all


tic
%% Load Test Data (half of the wing)
load vjstuff                          % load test
wing = VJ;

%% Build Half the wing
wing.mirror = 0;
[X,Y,Z] = dcBuild(wing,[]);

%% Start Timer
time(1) = toc;

%% Find the boundary of the wing half
npts = 100;
boundary = getborder(X,Y,Z,npts);
% figure; line(boundary(:,1),boundary(:,2),boundary(:,3),'marker','.')

%% Mesh half the wing
boundary(end,:) = [];
npts = size(boundary,1);
connectivity = [(1:npts)' [2:npts 1]'];
hdata.h0 = curvelength(boundary([1 2],1),boundary([1 2],2),boundary([1 2],3));
hdata.type = 3;
boundary = boundary(:,1:2);
[p,t] = mesh2d(boundary,connectivity,hdata,boundary);

%% Interpolate Z data onto 2D mesh
z = griddata(X,Y,Z,p(:,1),p(:,2),'v4');
p(:,3) = z;


time(2) = toc;

%% Reorder Mesh to optimize solver
[p,t] = renumber(p,t);

time(3) = toc;




%% Assign all elements as bi-directional carbon fiber
t(:,4) = 0;

%% Convert inches to meters
p = p*25.4/1000;

%% Show Mesh elements types
% showMesh(p,t);

%% Show Element Numbering
% xc = xc(i)*25.4/1000;
% yc = yc(i)*25.4/1000;
% figure; trisurf(t(:,1:3),p(:,1),p(:,2),p(:,3)*0,'facecolor','w')
% view(0,90)
% line(xc,yc,xc*0,'marker','o','markerfacecolor','r','linestyle','none');
% text(xc,yc,xc*0,num2str([1:length(xc)]'),'fontsize',5);
% return






%% Run AVL
name = 'testFEA';
Nchord = 8;
Nspan = 12;
load runData
runData.velocity = 13;
avl(name,wing,runData,Nspan,Nchord);
clc;

%% Read in dCP from AVL
[xAVL,yAVL,dCPo] = readPressure(name);
xAVL = xAVL*2.54/100;
yAVL = (yAVL + wing.chord)*2.54/100;

time(4) = toc;

%% Create artificial dCP
% n_nodes = size(p,1);
% dCP = (0.01/n_nodes)*ones(n_nodes,1);

%% Interpolate pressure onto mesh
dCP = griddata(xAVL,yAVL,dCPo(1:Nspan,1:Nchord),p(:,1),p(:,2),'v4');
% figure; 
% trisurf(t(:,1:3),p(:,1),p(:,2),p(:,3),dCP,'facecolor','interp');
% title('Delta Cp')
% axis equal
% colorbar

%% Convert pressure coefficients into pressure
V = runData.velocity;
rho = runData.density;
dynamPress = (1/2)*rho*V^2;
dP = dCP*dynamPress;

%% Convert Pressure to force
Fn = p2force(p,t,dP);

time(5) = toc;

%% Find all nodes along the root of the wing (X == Xmin)
%% Set Default boundary conditions
bc = find(p(:,1)==min(p(:,1)));
bc = bc(:);
bc(:,2) = 0;   % Z Displacement Condition
bc(:,3) = 0;   % X-axis Rotational Condition
bc(:,4) = 0;   % Y-axis Rotational Condition

[EDOF,F,bc] = cnect(t,p,3,Fn(:),bc);

time(6) = toc;

%% Build stiffness matrix
K = assemble(p,t,EDOF);
time(7) = toc;

%% Solve Equations
% prm=symrcm(K);
% K = K(prm,prm);
% F = F(prm);
% bc(:,1) = prm(bc(:,1));
[answer] = solveq(K,F,bc);
time(8) = toc;

%% Find Displacements/Rotations (3 Degrees of Freedom / node)
n = size(answer,1)/3;
answer = reshape(answer,3,n)';
W = answer(:,1);
tx = answer(:,2);
ty =  answer(:,3);


%% Plot deformation
figure
showDeform(p,t,W,W);
title(['Total Force: ' num2str(2*sum(Fn)) ' N']);
figure
showDeform(p,t,tx,W);
title('Tx');
figure
showDeform(p,t,ty,W);
title('Ty');

%% Show the sparcity of the stiffness matrix
figure
spy(K)
title('Global Stiffness Matrix')
xlabel(['Bandwidth: ' num2str(bandwidth(K))]);

dt = [time(1) diff(time) time(end)];
fprintf(1,'\n Load: %1.5f\n Mesh: %1.5f\n Renumber: %1.5f\n AVL: %1.5f\n Forces: %1.5f\n BCs: %1.5f\n Assemble: %1.5f\n Solve: %1.5f\n Total: %1.5f\n',dt)