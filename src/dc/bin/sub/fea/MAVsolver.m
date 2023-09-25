function [W,tx,ty,Q] = MAVsolver(K,F,bc),


%% Solve Equations
[answer,Q] = solveq(K,F,bc);

% Pre-allocate memory for output (Note, will have to change X)
n = [length(F)/3, 1];
W = zeros(n);
tx = zeros(n);
ty = zeros(n);

% Find Displacements/Rotations
answer = reshape(answer,3,size(answer,1)/3)';
W = answer(:,1)*100/2.54;
tx = answer(:,2);
ty = answer(:,3);