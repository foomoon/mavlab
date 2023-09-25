function step = scallop2step(hs,D),
% SCALLOP2STEP - Calculates stepover distance from the scallop height
%   step = scallop2step(hs,D) computes the stepover distance (step) from the
%   scallop height (hs) and tool diameter (D).

step = sqrt(-4*(hs^2-hs*D));