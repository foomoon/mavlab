function hs = step2scallop(step,D),
% STEP2SCALLOP - Calculates scallop height from stepover distance
%   hs = step2scallop(step,D) computes the scallop height (hs) from the
%   stepover distance (step) and tool diameter (D).

hs = (1/2)*(D - sqrt(D^2 - step^2));