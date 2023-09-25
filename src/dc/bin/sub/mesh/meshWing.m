function [cost,p,t,wing] = meshWing(wing),

%% INITIAL PARAMETERS
s = wing.span;
c = wing.chord;
w = wing.layup.h0;
k = w;

%% NODES THAT MUST BE INCLUDED IN MESH
% (1) DEFINITELY NEED BOUNDARY
boundary = wing.layup.boundary;
fix1 = boundary;
% (2) INCLUDE ANY NODES ALONG MATERIAL CHANGE BOUNDARIES (IE: LATEX/CARBON)
fix2 = wing.layup.flex;
% (3) INCLUDE ANY NODES ON BATTON BOUNDARIES
b = wing.layup.batton;
try
    [fix3] = fixNodes(b,k);
catch
    fix3 = [];
end
% (*) INCLUDE ALL FIXED NODES
fix = [fix1; fix2; fix3];

%% CONSTRUCT MESH PARAMETERS
hdata = [];
% CONCENTRATE MESH NEAR BATTONS
for i=1:size(wing.layup.batton,1),
    hdata(i).type  = 2;
    hdata(i).slope = .15;
    hdata(i).ref   = [b(i,[1 3]) b(i,[2 4])];
    hdata(i).h0    = k;
    hdata(i).lim   = max(s/2,c);
end
next = length(hdata);
% CONCENTRATE MESH NEAR FLEXIBLE REGION
for i=1:size(fix2,1)-1,
    ni = next + i;
    hdata(ni).type  = 2;
    hdata(ni).slope = .15;
    hdata(ni).ref   = [fix2(i,1) fix2(i,2) fix2(i+1,1) fix2(i+1,2)];
    hdata(ni).h0    = k;
    hdata(ni).lim   = max(s/2,c);
end
if isempty(hdata),
    hdata.type = 3;
    hdata.h0   = k;
end

%% CREATE MESH
n = size(boundary,1);
connectivity = [[1:n]' [2:n 1]'];
[p,t,cost] = mesh2d(boundary,connectivity,hdata,fix);

%% FIND EACH TYPE OF ELEMENT
% FIND FLEXIBLE REGION
FlexElements = findFlex(p,t,fix2);
% FIND BATTONS
BatElements  = findBat(p,t,b,w);

%% SET 4TH COLUMN OF TRIANGLE VECTOR TO ELEMENT TYPE
% Type 0 = Carbon Fiber
% Type 1 = Flexible Membrane
% Type 2 = Batton Reinforcement
t(:,4) = 0;
t(FlexElements,4) = 1;
t(BatElements,4) = 2;

%% INTERP ZDATA
wing.mirror = 0;
[Y,X,Z] = buildWing(wing);
p(:,3) = fitgrid(X,Y,Z,p(:,1),p(:,2),'cubic');


%% SUBROUTINES
function [out] = fixNodes(Bnodes,h0),
Nx = Bnodes(:,1:2);
Ny = Bnodes(:,3:4);
Dx = diff(Nx,1,2);
Dy = diff(Ny,1,2);
L  = sqrt(Dx.^2 + Dy.^2);
step   = h0./L;
nNodes = round((1./step) + 1);
X = [];
Y = [];
for i=1:size(Bnodes,1),
    x = linspace(Nx(i,1),Nx(i,2),nNodes(i))';
    y = linspace(Ny(i,1),Ny(i,2),nNodes(i))';
    X = [X; x];
    Y = [Y; y];
end
out = [X Y];


% function n = findBattons(b,p),
% k = 1;
% for i=1:length(b),
%     for j=1:length(b{i,1}),
%         n(k) = find(p(:,1) == b{i,1}(j) & p(:,2) == b{i,2}(j));
%         k = k + 1;
%     end
% end


function F = findFlex(p,t,fix),
if isempty(fix),
    F = [];
    return
end
[in on] = inpolygon(p(:,1),p(:,2),fix(:,1),fix(:,2));
onFlex = find(in+on); 
A = ismember(t,onFlex);
F = find(sum(A,2) == 3);


function F = findBat(p,t,b,w),
% THIS METHOD IS BOGUS
k = 0.2*w;
in = 0;
on = 0;
n = size(b,1);
for i=1:2:n,
    X = [b(i,1:2)-k b(i+1,1:2)+k b(i,1)-k];
    Y = [b(i,3)+k b(i,4)-k b(i+1,3)-k b(i+1,4)+k b(i,3)+k];
    [in0 on0] = inpolygon(p(:,1),p(:,2),X,Y);
    in = in + in0;
    on = on + on0;
end

onBat = find(in+on); 
A = ismember(t,onBat);
F = find(sum(A,2) == 3);
