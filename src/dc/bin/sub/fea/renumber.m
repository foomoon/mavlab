function [p,t] = renumber(p,t)

% Reorder Mesh to optimize solver
[m,n] = size(t);
x = reshape(p(t,1),m,3);
y = reshape(p(t,2),m,3);
xo = 0; yo = 0;
xc = 1/3*sum(x,2);
yc = 1/3*sum(y,2);
R = sqrt((xc-xo).^2 + (yc-yo).^2);
[i,i] = sort(R);
% i = randperm(length(i));  % Randomly choose order (Really bad idea)
t = t(i,:);
told = t;
t = t';
count = 1;
for j=1:numel(t),
    if t(j) >= count
        t(t==count) = count+length(p);
        t(t==t(j)) = count;
        count = count + 1;
    end
end
t = t';
p(t,:) = p(told,:);
t = triorder(t,p);


function tri = triorder(tri,p)

x = p(:,1); y = p(:,2); 

tx = x(tri);
ty = y(tri);

cw = -1;
ccw = 1;

% Need to re-order all triangles so they are counter-clockwise
i = 2;
dir = (tx(:,i) - tx(:,i-1)) .* (ty(:,i+1) - ty(:,i)) - (ty(:,i) - ty(:,i-1)) .* (tx(:,i+1) - tx(:,i));
dir = -dir./abs(dir);
dir = dir == ccw;
tri(dir,:) = fliplr(tri(dir,:));