function h = showDeform(p,t,W,d),

% Set all triangles to be clockwise
t = triorder(t(:,1:3),p);

hold on
h(1) = trisurf(t(:,1:3), p(:,1), p(:,2), p(:,3)+d, W,...
    'facecolor','interp','edgealpha',.3,'facealpha',.6,...
    'facelighting','gouraud','edgecolor','k'); 
h(2) = trisurf(t(:,1:3), -p(:,1), p(:,2), p(:,3)+d, W,...
    'facecolor','interp','edgealpha',.3,'facealpha',.6,...
    'facelighting','gouraud','edgecolor','k');
hold off

% Axes Labels
xlabel('X')
ylabel('Y')

% Axes properties
axis equal
axis tight
box on
view(-60,35)

% Add light 
light

% Colorbar
% pause(1)
% colorbar;




function tri = triorder(tri,p)

x = p(:,1); y = p(:,2); %z = p(:,3);

tx = x(tri);
ty = y(tri);
% tz = z(tri);

% Need to re-order all triangles so they are clockwise
i = 2;
dir = (tx(:,i) - tx(:,i-1)) .* (ty(:,i+1) - ty(:,i)) - (ty(:,i) - ty(:,i-1)) .* (tx(:,i+1) - tx(:,i));
dir = -dir./abs(dir);
dir = dir == -1;
% tx(dir,:) = fliplr(tx(dir,:));
% ty(dir,:) = fliplr(ty(dir,:));
% tz(dir,:) = fliplr(tz(dir,:));
tri(dir,:) = fliplr(tri(dir,:));