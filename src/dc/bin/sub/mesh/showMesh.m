function showMesh(p,t),

% FIND REGIONS
RigidElements = find(t(:,4) == 0);
FlexElements  = find(t(:,4) == 1);
BatElements   = find(t(:,4) == 2);

% CALCULATE POINTS ON OTHER HALF THE WING
p2 = p;  
p2(:,1) = -p2(:,1);

% INITIALIZE FIGURE
figure
hold on
axis equal
axis off
% RIGID REGION
patch('faces',t(RigidElements,1:3),'vertices',p,'facecolor','k','edgecolor','g'); 
patch('faces',t(RigidElements,1:3),'vertices',p2,'facecolor','k','edgecolor','g','handlevisibility','off'); 
% FLEXIBLE REGION
patch('faces',t(FlexElements,1:3),'vertices',p,'facecolor','w','edgecolor','k'); 
patch('faces',t(FlexElements,1:3),'vertices',p2,'facecolor','w','edgecolor','k','handlevisibility','off'); 
% BATTON REINFORCED REGION
patch('faces',t(BatElements,1:3),'vertices',p,'facecolor','k','edgecolor','r'); 
patch('faces',t(BatElements,1:3),'vertices',p2,'facecolor','k','edgecolor','r','handlevisibility','off'); 
% LEGEND
legend('Carbon Fiber','Latex','Battons')
