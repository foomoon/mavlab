function mav = img2mav(filename),
% IMG2MAV Image processor for bitmap wing geometry.

% Find edges using the Canny operator with hysteresis thresholds of 0.1
% and 0.2 with smoothing parameter sigma set to 1
im = imread(filename);
if isrgb(im),
    im = rgb2gray(im);
end
edgeim = edge(im,'canny', [0.1 0.2], 1);


% Link edge pixels together into lists of sequential edge points, one
% list for each edge contour.  Discard contours less than 10 pixels long.
[edgelist, labelededgeim] = edgelink(edgeim, 10);


% Fit line segments to the edgelists with the following parameters:
tol = 2;                            % Line segments are fitted with maximum deviation from original edge of 2 pixels.
angtol = .25;                       % Segments differing in angle by less than "angtol" radians
linkrad = 10;                       % and end points within "linkrad" pixels will be merged.


% Link Close line segments together and plot them
[m,n] = size(edgelist);
p = max([m,n]);
% figure
for i=1:p,
    [seglist, nedgelist] = lineseg2(edgelist{i}, tol, angtol, linkrad);
    [x,y]=linksegs(seglist);
%     subplot(3,3,i+3); plot(x,y,'.',[x; x(1)],[y; y(1)],'r'); title(num2str(i)); axis equal; axis off;
    mav{i} = [x,y];
end




function [x3,y3] = linksegs(seglist),

seg = [seglist(:,1:2); seglist(:,3:4)];
[seg(:,1),ind] = sort(seg(:,1));
seg(:,2) = seg(ind,2);
for i = 1:length(seg),
    mag(i,1) = norm(seg(i,:));
end
ang = atan((seg(:,1)./seg(:,2)));
[mag,ind2] = sort(mag);
mag = mag/max(mag);
ang = ang(ind2);
x = mag.*cos(ang);
y = mag.*sin(ang);

T = [0 1; -1 0];
pts = [x,y];
pts = [T*pts']';
x = pts(:,1);
y = pts(:,2);

% Slight rotational Adjustment
[maxx,maxi] = max(x);
[minx,mini] = min(x);
y = y-y(mini);
theta = atan(y(maxi)/maxx);
T = [cos(theta) sin(theta); -sin(theta) cos(theta)];
pts = [x,y];
pts = [T*pts']';
x = pts(:,1);
y = pts(:,2);

% Separate Part into Upper & Lower Halves
k = 1;
j = 1;
for i=1:length(y),
    if y(i) >= -.01,
        indpos(k) = i;
        k =  k+1;
    else
        indneg(j) = i;
        j = j+1;
    end
end
xup = x(indpos);
yup = y(indpos);
xlow = x(indneg);
ylow = y(indneg);

% Resort coordinates from left to right (x-direction)
[xup,i] = sort(xup);
yup = yup(i);
[xlow,i] = sort(xlow);
ylow = ylow(i);

% Reverse Direction of Lower Half
xlow = xlow(end:-1:1);
ylow = ylow(end:-1:1);

% Combine Halves into one continuous vector (Starting at the lower right
% and working clockwise around the part)
x2 = [xlow; xup;];
y2 = [ylow; yup;];

% Average Every two points (This should combine points who are extremely
% close together) Needs some work
for i=1:length(x2)-1,
    x3(i,1) = mean([x2(i) x2(i+1)]);
    y3(i,1) = mean([y2(i) y2(i+1)]);
end

% Rotate back to where it was
theta = -theta;
T = [cos(theta) sin(theta); -sin(theta) cos(theta)];
pts = [x3,y3];
pts = [T*pts']';
x3 = pts(:,1);
y3 = pts(:,2);

% Shift and Normalize
x3 = x3-min(x3);
mx3 = max(x3);
x3 = x3/mx3;
y3 = y3/mx3;

% Tweeking Rotations
theta = atan((y3(end)-y3(end-1))/(x3(end)-x3(end-1)));
if abs(theta) <= .05,
    T = [cos(theta) sin(theta); -sin(theta) cos(theta)];
    pts = [x3,y3];
    pts = [T*pts']';
    x3 = pts(:,1);
    y3 = pts(:,2);
end

theta = atan((y3(end)-y3(1))/(x3(end)-x3(1)))-pi/2;
if abs(theta) <= .1,
    T = [cos(theta) sin(theta); -sin(theta) cos(theta)];
    pts = [x3,y3];
    pts = [T*pts']';
    x3 = pts(:,1);
    y3 = pts(:,2);
end

% Shift and Normalize Again
x3 = x3-min(x3);
mx3 = max(x3);
x3 = x3/mx3;
y3 = y3/mx3;