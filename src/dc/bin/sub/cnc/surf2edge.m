function [e1 e2] = surf2edge(X,Y,Z)

global dcDebug
% ind = convhull(X,Y);
% x = X(ind);
% y = Y(ind);
% z = Z(ind);
x1 = fliplr(X(1,:)); x2 =  X(:,1)'; x3 = X(end,:); x4 = flipud(X(:,end))';
y1 = fliplr(Y(1,:)); y2 =  Y(:,1)'; y3 = Y(end,:); y4 = flipud(Y(:,end))';
z1 = fliplr(Z(1,:)); z2 =  Z(:,1)'; z3 = Z(end,:); z4 = flipud(Z(:,end))';
x = [x1(:); x2(:); x3(:); x4(:)];
y = [y1(:); y2(:); y3(:); y4(:)];
z = [z1(:); z2(:); z3(:); z4(:)];

[m,mini] = min(x);
[m,maxi] = max(x);

mm = sort([mini maxi]);
n = length(x);

edge1 = mm(1):mm(2);
edge2 = [mm(2):n 1:mm(1)];

c = [x(:) y(:) z(:)];
e1 = c(edge1,:);
e2 = c(edge2,:);

% If not left 2 right
if sum(diff(e1(:,1))) < 0,
    e1 = flipud(e1);
end
if sum(diff(e2(:,1))) < 0,
    e2 = flipud(e2);
end

% If e2 is on top, swtich e1,e2
if max(e1(:,2)) < max(e2(:,2))
    temp = e1;
    e1 = e2;
    e2 = temp;
end

% If end points are not the same Y value
% if e1(1,2) ~= e1(end,2)
%     if max(e1(:,2)) > max(e2(:,2))
%         e1(end,:) = [];
%         e2(1,:) = [];
%     else
%         e1(1,:) = [];
%         e2(end,:) = [];
%     end
% end
if e1(end,2) > e1(1,2)
    e1(1,:) = [];
elseif e1(end,2) < e1(1,2)
    e1(end,:) = [];
end
if e2(end,2) < e2(1,2)
    e2(1,:) = [];
elseif e1(end,2) > e1(1,2)
    e2(end,:) = [];
end


% Get rid of repeated X-points
e1([false; (diff(e1(:,1)) == 0)],:) = [];
e2([false; (diff(e2(:,1)) == 0)],:) = [];


if dcDebug
    fprintf(1,'(%s)\n',mfilename)
    % % For Debugging
    % figure
    % line(e1(:,1),e1(:,2),e1(:,3),'color','r','marker','.')
    % line(e2(:,1),e2(:,2),e2(:,3),'color','g','marker','.')
    % legend('edge 1','edge 2')
    % grid on
end

