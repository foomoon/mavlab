function A = triarea(X,Y)

% T = [X; Y; [1 1 1]];
% 
% A = (1/2)*abs(det(T));

A = 1/2*abs(X(1)*Y(2)-X(1)*Y(3)-Y(1)*X(2)+Y(1)*X(3)+X(2)*Y(3)-X(3)*Y(2));