function srf = srf2nrb(X,Y,Z)

% crv = struct;%(size(X,1),1);
for i=1:size(X,1)
    crv(i) = nrbinterp(X(i,:), Y(i,:), Z(i,:));
end


srf = nrbloft(crv);