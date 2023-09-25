function [srf]=nrbloft(crv)
% Loft Univariate NURBS curves into a NURBS surface

nn = length(crv);

% ensure both curves have a common degree
crv = commonDegree(crv);

% merge the knot vectors, to obtain a common knot vector
crv = mergeKnots(crv);

% Do knot insertion where necessary in each curve
for i=1:nn
    coefs(:,:,i) = crv(i).coefs;
end
uknots = crv(1).knots;

a.coefs = coefs;
[x,y,z]=nrbnet(a);

for i=1:size(x,1)
    v(i) = nrbinterp(x(i,:),y(i,:),z(i,:));
end
% ensure both curves have a common degree
v = commonDegree(v);

% merge the knot vectors, to obtain a common knot vector
v = mergeKnots(v);

vknots = v(1).knots;

for i=1:length(v)
    c(:,:,i) = v(i).coefs;
end

srf = nrbmak(c,{vknots uknots});
srf = nrbtransp(srf);


return





function crv = commonDegree(crv)
nn = length(crv);
% d = max(cell2vec({crv.order}));
for i=1:nn
    d(i) = crv(i).order;
end
d = max(d);
for i=1:nn,
    degree_raised = d - crv(i).order;
    if degree_raised > 0
        crv(i) = nrbdegelev(crv(i), degree_raised);
    end
end


function crv = mergeKnots(crv)
nn = length(crv);
k = ({crv.knots});
ku = unique(cell2vec(k));
n = length(ku);
ka = cell(1,nn);
for i = 1:n
    for j=1:nn
        ii(j) = length(find(k{j} == ku(i)));
    end
    m = max(ii);
    for j=1:nn
        ka{j} = [ka{j} ku(i)*ones(1,m-ii(j))];
    end
end

% Do knot insertion where necessary in each curve
for i=1:nn
    crv(i) = nrbkntins(crv(i), ka{i});
end



function v = cell2vec(c)
n = numel(c);
v = [];
for i=1:n
    t = c{i};
    t = t(:);
    v = [v; t];
end






