function L = curvelength(X,Y,Z)

L = sum(sqrt(diff(X).^2 + diff(Y).^2 + diff(Z).^2));