function nrb = pp2nrb(pp)

% pp.coefs = squeeze(pp.coefs);
b = pp2sp(pp);
nrb = nrbmak(b.coefs,b.knots);

