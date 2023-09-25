function nrbcrv(x,y,z)

if nargin == 2
    z = zeros(size(x));
end
ppx = csapi(x,x);
ppy = csapi(x,y);
ppz = csapi(x,z);

coefs = 