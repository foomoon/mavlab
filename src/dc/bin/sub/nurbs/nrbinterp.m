function nrb = nrbinterp(x,y,z)

switch nargin
    case 1,
        ctrlpts = x;
        if size(ctrlpts,1) == 2
            ctrlpts(3,:) = 0;
        end
    case 2
        z = 0*x;
        ctrlpts = [x(:)'; y(:)'; z(:)'];
    case 3
        ctrlpts = [x(:)'; y(:)'; z(:)'];
end


% sp = cscvn(ctrlpts);
sp = nintrp1(ctrlpts);
nrb = pp2nrb(sp);


function sp = nintrp1(ctrlpts)
n = size(ctrlpts,2);
x = linspace(0,1,n);
sp = csapi(x,ctrlpts);