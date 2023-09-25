function mav2avl(fid,wing,varargin)

N = varargin{3};
Nspan = N(1);
NbHalf = round(Nspan/2);
wing.mirror = 0;
[Y,X,Z] = dcBuild(wing,NbHalf);
X = -X;

mat2avl(fid,X,Y,Z,varargin{:});