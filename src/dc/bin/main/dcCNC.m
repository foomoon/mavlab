function dcCNC(wing)

if ~isfield(wing,'X') | isempty(wing.X)
    wing = dcBuild(wing);
end

X = wing.X;
Y = wing.Y;
Z = wing.Z;

[wing.X,wing.Y,wing.Z] = prep4cnc(X,Y,Z);

figure('menubar','none');
dcSurf(wing);
axis equal
cnctoolbar