function Dp = cfHooke(ply_angles)

thickness = .35/1000;
ply_angles = ply_angles * pi / 180;
N_layers = length(ply_angles);
lower_surface = (N_layers * thickness)/-2;
zbar = zeros(size(ply_angles));
Dp = zeros(3,3);

zbar = lower_surface + (thickness/2) + [0:N_layers-1]*thickness;

Q = [10700, 338,   0;
     338,   10700, 0;
     0,     0,     1632]*1E6;

for i = 1:N_layers
    c     = cos(ply_angles(i));
    s     = sin(ply_angles(i));
    T_eps = [c*c,s*s,c*s;s*s,c*c,-c*s;-2*c*s,2*c*s,c*c-s*s];
    Qbar  = T_eps'*Q*T_eps;
    Dp    = Dp + (thickness^3/12+thickness*zbar(i)^2)*Qbar;
end

