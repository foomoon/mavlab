function Db = battenHooke(N_batten_layers),

thickness = .142/1000;
lower_surface = (N_batten_layers * thickness)/-2;
zbar = zeros(1,N_batten_layers);
Db = zeros(3,3);

zbar = lower_surface + (thickness/2) + [0:N_batten_layers-1]*thickness;

% Q = [.1608,.00284,0;.00284,.1608,0;0,0,.0064]*1E12;
Q = [160.8,  2.84,   0;
     2.84,   160.8,  0;
     0,      0,      6.4]*1E9;
% Q = [6.7568e-012	1.4193e-033	0;...
%      1.4193e-033	1.0363e-010	0;...
%      0	0	2.1978e-010];


for i = 1:N_batten_layers,
    Db = Db + (thickness^3/12+thickness*zbar(i)^2)*Q;
end

