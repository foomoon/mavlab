function wing = mkWing

wing.span     = 1; % Wing Span
wing.chord    = 1; % Root Chord
wing.camber   = 0; % Maximum Camber [0..1]
wing.tip      = 0; % Scale Z position as function of span [0..1]
wing.twist    = 0; % Spanwise twist in degrees
wing.sweep    = 0; % Quarter chord sweep in degrees
wing.dihedral = 0; % Dihedral angle in degrees
wing.edgeref  = 0; % Use Leading edge as Z-reference or max camber as Z-ref (0,1)
wing.mirror   = 1; % Mirror wing (0,1)

wing.planform.Ledge.x = linspace(0,1,20)';
wing.planform.Ledge.y = ones(20,1);
wing.planform.Ledge.z = zeros(20,1);
wing.planform.Tedge.x = linspace(0,1,20)';
wing.planform.Tedge.y = -ones(20,1);
wing.planform.Tedge.z = zeros(20,1);

wing.foils.x = linspace(0,1,20)';
wing.foils.y = zeros(20,1);

wing.ind = [];

wing.X = [];
wing.Y = [];
wing.Z = [];