function [X,Y,dCP] = readPressure(name)%,velocity,density)

cdir = cd;
cd(fileparts(which('avl26.exe')))

% if nargin == 1
%     velocity = 1;
%     density = 1;
% end
force_file = textread([name '\' name '_forces.txt'],'%s','delimiter','\n');

% Get Number of panels
% '# Chordwise  =  Nchord   # Spanwise = Nspan   First strip  =  1'
panel_str = force_file{9};
panel_str = strrep(panel_str,'#','');
panel_str = strrep(panel_str,'=','');
panel_str = strrep(panel_str,'Chordwise','');
panel_str = strrep(panel_str,'Spanwise','');
panel_str = strrep(panel_str,'First strip','');
panels = str2num(panel_str); %#ok<ST2NM>

N_chord_panels = panels(1);
N_span_panels  = panels(2);

% N_size = N_chord_panels*N_span_panels;

X_AVL = zeros(N_span_panels,N_chord_panels);
Y_AVL = zeros(N_span_panels,N_chord_panels); 
dCP_AVL = zeros(N_span_panels,N_chord_panels); 

row_position = 30;
count = 1;

for i = 1:N_span_panels
    for j = 1:N_chord_panels
        Q = cell2mat(force_file(row_position));
        Q(1:4) = [];
        Q = str2num(Q); %#ok<ST2NM>
        X_AVL(i,j) = Q(1);
        Y_AVL(i,j) = Q(2);
        dCP_AVL(i,j) = Q(6);
        count = count+1;
        row_position = row_position+1;
    end
    row_position = row_position+11;
end

% Convert dCP to dP
% dynamic_pressure = (1/2)*density*velocity^2;
% dP = dCP_AVL*dynamic_pressure;
dCP = dCP_AVL;

% Switch X,Y because AVL has a different coordinate system
X = Y_AVL;
Y = -X_AVL;

cd(cdir)