function plane=findcg(plane),
% FINDCG CG balancing tool.

if ~nargin 
    plane = {'','X','Y','Z','Mass';...
            'airframe',[],[],[],[];...
            'motor',[],[],[],[];...
            'battery',[],[],[],[];...
            'reciever',[],[],[],[];...
            'servo1',[],[],[],[];...
            'servo2',[],[],[],[];...
            'camera',[],[],[],[];...
            'transmitter',[],[],[],[];
            'centroids',[],[],[],[]}    
       
    % Initialize
    phrase = {'x-cg'; 'y-cg'; 'z-cg'; 'mass'};
    for i=1:size(plane,1);
        for j=1:4
            plane{i,j+1}=input([char(plane{i,1}) '.' char(phrase(j)) ': ']);
        end
    end
else
    plane
end

% Calculate total Mass
total_mass = sum(cell2mat({plane{:,5}}));
x_cg = sum(cell2mat({plane{2:end-1,2}}').*cell2mat({plane{2:end-1,5}}'))/total_mass;
y_cg = sum(cell2mat({plane{2:end-1,3}}').*cell2mat({plane{2:end-1,5}}'))/total_mass;
z_cg = sum(cell2mat({plane{2:end-1,4}}').*cell2mat({plane{2:end-1,5}}'))/total_mass;

n=size(plane,1)+1;

plane{n,1} = 'centroids';
plane{n,2} = x_cg;
plane{n,3} = y_cg;
plane{n,4} = z_cg;
plane{n,5} = total_mass;