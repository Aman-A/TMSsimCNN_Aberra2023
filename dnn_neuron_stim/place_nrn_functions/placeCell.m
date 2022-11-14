% Input 1x3 vector3 of cell location/normal:
% cell = [x,y,z], cell_normal = [Nx,Ny,Nz]
% and Npoints x 3 vector of cell morphology coordinates in local
% coordinates (z-oriented somato-dendritic axis)
% and phi in degrees
% Outputs Cnew, coordinates translated and rotated to cell position and
% aligned to normal vector with random azimuthal orientation
function [Cnew,phi] = placeCell(cell_origin,cell_normal,C,phi)             
%% Rotate about z axis (aximuthal angle)
if nargin < 4 % Use random phi if no phi input
    phi = rand(1)*2*pi; % random number between 0 and 2*pi
    Rz = [cos(phi) -sin(phi) 0; % convert from degrees to radians
       sin(phi) cos(phi) 0;
       0 0 1];     
    Cnew = (Rz*C')';
elseif isempty(phi) || phi == 0 
    Cnew = C; % don't rotate cell
else        
    Rz = [cos(phi*pi/180) -sin(phi*pi/180) 0; % convert from degrees to radians
          sin(phi*pi/180) cos(phi*pi/180) 0;
          0 0 1];     
    Cnew = (Rz*C')';        
end
%% Rotate to align somato-dendritic axis with normal
if ~isempty(cell_normal)% check that cell_normal not empty
    z = [0,0,1]'; % Default alignment of somato-dendritic axis
    R = getRmatrix(z,cell_normal); % get rotation matrix between z-axis and cell_normal
    Cnew = (R*Cnew')';
end                  
%% Translate to new origin 
if ~isempty(cell_origin) % check that cell_origin not empty    
    Cnew = Cnew + cell_origin; % translate point cloud to cell_origin, assumes
                               % soma is at [0,0,0] in Cnew
end    
end