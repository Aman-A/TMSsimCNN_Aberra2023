function [C,sample_method_struct] = samplePts(sample_method_struct,varargin) 
%SAMPLEPTS Outputs coordinates of sampling points
%  
%   Inputs 
%   ------ 
%   Optional Inputs 
%   --------------- 
%   Outputs 
%   ------- 
%   C : N x 3 array
%       Coordinates of sampling points
%   sample_method_struct : struct
%                          Same as input, but missing fields are added with
%                          default values
%   Examples 
%   --------------- 

% AUTHOR    : Aman Aberra 
% get sampling points
in.print_level = 0; 
in = sl.in.processVarargin(in,varargin);
if nargin == 0   
   sample_method_struct.method = 'box';
   sample_method_struct.Nx = 9;
   sample_method_struct.Ny = 9;
   sample_method_struct.Nz = 9;
   sample_method_struct.lx = 2;
   sample_method_struct.ly = 2;
   sample_method_struct.lz = 2;
   sample_method_struct.rshift = [0 0 0]; 
   sample_method_struct.reorder = 0; 
end                 
if isfield(sample_method_struct,'rshift')    
    if isempty(sample_method_struct.rshift)
       sample_method_struct.rshift = [0,0,0];    
    elseif iscolumn(sample_method_struct.rshift)
        sample_method_struct.rshift = sample_method_struct.rshift';
    end    
else
    sample_method_struct.rshift = [0 0 0];    
end
if ~isfield(sample_method_struct,'reorder')
    sample_method_struct.reorder = 0; 
end
if strcmp(sample_method_struct.method,'box')
   Nx = sample_method_struct.Nx; 
   Ny = sample_method_struct.Ny; 
   Nz = sample_method_struct.Nz; 
   lx = sample_method_struct.lx;
   ly = sample_method_struct.ly;
   lz = sample_method_struct.lz;
   % make grid centered at [0 0 0], z along [0 0 1]
   [X,Y,Z] = meshgrid(linspace(-lx/2,lx/2,Nx),...
       linspace(-ly/2,ly/2,Ny),...
       linspace(-lz/2,lz/2,Nz));
   C = [X(:),Y(:),Z(:)];
   if sample_method_struct.reorder
       C = [0 0 0;C(~all(C == [0 0 0],2),:)]; % move (or add) [0 0 0] to first row for placeCell to use as reference
   end
   C = C + sample_method_struct.rshift;
   % translate and rotate to cell position
   %   Cs = placeCell(cell_origin,cell_normal,C,phi);
   %   Cs = [Cs(2:uint64(size(C,1)/2),:);Cs(1,:);Cs(uint64(size(C,1)/2 + 1):end,:)];  % move origin back to middle
   print_str = ['Generated sampling grid with sample_method = ''%s'',\n  Nx = %g, Ny = %g, Nz = %g, ',...
       'lx = %g, ly = %g, lz = %g\n  r_shift = [%g %g %g], reorder = %g\n'];
elseif strcmp(sample_method_struct.method,'sphere')
    % TODO finish sphere
    error('Not implemented yet'); 
elseif strcmp(sample_method_struct.method,'cylinder')          
    % Cylinder centered on cell body (r = 0, z= 0)    
    Nr = sample_method_struct.Nr; % number of radial points
    Nphi = sample_method_struct.Nphi; % number of azimuthal points
    Nz = sample_method_struct.Nz;  % number of vertical points
    lr = sample_method_struct.lr; % radius length
    lz = sample_method_struct.lz; % z length
    % make grid centered at [0 0 0], z along [0 0 1]   
    [R,PHI,Z] = meshgrid(linspace(0,lr,Nr),...
        linspace(0,2*pi-2*pi/Nphi,Nphi),... % leave out 2*pi redudant point
        linspace(-lz/2,lz/2,Nz));
    X = R.*cos(PHI); % cartesian coordinates
    Y = R.*sin(PHI);
    C = [X(:),Y(:),Z(:)];
    C = C + sample_method_struct.rshift; 
%     Cs = placeCell(cell_origin,cell_normal,C,phi);
    print_str = ['Generated sampling grid with sample_method = ''%s'',\n  Nr = %g, Nphi = %g, Nz = %g, ',...
                 'lr = %g, lz = %g\n  r_shift = [%g %g %g], reorder = %g\n']; 
end
if in.print_level > 0    
    vals = struct2cell(reformat(sample_method_struct)); 
    fprintf(print_str,vals{:});    
end
end
% ensure order of fields matches default for correct printing
function s_out = reformat(s_in)
s_out = struct();
s_out.method = s_in.method;
if strcmp(s_out.method,'box')
    s_out.Nx = s_in.Nx;
    s_out.Ny = s_in.Ny;
    s_out.Nz = s_in.Nz;
    s_out.lx = s_in.lx;
    s_out.ly = s_in.ly;
    s_out.lz = s_in.lz;
    s_out.rshift = s_in.rshift;
    s_out.reorder = s_in.reorder; 
elseif strcmp(s_out.method,'cylinder')
    s_out.Nr = s_in.Nr;
    s_out.Nphi = s_in.Nphi;
    s_out.Nz = s_in.Nz;
    s_out.lr = s_in.lr;
    s_out.lz = s_in.lz;
    s_out.rshift = s_in.rshift;
    s_out.reorder = s_in.reorder; 
end
end