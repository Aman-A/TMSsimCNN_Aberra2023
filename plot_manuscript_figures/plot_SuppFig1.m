%% L2/3 PCs
cell_ids = 6:10;
nrn_model_ver = 'maxH';
cell_origin = [0 0 0];
cell_normal = [0 0 1];
phi = 0; 
save_fig = 1;
sample_method_struct.method = 'box';
N = 9; 
l = 2e3; % microns
sample_method_struct.Nx = N;
sample_method_struct.Ny = N;
sample_method_struct.Nz = N;
sample_method_struct.lx = l;
sample_method_struct.ly = l;
sample_method_struct.lz = l;
sample_method_struct.rshift = [0 0 -0.4868]*1e3; % microns
shift_vec = [sample_method_struct.lx*1.5 0 0]; 
plotEgrid_pos_cells(cell_ids,nrn_model_ver,cell_origin,cell_normal,...
                            phi,sample_method_struct,shift_vec,save_fig)
%% L4 LBCs
cell_ids = 11:15;
nrn_model_ver = 'maxH';
cell_origin = [0 0 0];
cell_normal = [0 0 1];
phi = 0; 
save_fig = 1;
sample_method_struct.method = 'box';
N = 9; 
l = 1.5e3; % microns
sample_method_struct.Nx = N;
sample_method_struct.Ny = N;
sample_method_struct.Nz = N;
sample_method_struct.lx = l;
sample_method_struct.ly = l;
sample_method_struct.lz = l;
sample_method_struct.rshift = [0 0 0]; 
shift_vec = [sample_method_struct.lx*1.5 0 0]; 
plotEgrid_pos_cells(cell_ids,nrn_model_ver,cell_origin,cell_normal,...
                            phi,sample_method_struct,shift_vec,save_fig)
%% L5 PCs
cell_ids = 16:20;
nrn_model_ver = 'maxH';
cell_origin = [0 0 0];
cell_normal = [0 0 1];
phi = 0; 
save_fig = 1;
sample_method_struct.method = 'box';
N = 9; 
sample_method_struct.Nx = N;
sample_method_struct.Ny = N;
sample_method_struct.Nz = N;
l = 1.5e3; % microns 
sample_method_struct.lx = l;
sample_method_struct.ly = l;
sample_method_struct.lz = l;
sample_method_struct.rshift = [0 0 0]; 
shift_vec = [sample_method_struct.lx*1.5 0 0]; 
plotEgrid_pos_cells(cell_ids,nrn_model_ver,cell_origin,cell_normal,...
                            phi,sample_method_struct,shift_vec,save_fig)