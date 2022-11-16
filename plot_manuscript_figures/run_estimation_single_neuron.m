% Run threshold estimation for neuron at single position 
addPaths_dnn_neuron_stim;
% data = download and unzip data...
% test dataset parameters
pos_index = 1000; % index of neuron within population (1:4999 for L2/3
cell_id = 16; % 6-10 L2/3 PCs, 11-15 L4 LBCs, 16-20 L5 PCs
cell_layer = 4; % index of layer (2 = L2/3, 3 = L4, 4 = L5)
l = 1.5; % 2 mm for L2/3 PCs, 1.5 mm for L4 LBCs and L5 PCs
E_mode = '3D_cart'; % format of E-field vectors (3D cartesian)
if cell_layer == 4
    % naming convention for weights_file is different for L5 PCs
    weights_file = sprintf(['tms_thresh2_l1.5_N9_%s_cell_E_center_cell%g_' ...
                            'train_3dcnn_rs_adr_e2k_seed1_new'],...
                            E_mode,cell_id);
else
    weights_file = sprintf(['tms_thresh2_N9_%s_cell_E_center_cell%g_' ...
                            'train_3dcnn_rs_adr_e2k_seed1'],...
                            E_mode,cell_id);
end
layer_set_num = 1; % layer set
mesh_name = 'ernie_m2m'; 
roi_name = 'handknob'; % name of region of interest within mesh
nrn_pop_ind = 1; % population index (each neuron is rotated 30 deg at each 
                  % position from population 1 to 12)
nrn_model_ver = 'maxH'; % version of neuron model (myelinated, adult human)
Efield_name = 'M1_PA_MagVenture_MC_B70_ernie';
tms_mode = 1; % 1 (PA) or -1 (AP)

%%
% Efield/neuron settings
Efield = struct(); 
Efield.E_position = pos_index;
Efield.cell_layer = cell_layer;
Efield.layer_set_num = layer_set_num;
Efield.nrn_pop_name = sprintf('nrn_pop%g',nrn_pop_ind); 
Efield.E_file = Efield_name;
Efield.interp_method = 'scattered_interp'; % 'scattered_interp' or 'simnibs_mesh_interp'
Efield.Efield_table = 'SimNIBS_TMS_Efield_sims.csv';
Efield.use_scalar_potentials = 0; 
% Sampling grid settings
sample_method_struct = struct(); 
sample_method_struct.method = 'box';
sample_method_struct.Nx = 9;
sample_method_struct.Ny = 9;
sample_method_struct.Nz = 9;
sample_method_struct.lx = l;
sample_method_struct.ly = l;
sample_method_struct.lz = l;
sample_method_struct.rshift = [0 0 0]; 

E_out = interpEfield1(cell_id,nrn_model_ver,Efield,...
                      'sample_method_struct',sample_method_struct);
E_global = E_out{1};
cell_normal = E_out{2};
phi = E_out{3};
E_loc = reorientEfield(cell_normal,phi,E_global); % rotate E-field into local coordinate system
% Get middle index for normalizing by |E| at this point 
[C,sample_method_struct] = samplePts(sample_method_struct,'print_level',1);
[~,center_ind,~] = intersect(C-sample_method_struct.rshift,[0 0 0],'rows'); 
scale_factor = norm(E_loc(center_ind,:)); % magnitude of E-field at center grid point
% Normalize to |E| of center point and reshape to NxNxNx3 tensors for input to CNN
E_in = processEfieldSamples({E_loc/scale_factor},E_mode,sample_method_struct.Nx);
dnn_output = getDNNOutput(E_in,weights_file);
% Convert threshold from |E| of center grid point (soma) back to A/us 
thresh = dnn_output/scale_factor;
fprintf('Estimated threshold = %.2f A/us (%.2f V/m)\n',thresh,dnn_output)
