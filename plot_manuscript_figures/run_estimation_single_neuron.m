% Run threshold estimation for neuron at single position 
mat_dir = addPaths_dnn_neuron_stim;
% data = download and unzip data...
% test dataset parameters
pos_index = 1000; % index of neuron within population (1:4999 for L2/3
cell_id = 6; % 6-10 L2/3 PCs, 11-15 L4 LBCs, 16-20 L5 PCs
% note naming convention for weights_file is different for L5 PCs, see 
% plotCNNThreshLayer_ernie_test_L4_Fig.m
weights_file = sprintf('tms_thresh2_N9_3D_cart_cell_E_center_cell%g_train_3dcnn_rs_adr_e2k_seed1',cell_id);
cell_layer = 2; % index of layer (2 = L2/3, 3 = L4, 4 = L5)
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
Efield.interp_method = 'scattered_interp'; % 'scattered_interp', 'mesh_interp', or 'poly'
Efield.Efield_table = 'SimNIBS_TMS_Efield_sims.csv';
Efield.use_scalar_potentials = 0; 
% Sampling grid settings
sample_method_struct = struct(); 
sample_method_struct.method = 'box';
sample_method_struct.Nx = 9;
sample_method_struct.Ny = 9;
sample_method_struct.Nz = 9;
sample_method_struct.lx = 2;
sample_method_struct.ly = 2;
sample_method_struct.lz = 2;
sample_method_struct.rshift = [0 0 0]; 

E_out = interpEfield1(cell_id,nrn_model_ver,Efield,...
                      'sample_method_struct',sample_method_struct);
E_vecs = E_out{1}; 
dnn_output = getDNNOutput({E_vecs},weights_file);