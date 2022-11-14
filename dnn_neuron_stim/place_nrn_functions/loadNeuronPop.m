function NeuronPop = loadNeuronPop(layers_or_layer_set_num,nrn_pop_name,nrn_model_ver,...
                                    varargin)
if nargin == 0
    layers_or_layer_set_num = 1;
    nrn_pop_name = 'nrn_pop1';
    nrn_model_ver = 'maxH';
end
in.reposition_mode = 'off';
in = sl.in.processVarargin(in,varargin); 
mat_dir = addPaths_dnn_neuron_stim;
if isnumeric(layers_or_layer_set_num)
   layers = loadLayers(layers_or_layer_set_num);
else
   layers = layers_or_layer_set_num;
end
layer_set_name = layers(1).layer_set_name;
mesh_name = layers(1).mesh_name;
roi_name = layers(1).roi_name;
mesh_roi_name = [mesh_name '_' roi_name];
mesh_roi_folder = fullfile(mat_dir,'output_data','layer_data',mesh_name,mesh_roi_name);
layer_folder = fullfile(mesh_roi_folder,layer_set_name);
nrn_model_ver_folder = fullfile(layer_folder,nrn_model_ver);
nrn_pop_file = getNrnPopFileName(nrn_pop_name,nrn_model_ver,in.reposition_mode); 
NeuronPop_data = load(fullfile(nrn_model_ver_folder,[nrn_pop_file '.mat']));
fprintf('loading from %s \n',fullfile(nrn_model_ver_folder,[nrn_pop_file '.mat']));
NeuronPop = NeuronPop_data.NeuronPop;