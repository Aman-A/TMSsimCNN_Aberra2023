function MeshROI = loadMeshROI(layers_or_layer_set_num)
%LOADMESHROI Load MeshROI for used to create input layer/layer_set_num
%
%   Inputs
%   ------
%   Optional Inputs
%   ---------------
%   Outputs
%   -------
%   Examples
%   ---------------

% AUTHOR    : Aman Aberra
mat_dir = addPaths_dnn_neuron_stim;
if isnumeric(layers_or_layer_set_num)
   layers = loadLayers(layers_or_layer_set_num);
else
   layers = layers_or_layer_set_num;
end
mesh_name = layers(1).mesh_name;
roi_name = layers(1).roi_name;
mesh_roi_name = [mesh_name '_' roi_name];
mesh_roi_folder = fullfile(mat_dir,'output_data','layer_data',mesh_name,mesh_roi_name);
MeshROI_data = load(fullfile(mesh_roi_folder,[mesh_roi_name '.mat']));
MeshROI = MeshROI_data.MeshROI;