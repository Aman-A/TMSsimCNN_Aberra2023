main_dir = addPaths_dnn_neuron_stim;
download_test_dataset = 0; % set to 1 to download test dataset, consists of:
                           % SimNIBS ernie example mesh meshed with
                           % mri2mesh
                           % ernie_m2m_handknob.mat ->
                           % dnn_neuron_stim/output_data/layer_data/ernie_m2m/ernie_m2m_handknob
                           % layer mesh data (layer_set_1)
                           % 
% After installing python dependencies (tensorflow, keras, numpy, and scipy)
% specify path to python executable by setting python_exec_str in 
% dnn_neuron_stim/mat_util/python_exec.m
path_to_python_executable = python_exec; 
fprintf('Path to python executable successfully added to python_exec: %s\n',...
        path_to_python_executable)
% OPTIONAL: Download test dataset to regenerate Fig3 or run estimation on
% neurons positioned within test E-field distribution (using ernie SimNIBS
% example head mesh, re-meshed with mri2mesh pipeline)

if download_test_dataset
    datafile = 'Aberra2022_TMSsimCNN_dataset';
    target_dir = 'data'; % path to directory test dataset will be 
                         % downloaded to, requires >X GB storage space
    mat_dir = addPaths_dnn_neuron_stim;    
    % Download

    % Unzip
    unzip(fullfile(target_dir,[datafile '.zip']),target_dir);
    % test dataset parameters
    layer_set_num = 1;
    mesh_name = 'ernie_m2m';
    roi_name = 'handknob';
    nrn_pop_inds = 1:12; 
    nrn_model_ver = 'maxH';
    nrn_pop_names = arrayfun(@(x) sprintf('nrn_pop%g_%s',x,nrn_model_ver),...
                            nrn_pop_inds,'UniformOutput',0);
    Efield_name = 'M1_PA_MagVenture_MC_B70_ernie';
    % move layer data 
    copyfile(fullfile(target_dir,datafile,'layer_data/'),...
             fullfile(mat_dir,'output_data','layer_data/'));    
    % move neuron simulation data
    copyfile(fullfile(target_dir,datafile,'nrn_sim_data/'),...
             fullfile(mat_dir,'nrn_sim_data/')); 
    % move cell data
    copyfile(fullfile(target_dir,datafile,'cell_data/'),...
             fullfile(mat_dir,'cell_data/')); 
    % move estimation data (generated from trained models)
    copyfile(fullfile(target_dir,datafile,'est_data/'),...
             fullfile(mat_dir,'cnn_data/est_data/')); 
    % move weights files (from trained models)
    copyfile(fullfile(target_dir,datafile,'weights/'),...
             fullfile(mat_dir,'cnn_data/weights/'));     
    % move simnibs data
    copyfile(fullfile(target_dir,datafile,'simnibs/'),...
             fullfile(mat_dir,'..','simnibs/'));
end