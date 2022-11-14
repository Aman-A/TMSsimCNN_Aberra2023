function output = getDNNOutputLayer(Efield,weights_file,varargin)
%GETDNNOUTPUTLAYER Pass E-fields from single layer to python
%keras/tensorflow model and predict neural response output
%
%   Inputs
%   ------
%   Efield : cell array or ND double array
%            E-field tensor at each position in layer
%   weights_file : string
%                  model to load and use for prediction
%   Optional Inputs
%   ---------------
%   Outputs
%   -------
%   Examples
%   ---------------

% AUTHOR    : Aman Aberra
if nargin == 0
   [~,data_dir] = addPaths_dnn_neuron_stim;
   sample_method_struct.Nx = 9; sample_method_struct.lz = 1;
   Ecell = loadEcell(2,6,'nrn_pop1','maxH',1,'M1_PA_MagVenture_MC_B70_ernie',...
                    'simnibs_mesh_interp',0,data_dir,...
                    'sample_method','box','sample_method_struct',sample_method_struct);    
   Efield = processEfieldSamples(Ecell,'3D_cart',sample_method_struct.Nx);
   weights_file = 'tms_thresh2_N9_3D_cart_cell_E_center_cell6_train_3dcnn_rs_adr_e2k_seed1';
end
in.conda_env = 'ml_env37';
if ismac
%     in.python_path = sprintf('/Users/$USER/opt/miniconda3/envs/%s/bin/python',in.conda_env);
    in.python_path = sprintf('/Users/$USER/miniforge3/envs/%s/bin/python',in.conda_env);
else
    in.python_path = 'python'; % load module before calling
end
in.py_func = 'runNetPredict.py';
in = sl.in.processVarargin(in,varargin);
[dnn_dir,dnn_data_dir] = addPaths_dnn_neuron_stim;
weights_dir = fullfile(dnn_data_dir,'cnn_data','weights');
% Convert Efield to ND tensor if cell array
if iscell(Efield)
    % Get invalid indices
    nan_inds = cellfun(@(x) any(isnan(x),'all'),Efield,'UniformOutput',1);
    Efield = cat(5,Efield{:});
    Efield = permute(Efield,[5,1:4]);
else
    % Get invalid indices
    nan_inds = any(isnan(Efield),2:ndims(Efield));
end
% Save to temp file
tmp_efield_file = [tempname '.hdf5'];
save(tmp_efield_file,'Efield','-v7.3');
fprintf('Saved E-field to tmp file %s\n',tmp_efield_file);
% Run in python
weights_file_path = fullfile(weights_dir,[weights_file '.hdf5']);
% Call python function
current_dir = pwd; % go to DNN_neuron_stim/
cd(dnn_dir);
sys_command = sprintf('%s %s %s %s',in.python_path,in.py_func,weights_file_path,tmp_efield_file);
res = system(sys_command);
cd(current_dir); % back to starting dir
if res ~= 0
    delete(tmp_efield_file);
    error('There was an error running %s\n',in.py_func);
end
% Load output
[out_dir,out_name,out_ext] = fileparts(tmp_efield_file);
output_file = fullfile(out_dir,[out_name, '_response',out_ext]);
output_data = h5read(output_file,'/NeuralResponse')'; % transpose output
output = nan(size(Efield,1),size(output_data,2));
output(~nan_inds,:) = output_data; % populate in correct indices (non-nans)

delete(tmp_efield_file);
delete(output_file);
end