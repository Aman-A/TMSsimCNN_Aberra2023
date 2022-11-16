function [dnn_data,scale_factors] = dnnEstLayer(cell_ids,weights_files,nrn_pop_names,...
                                                E_mode,layer_num,layer_set_num,Efield_name,...
                                                Eopts,nrn_model_ver,mode,interp_method,...
                                                sample_method_struct,output_var,varargin)
%DNNESTLAYER Run DNN estimation for multiple cells and populatiosn
%(rotations) in single layer. Loads if exists/saves after generating. 
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
[~,dnn_data_dir] = addPaths_dnn_neuron_stim;
in.use_scalar_potentials = 0;
in.amp_mode = 'E_center';
in.stim_amp = 100; % default 100 A/us
in.est_dir = fullfile(dnn_data_dir,'cnn_data','est_data');
in.sim_data_dir = dnn_data_dir;
in = sl.in.processVarargin(in,varargin); 
%%
num_weights = length(weights_files);
if ~iscell(nrn_pop_names)
    nrn_pop_names = {nrn_pop_names};
end
num_pops = length(nrn_pop_names);
if strcmp(output_var,'threshold')
    pre = 'thresh';
elseif strcmp(output_var,'probability')
    pre = 'prob';
elseif strcmp(output_var,'polarization')
    pre = 'pol';
else
    error('output_var = ''%s'' not valid\n',output_var); 
end
dnn_est_data_filename = fullfile(in.est_dir,sprintf('%s_w%g_%s_ls%g_%s',...
                                                    pre,mode,Efield_name,layer_set_num,weights_files{1}));
if num_pops > 1
   dnn_est_data_filename = sprintf('%s_%s-%s',dnn_est_data_filename,nrn_pop_names{1},...
                                                nrn_pop_names{end});
else
    dnn_est_data_filename = sprintf('%s_%s',dnn_est_data_filename,nrn_pop_names{1});
end
if num_weights > 1
    dnn_est_data_filename = sprintf('%s_cell%g-%g',dnn_est_data_filename,...
                                                cell_ids(1),cell_ids(end));
else
    dnn_est_data_filename = sprintf('%s_cell%g',dnn_est_data_filename,...
                                    cell_ids(1));
end
dnn_est_data_filename = [dnn_est_data_filename '.mat'];
if exist(dnn_est_data_filename,'file')
    dnn_saved_data = load(dnn_est_data_filename);
    dnn_data = dnn_saved_data.dnn_data;
    if strcmp(output_var,'threshold')
        scale_factors = dnn_saved_data.scale_factors;
    end
    assert(any(strcmp(weights_files,dnn_saved_data.weights_files)),...
            'input weights_files do not match saved weights_files');
    fprintf('Loaded saved DNN prediction data from %s\n',dnn_est_data_filename);
else
    layers = loadLayers(layer_set_num); 
    num_pos = layers(layer_num).num_elem;     
    fprintf('Running estimation...\n');
    dnn_data = cell(1,num_weights);
    scale_factors = cell(1,num_weights);
    if strcmp(in.amp_mode,'E_center')
       [C,sample_method_struct] = samplePts(sample_method_struct,'print_level',1);
       [~,center_ind,~] = intersect(C-sample_method_struct.rshift,[0 0 0],'rows'); % get middle index for normalizing by |E| at this point 
    end
    for j = 1:num_weights
        weights_filej = weights_files{j};
        dnn_dataj = zeros(num_pos,num_pops);
        for i = 1:num_pops
            nrn_popi = nrn_pop_names{i};
            fprintf('Loading Efield for %s\n',nrn_popi);
            if mode < 0
               reverse_E = 1;
            else
               reverse_E = 0;
            end
            [Ecellij,file_exists,sample_method_struct] = loadEcell(layer_num,cell_ids(j),nrn_popi,nrn_model_ver,...
                              layer_set_num,Efield_name,interp_method,in.use_scalar_potentials,...
                              in.sim_data_dir,'sample_method',sample_method_struct.method,...
                              'sample_method_struct',sample_method_struct,'reverse_E',reverse_E);
            if ~file_exists
               error('Ecell not loaded')
            end
            if strcmp(output_var,'probability')
                Ecellij = cellfun(@(x) x*in.stim_amp,Ecellij,'UniformOutput',0);
            elseif strcmp(output_var,'threshold')
                if strcmp(in.amp_mode,'E_center')
                    % Ecell always in cartesian coordinates                    
                    scale_factorsj = cellfun(@(x) norm(x(center_ind,:)),Ecellij,'UniformOutput',1); % mag of center E-field
                    Ecellij = cellfun(@(x,y) x/y,Ecellij,num2cell(scale_factorsj),'UniformOutput',0); 
                end
            end
            Efieldij = processEfieldSamples(Ecellij,E_mode,sample_method_struct.Nx,Eopts);
            %% Get estimate from DNN
            if strcmp(output_var,'probability')
                firing_probs_dnn = getDNNOutput(Efieldij,weights_filej);
                dnn_dataj(:,i) = firing_probs_dnn;
            elseif strcmp(output_var,'threshold')
                threshEs_dnn = getDNNOutput(Efieldij,weights_filej);
                threshEs_dnn = threshEs_dnn./scale_factorsj; % convert back to dI/dt threshold by dividing by |E|_center in V/m per A/�s
                dnn_dataj(:,i) = threshEs_dnn;
                if i == 1
                   scale_factors{j} = scale_factorsj; % same within layer
                end
            end
        end
        dnn_data{j} = dnn_dataj;
    end
    save_data = struct();
    save_data.dnn_data = dnn_data;
    save_data.output_var = output_var;
    save_data.stim_amp = in.stim_amp;
    save_data.weights_files = weights_files;        
    if strcmp(output_var,'threshold')
        save_data.scale_factors = scale_factors;
    end
    save(dnn_est_data_filename,'-STRUCT','save_data');
    fprintf('Saved estimation data to %s\n',dnn_est_data_filename); 
end