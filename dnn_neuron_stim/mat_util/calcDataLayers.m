function data_layer = calcDataLayers(data,cell_model_names,cell_ids,varargin)
%CALCDATALAYERS Calculate summary statistic on data from all cells in each
%layer at each position
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
in.func = 'median'; 
in.norm_mode = 'none';
in.cutoff_quantile = 0.5; % default cutoff quantile if using norm_mode = 'bin_quantile'
in = sl.in.processVarargin(in,varargin); 
num_layers = length(cell_ids); 
data_layer = cell(num_layers,1); % data combined within each layer
mat_dir = addPaths_dnn_neuron_stim; 
%% apply function to all cells within each layer (e.g. median)
for i = 1:num_layers
    if ~isempty(cell_ids{i})
        cell_model_names_i = cellModelNames(cell_ids{i},'mat_dir',mat_dir); % cell names in layer
        [~,~,data_inds] = intersect(cell_model_names_i,cell_model_names); % get indices of layer cells in threshEs        
        data_layer{i} = apply_func(cell2mat(data(data_inds)),in.func);        
    else
        data_layer{i} = 1e6; % placeholder large number to allow UniformOutput in cellfun below, no cells in this layer
    end    
end
if ~strcmp(in.func,'none')
    fprintf('Computed %s in each layer\n',in.func)
end
%% Normalize
if strcmp(in.norm_mode,'min') % divide by minimum (min = 1)
    global_min = min(cellfun(@(x) min(x(:),'omitnan'),data_layer),'omitnan');
    data_layer = cellfun(@(x) x/global_min,data_layer,'UniformOutput',0);
elseif strcmp(in.norm_mode,'max') % divide by maximum 
    global_max = max(cellfun(@(x) max(x(:),'omitnan'),data_layer),'omitnan'); 
    data_layer = cellfun(@(x) x/global_max,data_layer,'UniformOutput',0);
elseif strcmp(in.norm_mode,'exc') % convert to excitability (1/threshold)
    data_layer = cellfun(@(x) 1./x,data_layer,'UniformOutput',0); % 1/thresh
elseif strcmp(in.norm_mode,'max_exc') % convert to excitability and divide by max (max = 1)
    data_layer = cellfun(@(x) 1./x,data_layer,'UniformOutput',0); % 1/thresh
    global_max = max(cellfun(@(x) max(x(:),'omitnan'),data_layer),'omitnan'); % max excitability (min threshold)
    data_layer = cellfun(@(x) x./global_max,data_layer,'UniformOutput',0); 
elseif strcmp(in.norm_mode,'min_layer') % divide by minimum in each layer (min = 1)
    data_layer = cellfun(@(x) x/min(x(:),'omitnan'),data_layer,'UniformOutput',0); 
elseif strcmp(in.norm_mode,'max_layer') % divide by maximum in each layer (max = 1)
    data_layer = cellfun(@(x) x/max(x(:),'omitnan'),data_layer,'UniformOutput',0); 
elseif strcmp(in.norm_mode,'maxabs_layer') % divide by maximum of abs value in each layer (max(abs) = 1)
    data_layer = cellfun(@(x) x/max(abs(x(:)),'omitnan'),data_layer,'UniformOutput',0);     
elseif strcmp(in.norm_mode,'bin_quantile')
    data_layer = cellfun(@(x) double(x > quantile(x,in.cutoff_quantile)),data_layer,'UniformOutput',0); 
end
if ~strcmp(in.norm_mode,'none')
    fprintf('Applied normalization: %s\n',in.norm_mode)
end
end
