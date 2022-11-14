function cell_model_names = cellModelNames(cell_ids,varargin)
% CELLMODELNAMES opens table with names of neuron models
%input scalar or vector of cell_ids, outputs string or cell array,
%respectively, of cell_model_names
in.mat_dir = '';
in.cell_models_file = 'cell_models'; % 25 cells in Aberra 2018
% in.cell_models_file = 'cell_models_all'; % all 1035 BB cells
in = sl.in.processVarargin(in,varargin);
if isempty(in.mat_dir)
   in.mat_dir = addPaths_dnn_neuron_stim;
end
cell_models_data_file = fullfile(in.mat_dir,'cell_data',[in.cell_models_file '.mat']);
cell_models_data = load(cell_models_data_file);
T = cell_models_data.T;
if nargin == 0 || isempty(cell_ids)
   cell_ids = 1:size(T,1); % output all cell ids
end
[~,inds] = ismember(cell_ids,T.cell_id);
if ~isempty(inds) && ~any(inds == 0)
    cell_model_names = T.Properties.RowNames(inds);
%     cell_model_names = T.Row(inds);
    if length(cell_model_names) == 1
        cell_model_names = cell_model_names{1};
    end
else
   error('All cell_ids should be between %g and %g',T.cell_id(1),T.cell_id(end))
end
end