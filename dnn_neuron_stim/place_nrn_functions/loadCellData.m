function cell_data = loadCellData(cell_id,nrn_model_ver,data_files,varargin)
% cellDataLoader loads cell morphology data for corresponding cell_id and
% nrn_model_ver
%
% cellDataLoader(cell_id,nrn_model_ver,varargin)
%
% Optional inputs
% --------------
% 'coordinates' or 'C'
% 'comp_types'
% 'parent_inds'
% 'areas'
% 'branch_orders'
% 'diams'
% 'secnames'
% 'sectypes'
% 'graph'
% Description of data structures
% comp_types : num_comp x 1 vector
%              compartment type of each compartment (default is to load):
%              0 - soma, 1 - axon, 2 - node, 3 - myelin, 4 - unmyelin,
%              5 - basal, 6 - apic used for coloring if vals is empty
% sectypes : num_comp x 1 vector
%           1 = soma (no parent)
%           2 = termination from non-bifurcation (parent is a 1,3, or 6) OR far from bifurcation (parent is 4 or 7, L > length constant) ** contains termination point (1) **
%           3 = intermediate from non-bifurcation (parent is 1,3, or 6)
%           4 = parent side of bifurcation (child is a 5,6, or 7) (parent is 3 OR far from bifurcation (4 or 7)) ** contains parent side bifurcation point (1)**
%           5 = child termination from bifurcation (parent is a 1 or 4 and L < length constant) ** contains termination point (1)**
%           6 = child side of bifurcation (parent is a 1 or 4) **contains child side bifurcation point (0)**
%           7 = parent bifurcation and child of bifurcation (parent is a 1 or 4, child is a 5,6, or 7) ** contains parent side bifurcation point (1) **
%
if nargin < 2
   nrn_model_ver = 'maxH'; % default adult, myelinated human model version
end
if nargin < 3
   data_files =  {'coordinates','comp_types','parent_inds','sectypes'};
end
in.cell_models_file = 'cell_models';
in.name_method_old = 1;
in = sl.in.processVarargin(in,varargin);
% if ~strcmp(nrn_model_ver,'maxH')
%     fprintf('WARNING: Make sure to regenerate cell_data files to ensure soma is centered at origin\n'); 
% end
mat_dir = addPaths_dnn_neuron_stim; % use default mat_dir
cell_data_fold = fullfile(mat_dir,'cell_data',nrn_model_ver);
cell_data = struct();
for i = 1:length(data_files)
    data_type = data_files{i};
    if strcmp(data_type,'C')
        data_type = 'coordinates'; % replace with coordinates
    end
    if strcmp(in.cell_models_file,'cell_models') && in.name_method_old
        cell_filename = sprintf('%s%g',data_type,cell_id);
    else
        cell_filename = sprintf('%s_%s',data_type,...
            cellModelNames(cell_id,'cell_models_file',in.cell_models_file));
    end
    data_file_path = fullfile(cell_data_fold,data_type,[cell_filename '.mat']);
    if ~exist(data_file_path,'file')       
       % generate cell data files
       fprintf('No cell data found for cell_id %g, nrn_model_ver: ''%s'', generating now...\n',...
                cell_id,nrn_model_ver); 
       saveCellData(cell_id,nrn_model_ver,'cell_models_file',in.cell_models_file);
    end
    data = load(data_file_path);
    var_name = fieldnames(data);
    cell_data.(var_name{1}) = data.(var_name{1});
end
if isempty(fieldnames(cell_data))
   error('No data was loaded, check input arguments');
end
end
