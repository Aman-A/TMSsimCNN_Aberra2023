function data = mapEstLayers(Maps,layersE,cell_ids,nrn_pop_names,varargin)
%MAPEST Estimate unif E approx for threshold/polarization data for
%arbitrary layersE, cell_ids, and nrn_pop_names input
%  
%   Inputs 
%   ------ 
%   Ein : input E-field, can be different formats
%   cell_ids : cell array
%              num_layers x 1 array of vectors where each element
%              represents a layer, containing a vector of cell indices
%              contained in that layer
%   Optional Inputs 
%   --------------- 
%   Outputs 
%   ------- 
%   Examples 
%   --------------- 
% TODO: Add option to input phis directly rather than having to use nrn_pop
% phis saved in Maps structure
% AUTHOR    : Aman Aberra 
in.reverse_E = 0; % set to 1 to flip direction E-field vectors, e.g. tms waveform mode = -1 or tES current polarity is reversed
in.scale_E = 1; % scale E-field amp, e.g. for different tDCS stim current 
in = sl.in.processVarargin(in,varargin); 
if in.scale_E ~= 1
   fprintf('Scaling layers E-field by %g\n',in.scale_E);  
end
if in.reverse_E
   fprintf('Reversing layers E-field direction\n');  
end
% Parameters of maps
if isfield(Maps,'threshEs')
    map_type = 'threshold';
elseif isfield(Maps,'deltaVms')
    map_type = 'polarization';
end       
if ~iscell(nrn_pop_names)
    nrn_pop_names = {nrn_pop_names};
end
num_layers = length(cell_ids);        
num_pops = length(nrn_pop_names); 
if isfield(Maps,'nrn_phis_all')
   nrn_phis_all = Maps.nrn_phis_all;    
else
    fprintf('Loading phis and saving to Maps...\n'); 
    % Load phis from NeuronPops
    nrn_phis_all = cell(num_pops,1);
    layer_set_num = Maps.layersE(1).layer_set_num;
    if isfield(Maps,'nrn_model_ver')
       nrn_model_ver = Maps.nrn_model_ver; 
    else
        nrn_model_ver = 'maxH'; % for Maps created before adding nrn_model_ver to struct, assume 'maxH' 
    end
    % check phi step from NeuronPops
    for pop = 1:num_pops
        NeuronPop = loadNeuronPop(layer_set_num,nrn_pop_names{pop},nrn_model_ver);
        nrn_phis_all{pop} = NeuronPop.phis;
    end    
    Maps.nrn_phis_all = nrn_phis_all; 
    [~,data_dir] = addPaths_simTBSnrn;    
    % Resave
    num_map_pops = length(Maps.nrn_pop_names);
    if num_map_pops == 1
        pop_str = Maps.nrn_pop_names{1};
    else
        pop_str = [Maps.nrn_pop_names{1} '-' Maps.nrn_pop_names{end}];
    end    
    if isfield(Maps,'deltaVms_weighted')
        method = 'cell'; % cell method has weighted field, stat doesn't
    else
        method = 'stat';
    end
    if strcmp(map_type,'polarization')
        if isnumeric(Maps.func)
            func_str = num2str(Maps.func,'%.3f_quant');
        else
            func_str = Maps.func;
        end
        map_filename = sprintf('%s_ls_%g_E_%s_P_%s_%s_%s_%s_%s',model_prefix,layer_set_num,...
                Efield_name,pop_str,method,...
                func_str,comp_type,sectype);        
    else 
        map_filename = sprintf('%s_ls_%g_E_%s_P_%s',model_prefix,layer_set_num,...
                Efield_name,pop_str);       
    end    
    map_fold = fullfile(data_dir,'nrn_sim_data','map_data');
    save(fullfile(map_fold,[map_filename '.mat']),'Maps'); 
    fprintf('Resaved Maps to %s\n',map_filename); 
end
switch map_type
    case 'threshold'
        fprintf('Interpolating thresholds using E at somas\n');
        thresh_int_grids = Maps.thresh_int_grids; 
        threshEs = cell(1,length([cell_ids{:}]));
        cell_cnt = 1; % initialize to 1
        
        for li = 1:num_layers
            cells_in_layer = length(cell_ids{li});
            if size(layersE(li).Efield,2) == 6 % old convention
               Ei = layersE(li).Efield(:,4:6); % Efield vectors in fem coordinate space
            else
               Ei = layersE(li).Efield; % Efield vectors in fem coordinate space 
            end
            if in.reverse_E
                Ei = -1*Ei;                
            end
            cell_normals = layersE(li).cell_normals;
            for cj = 1:cells_in_layer
                cell_model_name = cellModelNames(cell_ids{li}(cj));
                [~,~,cell_ind] = intersect(cell_model_name,Maps.cell_model_names);
                threshEs{cell_cnt} = zeros(size(Ei,1),num_pops);
                % transform E to local neuron space (s-d axis aligned to [0 0 1])
                thresh_int_grid_ij = thresh_int_grids{cell_ind};
                for pop = 1:num_pops
                    [~,~,pop_ind] = intersect(nrn_pop_names{pop},Maps.nrn_pop_names);
                    Ei_loc = zeros(size(Ei));
                    for posk = 1:size(Ei,1)
                        Ei_loc(posk,:) = reorientEfield(cell_normals(posk,:),...
                                                        nrn_phis_all{pop_ind}{li}{cj}(posk),Ei(posk,:));
                    end
                    [Etheta_fem,Ephi_fem,Emag_fem] = cart2sphD(Ei_loc(:,1),...
                                                               Ei_loc(:,2),...
                                                               Ei_loc(:,3),...
                                                               'phi_mode',2); % sets to -180 to 180 deg                                        
                    threshEs{cell_cnt}(:,pop) = (1./(Emag_fem*in.scale_E)).*thresh_int_grid_ij(Etheta_fem,Ephi_fem);
                end
                cell_cnt = cell_cnt + 1; 
            end
        end
        data = threshEs; 
        
    case 'polarization'
        fprintf('Interpolating polarizations using E at somas. func: ''%s'', comp_type: ''%s'', sectype: ''%s''\n',...
                func_str,Maps.comp_type,Maps.sectype);
        deltaVm_int_grids = Maps.deltaVm_int_grids; 
        deltaVms = cell(1,length([cell_ids{:}]));
        cell_cnt = 1; % initialize to 1        
        for li = 1:num_layers
            cells_in_layer = length(cell_ids{li});
            if size(layersE(li).Efield,2) == 6 % old convention
               Ei = layersE(li).Efield(:,4:6); % Efield vectors in fem coordinate space
            else
               Ei = layersE(li).Efield; % Efield vectors in fem coordinate space 
            end
            if in.reverse_E
                Ei = -1*Ei;                
            end
            cell_normals = layersE(li).cell_normals;
            for cj = 1:cells_in_layer
                cell_model_name = cellModelNames(cell_ids{li}(cj));
                [~,~,cell_ind] = intersect(cell_model_name,Maps.cell_model_names);
                deltaVms{cell_cnt} = zeros(size(Ei,1),num_pops);
                % transform E to local neuron space (s-d axis aligned to [0 0 1])
                deltaVm_int_grid_ij = deltaVm_int_grids{cell_ind};
                for pop = 1:num_pops
                    [~,~,pop_ind] = intersect(nrn_pop_names{pop},Maps.nrn_pop_names);
                    Ei_loc = zeros(size(Ei));
                    for posk = 1:size(Ei,1)
                        Ei_loc(posk,:) = reorientEfield(cell_normals(posk,:),...
                                                        nrn_phis_all{pop_ind}{li}{cj}(posk),Ei(posk,:));
                    end
                    [Etheta_fem,Ephi_fem,Emag_fem] = cart2sphD(Ei_loc(:,1),...
                                                               Ei_loc(:,2),...
                                                               Ei_loc(:,3),...
                                                               'phi_mode',2); % sets to -180 to 180 deg                                        
                    deltaVms{cell_cnt}(:,pop) = (Emag_fem*in.scale_E).*deltaVm_int_grid_ij(Etheta_fem,Ephi_fem);
                end
                cell_cnt = cell_cnt + 1; 
            end
        end
        data = deltaVms; 
end