function Maps = getMapsSomaE(model_prefix,layer_set_num,Efield_name,nrn_pop_names,...
                         nrn_model_ver,map_type,cell_ids,varargin)
%GETMAPSSOMAE Use uniform E threshold-direction or polarization-direction map
%and Efield distribution to to estimate thresholds or polarizations of
%neurons, respectively
%Uses E at cell body
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
in.save_maps = 1;
in = sl.in.processVarargin(in,varargin);
[mat_dir,data_dir] = addPaths_dnn_neuron_stim;
map_fold = fullfile(mat_dir,'nrn_sim_data','map_data');
if isunix && ~ismac % linux (on cluster)
    data_dir = mat_dir;     
end
layersE = loadLayers(layer_set_num,'opt','layersE','Efield_name',Efield_name);
if iscell(nrn_pop_names)
    num_pops = length(nrn_pop_names);
elseif ischar(nrn_pop_names)
    nrn_pop_names = {nrn_pop_names};
    num_pops = 1;
end
if num_pops == 1
    map_filename = sprintf('%s_ls_%g_E_%s_P_%s',model_prefix,layer_set_num,...
                            Efield_name,nrn_pop_names{1});
else
    map_filename = sprintf('%s_ls_%g_E_%s_P_%s-%s_all',model_prefix,layer_set_num,...
                            Efield_name,nrn_pop_names{1},nrn_pop_names{end});
end
map_filename = fullfile(map_fold,[map_filename '.mat']); % full path
if exist(map_filename,'file')
    fprintf('Loading maps from %s\n',map_filename); 
    Maps_data = load(map_filename); 
    Maps = Maps_data.Maps; 
else
    fprintf('Creating maps...\n'); 
    nrn_phis_all = cell(num_pops,1);
    for pop = 1:num_pops
        NeuronPop = loadNeuronPop(layer_set_num,nrn_pop_names{pop},nrn_model_ver);
        nrn_phis_all{pop} = NeuronPop.phis;
    end
    fprintf('Loaded phis from %g NeuronPops\n',num_pops);
    switch map_type
        case 'polarization'
            deltaVmMap_data = load(fullfile(data_dir,'nrn_sim_data',[model_prefix '.mat']));
            cell_model_names = deltaVmMap_data.cell_model_names;
            deltaVmU = deltaVmMap_data.deltaVms; % uniform E threshold-direction maps
            thetasU = deltaVmMap_data.thetas;
            phisU = deltaVmMap_data.phis;
            deltaVmUsoma  = processNrnCompData(deltaVmU,getCellIds(cell_model_names),...
                nrn_model_ver,'max','soma','all');
            deltaVm_int_grids = makeMWgriddedInterp(thetasU,phisU,deltaVmUsoma);
            num_cells = length(cell_model_names);
            deltaVms = cell(1,num_cells);
            num_layers =  length(cell_ids);
            fprintf('Interpolating thresholds using E at somas\n');
            E_vals = cell(length(num_layers),length(num_pops));
            fprintf('Interpolating polarizations using E at somas\n');
            for li = 1:num_layers
                cells_in_layer = length(cell_ids{li});
                Ei = layersE(li).Efield(:,4:6); % Efield vectors in fem coordinate space
                cell_normals = layersE(li).cell_normals;
                for cj = 1:cells_in_layer
                    cell_model_name = cellModelNames(cell_ids{li}(cj));
                    [~,~,cell_ind] = intersect(cell_model_name,cell_model_names);
                    deltaVms{cell_ind} = zeros(size(Ei,1),num_pops);
                    deltaVm_int_grid_ij = deltaVm_int_grids{cell_ind};
                    % transform E to local neuron space (s-d axis aligned to [0 0 1])
                    for pop = 1:num_pops
                        Ei_loc = zeros(size(Ei));
                        for posk = 1:size(Ei,1)
                            Ei_loc(posk,:) = reorientEfield(cell_normals(posk,:),nrn_phis_all{pop}{li}{cj}(posk),Ei(posk,:));
                            E_vals{li}{cj} = Ei_loc;
                        end
                        [Etheta_fem,Ephi_fem,Emag_fem] = cart2sphD(Ei_loc(:,1),...
                                                                   Ei_loc(:,2),...
                                                                   Ei_loc(:,3),...
                                                                   'phi_mode',2); % sets to -180 to 180 deg                      
                        deltaVms{cell_ind}(:,pop) = Emag_fem.*deltaVm_int_grid_ij(Etheta_fem,Ephi_fem);
                    end
                end
            end
            Maps.model_prefix = model_prefix;
            Maps.layersE = layersE;
            Maps.Ei_loc = E_vals;
            Maps.Efield_name = layersE(1).Efield_name;
            Maps.nrn_pop_names = nrn_pop_names;
            Maps.nrn_phis_all = nrn_phis_all; 
            Maps.cell_ids = cell_ids; 
            Maps.cell_model_names = cell_model_names;
            Maps.deltaVm_int_grids = deltaVm_int_grids;
            Maps.deltaVms = deltaVms;
        case 'threshold'
            threshmap_data = load(fullfile(data_dir,'nrn_sim_data',[model_prefix '.mat']));
            cell_model_names = threshmap_data.cell_model_names;
            threshEsU = threshmap_data.threshEs; % uniform E threshold-direction maps
            thetasU = threshmap_data.thetas;
            phisU = threshmap_data.phis;
            thresh_int_grids = makeMWgriddedInterp(thetasU,phisU,threshEsU);
            num_cells = length(cell_model_names);
            threshEs = cell(1,num_cells);
            num_layers =  length(cell_ids);
            fprintf('Interpolating thresholds using E at somas\n');
            for li = 1:num_layers
                cells_in_layer = length(cell_ids{li});
                if size(layersE(li).Efield,2) == 6 % old format [x,y,z,Ex,Ey,Ez]
                    Ei = layersE(li).Efield(:,4:6); % Efield vectors in fem coordinate space
                else
                    Ei = layersE(li).Efield; % Efield vectors in fem coordinate space 
                end
                cell_normals = layersE(li).cell_normals;
                for cj = 1:cells_in_layer
                    cell_model_name = cellModelNames(cell_ids{li}(cj));
                    [~,~,cell_ind] = intersect(cell_model_name,cell_model_names);
                    threshEs{cell_ind} = zeros(size(Ei,1),num_pops);
                    thresh_int_grid_ij = thresh_int_grids{cell_ind};
                    % transform E to local neuron space (s-d axis aligned to [0 0 1])
                    for pop = 1:num_pops
                        Ei_loc = zeros(size(Ei));
                        for posk = 1:size(Ei,1)
                            Ei_loc(posk,:) = reorientEfield(cell_normals(posk,:),nrn_phis_all{pop}{li}{cj}(posk),Ei(posk,:));
                        end
                        [phis,lambdas,Emag_fem] = cart2sph(Ei_loc(:,1),Ei_loc(:,2),Ei_loc(:,3)); % convert to spherical coords
                        Etheta_fem = 90-lambdas*180/pi; % convert to polar angle (angle from z) in deg
                        Ephi_fem = phis*180/pi; % convert to deg
                        Ephi_fem(Ephi_fem >= 180) = Ephi_fem(Ephi_fem >= 180) - 360; % set to -180 to 180 deg                    
                        threshEs{cell_ind}(:,pop) = (1./Emag_fem).*thresh_int_grid_ij(Etheta_fem,Ephi_fem);
                    end
                end
            end
            Maps.model_prefix = model_prefix;
            Maps.layersE = layersE;
            Maps.Efield_name = layersE(1).Efield_name;
            Maps.nrn_pop_names = nrn_pop_names;
            Maps.nrn_phis_all = nrn_phis_all; 
            Maps.cell_ids = cell_ids; 
            Maps.cell_model_names = cell_model_names;
            Maps.thresh_int_grids = thresh_int_grids;
            Maps.threshEs = threshEs;
    end
    if in.save_maps
        save(map_filename,'Maps');
        fprintf('Saved Maps to %s\n',map_filename);
    end
end