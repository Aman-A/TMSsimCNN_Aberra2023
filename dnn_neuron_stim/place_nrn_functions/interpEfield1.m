function [v_or_E_out] = interpEfield1(cell_id,nrn_model_ver,Efield,varargin)
%INTERPEFIELD1 Interpolates E-field at all neuron compartments of 1 neuron
%in population
%
%   Inputs
%   ------
%   Efield : struct
%            defined in loadParams.m, includes fields:
%

%   Optional Inputs
%   ---------------
%   Outputs
%   -------
%   Examples
%   ---------------

% AUTHOR    : Aman Aberra
mat_dir = addPaths_dnn_neuron_stim;
if nargin == 0
    cell_id = 16;
    nrn_model_ver = 'maxH';
    Efield.E_position = 1;
    Efield.cell_layer = 4;
    Efield.layer_set_num = 1;
    Efield.nrn_pop_name = 'nrn_pop1';
    Efield.E_file = 'M1_PA_MagVenture_MC_B70';
    Efield.interp_method = 'simnibs_mesh_interp'; % 'scattered_interp', 'mesh_interp', or 'poly'
    Efield.Efield_table = 'SimNIBS_TMS_Efield_sims.csv';
    Efield.use_scalar_potentials = 0;     
end
in.phi = 0;
in.out_fill = 'nan'; % sets what to assign coordinates outside of mesh 
                     % for interp_method = 'simnibs_mesh_interp'
in.knn = 10; % number of nearest numbers to use for 
             % interp_method = 'scattered_interp'
in.ROI_expand_factor = 0.5;
in.comp_type = 'all';
in.sectype = 'all';
in.cell_models_file = 'cell_models';
in.data_dir = mat_dir; 
in.sample_method_struct = []; % input to interp E-field on sampling grid instead of cell
in = sl.in.processVarargin(in,varargin);
% Get Efield parameters
layer_set_num = Efield.layer_set_num;
pos_ind = Efield.E_position;
nrn_pop_name = Efield.nrn_pop_name;
Efield_name = Efield.E_file;
cell_layer = Efield.cell_layer;
interp_method = Efield.interp_method;
Efield_table = Efield.Efield_table;
use_scalar_potentials = Efield.use_scalar_potentials;

%  Get cell coordinate file names
input_folder = fullfile(mat_dir,'input_data');
layers = loadLayers(layer_set_num);
NeuronPop = loadNeuronPop(layers,nrn_pop_name,nrn_model_ver);
cell_ids = NeuronPop.cell_ids;
% check that cell is placed in input layer (to match phi)
if ~any(cell_ids{cell_layer} == cell_id)
   fprintf('WARNING: cell_id %g was not placed in layer %g in this NeuronPop',cell_id,cell_layer);
   phi = in.phi;
else
   % get phi from already generated NeuronPop
   phis = NeuronPop.phis;
   phi = phis{cell_layer}{cell_ids{cell_layer} == cell_id}(pos_ind);
end
cell_origin = layers(cell_layer).cell_origins(pos_ind,:);    
cell_normal = layers(cell_layer).cell_normals(pos_ind,:);
if isempty(in.sample_method_struct)
    Cdata = loadCellData(cell_id,nrn_model_ver,...
                                 {'coordinates','comp_types','sectypes'},...
                                  'cell_models_file',in.cell_models_file);
    C = Cdata.C*1e-3; % coordinates of cell (local coordinates) in mm
    if ~(strcmp(in.comp_type,'all') && strcmp(in.sectype,'all')) % not all for both
        inds = getCellCompInds(Cdata,in.comp_type,in.sectype);
        C = C(inds,:);
    end
else
   C = samplePts(in.sample_method_struct,'print_level',1);     	
end
%% Set up folders for loading/saving E data
% Folders for loading E data
Efield_folder = fullfile(input_folder,'fem_efield_data',layers(1).mesh_name);
Efield_file = fullfile(Efield_folder,[Efield_name '.mat']); % full E-field, only used if .msh solution file doesn't exist
Efield_ROI_folder = fullfile(Efield_folder,layers(1).roi_name);
if exist(Efield_ROI_folder,'dir') == 0
   mkdir(Efield_ROI_folder);
end
%% Interpolate E field and save to Efield_pop_fold using either method
switch interp_method
    case 'scattered_interp' % interpolate in MATLAB using ScatteredInterpolant
        fprintf('Interpolating E-field using ScatteredInterpolant method in MATLAB (scattered_interp)\n')
        % Save/load E-field within relevant ROI
        Efield_ROI_file = fullfile(Efield_ROI_folder,[Efield_name '_' layers(1).roi_name sprintf('_%g',in.ROI_expand_factor) '.mat']);
        out = loadEfieldROI(layers,Efield_ROI_file,Efield_file,Efield_table,...
                            interp_method,'ROI_expand_factor',in.ROI_expand_factor,...
                            'data_dir',in.data_dir);
        E = out{1};
        if length(out) == 2
            include_v = 1; % voltages loaded
            v = out{2};
        else
            include_v = 0;
        end
        % Get neuron model coordinates
        Cij = placeCell(cell_origin,cell_normal,C,phi); % get coordinates of all compartments for kth cell (in mm)
        if use_scalar_potentials
            % Interpolate scalar potential at neuron model compartment coordinates
            if include_v
                vpts = v(:,1:3);
                vvals = v(:,4);
                % Interpolate v on compartments
                inds = knnsearch(vpts,Cij,'k',in.knn); % find k nearest points in v
                unique_inds = unique(inds); % extract unique points
                pts_near_Cij = vpts(unique_inds,:); % extract coordinates
                v_near_i = vvals(unique_inds);
                % make scattered interpolants for each component of E
                v_int = scatteredInterpolant(pts_near_Cij(:,1),pts_near_Cij(:,2),pts_near_Cij(:,3),v_near_i,'linear');
                % function output
                v_or_E_out = v_int(Cij(:,1),Cij(:,2),Cij(:,3))*1e3; % potentials in mV
%               fprintf('Outputting v (in mV) on compartments\n');
            else
               error('scalar potential data not available\n');
            end
        else % Interpolate E-field at neuron model compartment coordinates
            Epts = E(:,1:3);
            Evec = E(:,4:6);
            % Interpolate Efield on compartments
            inds = knnsearch(Epts,Cij,'k',in.knn); % find k nearest points in E
            unique_inds = unique(inds); % extract unique points
            pts_near_Cij = Epts(unique_inds,:); % extract coordinates
            E_near_ij = Evec(unique_inds,:);
            % make scattered interpolants for each component of E
            Ex_int = scatteredInterpolant(pts_near_Cij(:,1),pts_near_Cij(:,2),pts_near_Cij(:,3),E_near_ij(:,1),'linear');
            Ey_int = scatteredInterpolant(pts_near_Cij(:,1),pts_near_Cij(:,2),pts_near_Cij(:,3),E_near_ij(:,2),'linear');
            Ez_int = scatteredInterpolant(pts_near_Cij(:,1),pts_near_Cij(:,2),pts_near_Cij(:,3),E_near_ij(:,3),'linear');
            % Get full vector
            Eint = [Ex_int(Cij(:,1),Cij(:,2),Cij(:,3)),...
                Ey_int(Cij(:,1),Cij(:,2),Cij(:,3)),...
                Ez_int(Cij(:,1),Cij(:,2),Cij(:,3))]; % get interpolated field components
            % function output
            v_or_E_out = {Eint,cell_normal,phi};
        end
    case 'simnibs_mesh_interp' % interpolate
        fprintf('Interpolating E-field using SimNIBS FEM interpolation method (simnibs_mesh_interp), cropped mesh\n')
        fprintf('Using out_fill = %s\n',in.out_fill);
        Efield_ROI_file = fullfile(Efield_ROI_folder,[Efield_name '_' layers(1).roi_name sprintf('_%g',in.ROI_expand_factor) '.msh']);
        loadEfieldROI(layers,Efield_ROI_file,Efield_file,Efield_table,...
                      interp_method,'ROI_expand_factor',in.ROI_expand_factor,...
                      'data_dir',in.data_dir);
        % Interpolate E-field at neuron model compartment coordinates
        Cij = placeCell(cell_origin,cell_normal,C,phi); % get coordinates of all compartments for kth cell (in mm)
        Edata = get_fields_at_coordinates(Efield_ROI_file,Cij,in.out_fill);
        Edata_str_arr = [Edata{:}];
        names = {Edata_str_arr.name};

        if use_scalar_potentials
            v_data_ind = strcmp(names,'v');
           if any(v_data_ind) == 1
              ve = Edata_str_arr(v_data_ind).data*1e3; % potentials in mV
              fprintf('V: %g points outside mesh, assigned with %s\n',sum(isnan(ve)),in.out_fill)
              v_or_E_out = ve;
%               fprintf('Outputting v (in mV) on compartments\n');
           else
               fprintf('v does not exist in this solution file\n');
           end
        else
            Eint = Edata_str_arr(strcmp(names,'E')).data;
            fprintf('E: %g points outside mesh, assigned with %s\n',sum(isnan(Eint(:,1))),in.out_fill);
            v_or_E_out = {Eint,cell_normal,phi};
%             fprintf('Outputting E (in V/m) on compartments\n');
        end
    case 'simnibs_mesh_interp2'
        fprintf('Interpolating E-field using SimNIBS FEM interpolation method (simnibs_mesh_interp2), full mesh\n')
        T = readtable(Efield_table);
        Eind = strcmp(T.Efield_name,Efield_name);
        if ~any(Eind)
            err_str = ['%s is not a valid Efield_name in %s, must be:\n',repmat('%s \n',1,length(Eind))];
            error(err_str,Efield_name,Efield_table,T.Efield_name{:});
        end
        mesh_name = T.mesh_name{Eind};
        sol_file = T.sol_file{Eind};
        sim_struct_name = T.sim_struct_name{Eind};
        mesh_sol_file_path = fullfile(mat_dir,'..','simnibs',mesh_name,sim_struct_name,[sol_file '.msh']);
        % Interpolate E-field at neuron model compartment coordinates
        Cij = placeCell(cell_origin,cell_normal,C,phi); % get coordinates of all compartments for kth cell (in mm)
        Edata = get_fields_at_coordinates(mesh_sol_file_path,Cij,in.out_fill);
        Edata_str_arr = [Edata{:}];
        names = {Edata_str_arr.name};
        if use_scalar_potentials
            v_data_ind = strcmp(names,'v');
           if any(v_data_ind) == 1
              v_or_E_out{1} = Edata_str_arr(v_data_ind).data*1e3; % potentials in mV
%               fprintf('Outputting v (in mV) on compartments\n');
           else
               fprintf('v does not exist in this solution file\n');
           end
        else
            Eint = Edata_str_arr(strcmp(names,'E')).data;
            v_or_E_out = {Eint,cell_normal,phi};
        end
    case 'poly_fit'
        fprintf('Using polynomial for local interpolation\n')
        Edata_dir = fullfile(mat_dir,'output_data','nrn_efields',...
            sprintf('layer_set_%g',layer_set_num),Efield_name);
        nrn_pop_name_full = [nrn_pop_name '_' nrn_model_ver];
        if strncmp(Efield.poly_interp_method,'simnibs_mesh_interp',length('simnibs_mesh_interp'))
            fit_file_name = 'fits_s';
        else % scattered interp
            fit_file_name = 'fits';
        end
        if use_scalar_potentials % load scalar potentials (tES only)
            fit_file_name = [fit_file_name 'v'];
        end
        fit_file_name = [fit_file_name '_p' num2str(Efield.poly_ord)];
        fit_data = load(fullfile(Edata_dir,nrn_pop_name_full,[fit_file_name '.mat']));
        fprintf('Loaded coefficients for polynomial of order %g, %s coordinates\n',Efield.poly_ord,fit_data.coords);
        fit_cell_ids = fit_data.cell_ids;
        coeffs_cell = fit_data.coeffs_all{fit_cell_ids == cell_id};
        poly_model.ModelTerms = fit_data.model_terms;
        % Use polynomial coefficients to get E/v at cell compartments in
        % global coords
        if use_scalar_potentials == 0 % E
            Eint = zeros(size(C)); % local E
            for j = 1:3
               poly_model.Coefficients = coeffs_cell(pos_ind,:,j); % get coefficients
               Eint(:,j) = polyvaln(poly_model,C); % E in local coordinates
            end
            if strcmp(fit_data.coords,'spherical') % convert to cartesian coordinates
                [Ex,Ey,Ez] = sph2cart(Eint(:,2)*pi/180,pi*(90-Eint(:,1))/180,Eint(:,3),'phi_mode',fit_data.phi_mode);
                Eint = [Ex,Ey,Ez];
                fprintf('Converted E vectors to cartesian coordinates\n')
            end
            % coefficients fit for local coordinates, need to convert to
            % global
            Eint = placeCell([],cell_normal,Eint,phi); % rotate E vectors to global coordintes (treat as cell coordinates)
            v_or_E_out = {Eint,cell_normal,phi};
        else % v
            poly_model.Coefficients = coeffs_cell(pos_ind,:); % get coefficients
            v_or_E_out = polyvaln(poly_model,C);
        end
    otherwise
        error('%s not a valid interp_method',interp_method);
end
fprintf('Finished E-field interpolation\n'); 
end
