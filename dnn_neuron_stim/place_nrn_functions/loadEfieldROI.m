function out = loadEfieldROI(layers,Efield_ROI_file,Efield_file,...
                             Efield_table,interp_method,varargin)
%LOADEFIELDROI ...
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
in.ROI_expand_factor = 0.5;
in.load_msh = 0; % default dont load/output .msh since it will be reloaded anyway
in.data_dir = mat_dir;
in = sl.in.processVarargin(in,varargin);
[~,Efield_name] = fileparts(Efield_file);
Efield_ROI_folder = fileparts(Efield_ROI_file);
switch interp_method
    case 'scattered_interp'
        if exist(Efield_ROI_file,'file') == 0
            if exist(Efield_file,'file')
                % get Efield from .mat saved in Efield_folder
                fprintf('Loading full Efield from %s\n',Efield_file);
                Edata = load(Efield_file);
                E_all = Edata.E;
                if isfield(Edata,'E_wm')
                    E_all = [E_all;Edata.E_wm];
                    fprintf('Incorporating E-field in WM\n');
                end
                if isfield(Edata,'v')
                    v = Edata.v;
                    include_v = 1;
                else
                    include_v = 0;
                end
                % potentials not saved in full .mat file, use .msh approach
                % below
            else
                % get Efield from .msh and extract ROI E-field values                
                T = readtable(Efield_table);
                Eind = strcmp(T.Efield_name,Efield_name);
                if ~any(Eind)
                    err_str = ['%s is not a valid Efield_name in %s, must be:\n',repmat('%s \n',1,length(Eind))];
                    error(err_str,Efield_name,Efield_table,T.Efield_name{:});
                end
                sol_file = T.sol_file{Eind};
                sim_struct_name = T.sim_struct_name{Eind};
                mesh_sol_file_path = fullfile(in.data_dir,'..','simnibs',layers(1).mesh_name,sim_struct_name,[sol_file '.msh']);
                fprintf('Loading E-field solution file from %s\n',mesh_sol_file_path); 
                m = mesh_load_gmsh4(mesh_sol_file_path);
                tet_centers = mesh_get_tetrahedron_centers(m);
                gm_ind = 2; wm_ind = 1;
                field_idx = get_field_idx(m,'E','element');
                % Get all E points and vectors at tet centers
                E_all = [tet_centers(m.tetrahedron_regions==gm_ind | m.tetrahedron_regions==wm_ind,:),...
                    m.element_data{field_idx}.tetdata(m.tetrahedron_regions==gm_ind | m.tetrahedron_regions==wm_ind,:)];
                if length(m.node_data) == 1
                    nodes = m.nodes;
                    v = [nodes, m.node_data{1}.data];
                    include_v = 1;
                else
                    include_v = 0;
                end
            end
            MeshROI = loadMeshROI(layers);
            ROI = MeshROI.ROI;
            for i = [1,3,5] % expand ROI to include points just outside
                ROI(i:i+1) = ROI(i:i+1)+in.ROI_expand_factor*diff(ROI(i:i+1))*[-1,1];
            end
            % Use ROI to filter E-field solution to just points of interest
            E = getEROI(E_all,ROI); % [ x, y, z, Ex, Ey, Ez]
            if exist(Efield_ROI_folder,'dir') == 0
                mkdir(Efield_ROI_folder);
                fprintf('Made directory %s\n',Efield_ROI_folder);
            else
                fprintf('Saving to %s\n',Efield_ROI_folder);
            end
            if include_v
                v = getEROI(v,ROI);
                save(Efield_ROI_file,'E','v');
                fprintf('Saved %s Efield (E) and potentials (v) in ROI\n',Efield_ROI_file)
            else
                save(Efield_ROI_file,'E');
                fprintf('Saved %s Efield in ROI\n',Efield_ROI_file)
            end
        else
            fprintf('Loading already created E-ROI file %s\n',Efield_ROI_file)
            Edata = load(Efield_ROI_file);
            E = Edata.E;
            if isfield(Edata,'v')
                v = Edata.v;
                include_v = 1;
            else
                include_v = 0;
            end
        end
        out = {E};
        if include_v
           out = [out,v];
        end
    case {'simnibs_mesh_interp'}
        if exist(Efield_ROI_file,'file') == 0
            T = readtable(Efield_table);
            Eind = strcmp(T.Efield_name,Efield_name);
            if ~any(Eind)
                err_str = ['%s is not a valid Efield_name in %s, must be:\n',repmat('%s \n',1,length(Eind))];
                error(err_str,Efield_name,Efield_table,T.Efield_name{:});
            end
            mesh_name = T.mesh_name{Eind};
            sol_file = T.sol_file{Eind};
            sim_struct_name = T.sim_struct_name{Eind};
            mesh_sol_file_path = fullfile(in.data_dir,'..','simnibs',mesh_name,sim_struct_name,[sol_file '.msh']);
            m_full = mesh_load_gmsh4(mesh_sol_file_path);
            MeshROI = loadMeshROI(layers);
            ROI = MeshROI.ROI;
            for i = [1,3,5] % expand ROI to include tets just outside
                ROI(i:i+1) = ROI(i:i+1)+in.ROI_expand_factor*diff(ROI(i:i+1))*[-1,1];
            end
            [m,include_v] = clipMeshROI_gmsh(m_full,ROI,1);
            mesh_save_gmsh4(m,Efield_ROI_file);
            if in.load_msh
                out = {m,include_v};
            else
                out = include_v;
            end
        else
            fprintf('Loading already created E-ROI msh file %s\n',Efield_ROI_file)
            if in.load_msh
                m = mesh_load_gmsh4(Efield_ROI_file);
                if isempty(m.node_data)
                    include_v = 0;
                else
                    include_v = 1;
                end
                out = {m,include_v};
            else
                out = nan;
            end
        end
end
