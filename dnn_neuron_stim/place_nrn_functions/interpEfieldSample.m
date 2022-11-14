function interpEfieldSample(layer_set_num,nrn_pop_name,nrn_model_ver,...
                        Efield_name,over_write,interp_method,Efield_table,...
                        sample_method_struct,varargin)
%INTERPEFIELDSAMPLE Interpolates E-field at grid centered on neuron
%locations
%
%   Inputs
%   ------
%   layer_set_num : integer
%                   specify layer set to populate with neurons
%   nrn_pop_name : string
%                  name of neuron population (set of neuron models in
%                  layer and their azimuthal orientations)
%   nrn_model_ver : string
%                   name corresponding to set of neuron model parameters,
%                   specified in outputMorphParams
%   Efield_name : string
%                     name of .mat file with E-field vectors from FEM
%   over_write : 1 or 0
%               set to 1 to overwrite existing data, otherwise skips
%   method : string
%           set to 'mesh_interp' to interpolate within FEM mesh using
%           SimNIBs functions or 'scattered_interp' to interpolate using
%           scatteredInterpolant in MATLAB
%   Efield_table : string
%           csv file to load Efield information from
%   sample_method_struct : struct
%                         parameters for sampling method, fields:
%                         method : 'box', 'sphere', or 'cylinder'
%                         lx,ly,lz : x, y, z length in mm for 'box'
%                         Nx, Ny, Nz : number of points in x, y, z for
%                         'box'
%                         lr : sphere radius in mm for 'sphere'
%                         Nr, Nth, Nph : number of points in r, theta, phi
%                         (spherical coordinates)
%                         lrho,lz : cylinder radius and
%                         height in mm for 'cylinder'
%                         Nrho,Nz : number of points in rho, z
%
%   Optional Inputs
%   ---------------
%   cell_ids : cell array
%              cell ids to override default cells in NeuronPop, must be
%              included in NeuronPop.cell_ids
%   Outputs
%   -------
%   Examples
%   ---------------

% AUTHOR    : Aman Aberra
in.cell_ids = {};
in.ROI_expand_factor = 0.5;
in.out_fill = 'nan';
in.knn = 10;
in = sl.in.processVarargin(in,varargin);
if nargin == 0
   layer_set_num = 1;
   nrn_pop_name = 'nrn_pop1';
   nrn_model_ver = 'maxH';
   Efield_name = 'M1_PA_MagVenture_MC_B70_ernie';
   over_write = 0;
   interp_method = 'simnibs_mesh_interp'; % 'scattered_interp' or 'simnibs_mesh_interp'
   Efield_table = 'SimNIBS_TMS_Efield_sims.csv';
   sample_method_struct.method = 'box';
   sample_method_struct.Nx = 9;
   sample_method_struct.Ny = 9;
   sample_method_struct.Nz = 9;
   sample_method_struct.lx = 1;
   sample_method_struct.ly = 1;
   sample_method_struct.lz = 1;
end
% Get cell coordinate file names
mat_dir = addPaths_dnn_neuron_stim;
output_folder = fullfile(mat_dir,'output_data');
input_folder = fullfile(mat_dir,'input_data');
layers = loadLayers(layer_set_num);
layer_set_name = layers(1).layer_set_name;
% cell coordinates placed in ROI
mesh_name = layers(1).mesh_name;
nrn_pop_name_full = [nrn_pop_name '_' nrn_model_ver];
NeuronPop = loadNeuronPop(layers,nrn_pop_name,nrn_model_ver);
phis = NeuronPop.phis;
cell_ids = NeuronPop.cell_ids;
if ~isempty(in.cell_ids)
   % check that cell ids in in.cell_ids exist in NeuronPop.cell_ids
   check_diff_ids = cellfun(@(x,y) setdiff(x,y),...
                            in.cell_ids,cell_ids,'UniformOutput',0);
   if any(~cellfun(@isempty,check_diff_ids))
      error('input cell_ids should only include cell models that are included in the NeuronPop (%s)',nrn_pop_name_full);
   else
       cell_ids = in.cell_ids;
   end
   cellid_str = num2str([cell_ids{:}],'%g, ');
   fprintf('Running cells %s\n',cellid_str(1:end-1))
end
% MeshROI = loadMeshROI(layer_set_num);
num_layers = length(layers);
num_cells = sum(cellfun(@length,cell_ids));
%% Set up folders for loading/saving E data
% Folders for loading E data
Efield_folder = fullfile(input_folder,'fem_efield_data',mesh_name);
Efield_file = fullfile(Efield_folder,[Efield_name '.mat']); % full E-field, only used if .msh solution file doesn't exist
Efield_ROI_folder = fullfile(Efield_folder,layers(1).roi_name);
% Folders for saving E data
output_dir = fullfile(output_folder,'nrn_efields',layer_set_name);
if exist(output_dir,'dir') == 0
    mkdir(output_dir);
end
nrn_efield_fold = fullfile(output_dir,Efield_name);
if exist(nrn_efield_fold,'dir') == 0
    mkdir(nrn_efield_fold);
end
Efield_pop_fold = fullfile(nrn_efield_fold,nrn_pop_name_full);
if exist(Efield_pop_fold,'dir') == 0
    fprintf('Creating output directory %s\n',Efield_pop_fold);
    mkdir(Efield_pop_fold)
else
    fprintf('Saving to existing output directory %s\n',Efield_pop_fold);
end
% Get sampling coordinates
% overwrites sample_method_struct with defaults for any missing fields
[C,sample_method_struct] = samplePts(sample_method_struct,'print_level',1); 
%% Interpolate E field and save to Efield_pop_fold using either method
switch interp_method
    case 'scattered_interp' % interpolate in MATLAB using ScatteredInterpolant
        fprintf('Using ScatteredInterpolant method in MATLAB\n')
        % Save/load E-field within relevant ROI
        Efield_ROI_file = fullfile(Efield_ROI_folder,...
            [Efield_name '_' layers(1).roi_name sprintf('_%g',in.ROI_expand_factor) '.mat']);
        out = loadEfieldROI(layers,Efield_ROI_file,Efield_file,Efield_table,...
                            interp_method,'ROI_expand_factor',in.ROI_expand_factor);
        E = out{1};
        Epts = E(:,1:3);
        Evec = E(:,4:6);
        if length(out) == 2
            include_v = 1; % voltages loaded
            v = out{2};
            vpts = v(:,1:3);
            vvals = v(:,4);
            fprintf('Potential data included\n')
        else
            include_v = 0;
            vpts = []; 
            vvals = []; 
            fprintf('No potential data\n')
        end
        %% Interpolate E-field at neuron model compartment coordinates
        % Set up parpool
        numCPUs = feature('numCores');
        fprintf('Number of CPUs requested = %g\n',numCPUs);        
        slurm_job_id = getenv('SLURM_JOB_ID'); % check if part of slurm job
        if ~isempty(slurm_job_id)
            pc_storage_dir = fullfile(mat_dir,'pc_storage',getenv('SLURM_JOB_ID'));
            mkdir(pc_storage_dir);
            pc = parcluster('local');
            pc.JobStorageLocation =  pc_storage_dir;
            poolobj = parpool(pc,numCPUs-1);
        else % use default
            pc = parcluster('local');
            poolobj = parpool(pc,numCPUs);
        end
        knn = in.knn;
        for i = 1:num_layers
            num_cells_in_layer = length(cell_ids{i});
            num_positions = layers(i).num_elem;
            % Get cell origins and normals
            cell_originsi = layers(i).cell_origins;
            cell_normalsi = layers(i).cell_normals;
            for j = 1:num_cells_in_layer
                cell_id = cell_ids{i}(j);
                Ecell_fileij = getEcellFileName(mat_dir,layer_set_num,Efield_name,...
                                              nrn_pop_name,nrn_model_ver,i,...
                                              cell_id,interp_method,...
                                              0,'sample_method',sample_method_struct.method,...
                                              'sample_method_struct',sample_method_struct);
                vcell_fileij = getEcellFileName(mat_dir,layer_set_num,Efield_name,...
                                              nrn_pop_name,nrn_model_ver,i,...
                                              cell_id,interp_method,...
                                              1,'sample_method',sample_method_struct.method,...
                                              'sample_method_struct',sample_method_struct);
                if ~exist(Ecell_fileij,'file') || (~exist(vcell_fileij,'file') && include_v) || over_write
                    %  get azimuthal rotations for jth cell in layer
                    phisij = phis{i}{j}; % numPositions x 1 vector of azimuthal rotations
                    Call = cell(num_positions,1); % sampling coordinates in global coord system
                    Ecell_global = cell(num_positions,1); % E vecs in global coord system
                    Ecell = cell(num_positions,1); % E vecs in local coord system
                    if include_v
                        vcell = cell(num_positions,1);
                    end
                    fprintf('Looping through all %g positions...\n',num_positions);
                    tic;
                    for k = 1:num_positions
                        %             fprintf('Position %g\n',j)
                        % Extract coordinates of jth position
%                                     Cij = Call(1+(j-1)*numComp:j*numComp,:); % ith cell, jth position coordinates
                        % get sample points for kth cell position of jth
                        % cell in ith layer                        
                        Cijk = placeCell(cell_originsi(k,:),cell_normalsi(k,:),C,phisij(k)); 
%                         Cij = placeCell(cell_originsi(k,:),cell_normalsi(k,:),C,phisij(k)); % get coordinates of all compartments for kth cell (in mm)
                        inds = knnsearch(Epts,Cijk,'k',knn); % find 10 nearest points in E
                        unique_inds = unique(inds); % extract unique points
                        pts_near_Cij = Epts(unique_inds,:); % extract coordinates
                        E_near_i = Evec(unique_inds,:);
                        % make scattered interpolants for each component of E
                        Ex_int = scatteredInterpolant(pts_near_Cij(:,1),pts_near_Cij(:,2),pts_near_Cij(:,3),E_near_i(:,1),'linear');
                        Ey_int = scatteredInterpolant(pts_near_Cij(:,1),pts_near_Cij(:,2),pts_near_Cij(:,3),E_near_i(:,2),'linear');
                        Ez_int = scatteredInterpolant(pts_near_Cij(:,1),pts_near_Cij(:,2),pts_near_Cij(:,3),E_near_i(:,3),'linear');

                        Eint = [Ex_int(Cijk(:,1),Cijk(:,2),Cijk(:,3)),...
                            Ey_int(Cijk(:,1),Cijk(:,2),Cijk(:,3)),...
                            Ez_int(Cijk(:,1),Cijk(:,2),Cijk(:,3))]; % get interpolated field components
                        if size(Eint) ~= size(Cijk)
                            error('E-field has different number of elements from cell-coordinates');
                        end
                        % Add to full array of Efield vectors (global)
                        Ecell_global{k} = Eint;
                        % rotate E-field to standard cell orientation
                        Eint = reorientEfield(cell_normalsi(k,:),phisij(k),Eint);
                        % Add to full array of Efield vectors (local)
                        Ecell{k} = Eint;
                        Call{k} = Cijk;
                        if include_v
                            % Interpolate v on compartments
                            inds = knnsearch(vpts,Cijk,'k',knn); % find k nearest points in v
                            unique_inds = unique(inds); % extract unique points
                            pts_near_Cijk = vpts(unique_inds,:); % extract coordinates
                            v_near_i = vvals(unique_inds);
                            % make scattered interpolants for each component of E
                            v_int = scatteredInterpolant(pts_near_Cijk(:,1),pts_near_Cijk(:,2),pts_near_Cijk(:,3),v_near_i,'linear');
                            % function output
                            vcell{k} = v_int(Cijk(:,1),Cijk(:,2),Cijk(:,3))*1e3; % potentials in mV
                            %               fprintf('Outputting v (in mV) on compartments\n');
                        else
%                             fprintf('WARNING: scalar potential data not availabile\n')
                        end
                        %             hold on;
                        %             scale_factor = 0.1;
                        %             plot3(pts(j,1),pts(j,2),pts(j,3),'Marker','o','MarkerSize',40,'Color','g')
                        %             quiver3(pts_near_Cij(:,1),pts_near_Cij(:,2),pts_near_Cij(:,3),E_near_i(:,1)*scale_factor,E_near_i(:,2)*scale_factor,E_near_i(:,3)*scale_factor,'k','Autoscale','off')
                        %             quiver3(pts(j,1),pts(j,2),pts(j,3),Eint(j,1)*scale_factor,Eint(j,2)*scale_factor,Eint(j,3)*scale_factor,'r','Autoscale','off');
                        %             drawnow;
                    end
                    fprintf('Done\n');
                    toc                    
                    % Save cell array in .mat file
                    Ecell_dataij = struct(); 
                    Ecell_dataij.Ecell = Ecell; 
                    Ecell_dataij.Ecell_global = Ecell_global; 
                    Ecell_dataij.Call = Call; 
                    Ecell_dataij.sample_method_struct = sample_method_struct; 
                    Ecell_dataij.C = C;
                    save(Ecell_fileij,'-STRUCT','Ecell_dataij');
                    if include_v
                        vcell_dataij = struct(); 
                        vcell_dataij.vcell = vcell; 
                        vcell_dataij.Call = Call;
                        vcell_dataij.sample_method_struct;
                        vcell_dataij.C = C; 
                        save(vcell_fileij,'-STRUCT','vcell_dataij');
                        fprintf('L%g of %g: Saved E and v for cell %g of %g (%g total cells)\n',i,num_layers,j,num_cells_in_layer,num_cells);
                    else
                        fprintf('L%g of %g: Saved E for cell %g of %g (%g total cells)\n',i,num_layers,j,num_cells_in_layer,num_cells);
                    end
                else
                    fprintf('file exists for cell %g, over_write = %g\n',cell_id,over_write);
                end
            end
        end
        delete(poolobj)
        if exist(pc_storage_dir,'dir')
            rmdir(pc_storage_dir,'s'); % delete parfor temp files if necessary
        end
   case 'simnibs_mesh_interp' % interpolate
        fprintf('Using FEM interpolation method with SimNIBs functions, cropped mesh\n')
        fprintf('Using out_fill = %s\n',in.out_fill);
        Efield_ROI_file = fullfile(Efield_ROI_folder,[Efield_name '_' layers(1).roi_name sprintf('_%g',in.ROI_expand_factor) '.msh']);
        out = loadEfieldROI(layers,Efield_ROI_file,Efield_file,Efield_table,interp_method);
        if isnan(out)
            m = mesh_load_gmsh4(Efield_ROI_file);
            if isempty(m.node_data)
                include_v = 0;
                fprintf('No potentials loaded\n')
            else
                include_v = 1;
                fprintf('Includes potentials\n')
            end
        elseif out == 1 || out == 0
            include_v = out;
        end
        numCPUs = feature('numCores');
        fprintf('Number of CPUs requested = %g\n',numCPUs);
        slurm_job_id = getenv('SLURM_JOB_ID'); 
        if ~isempty(slurm_job_id)
            pc_storage_dir = fullfile(mat_dir,'pc_storage',slurm_job_id);
            mkdir(pc_storage_dir);
            pc = parcluster('local');
            pc.JobStorageLocation =  pc_storage_dir;
            poolobj = parpool(pc,numCPUs-1);
        else % use default
            pc = parcluster('local');
            poolobj = parpool(pc,numCPUs);
        end
        out_fill = in.out_fill;
        fprintf('using out_fill mode = %s\n',out_fill);
        for i = 1:num_layers
            num_cells_in_layer = length(cell_ids{i});
            num_positions = layers(i).num_elem;
            % Get cell origins and normals
            cell_originsi = layers(i).cell_origins;
            cell_normalsi = layers(i).cell_normals;
            for j = 1:num_cells_in_layer
                cell_id = cell_ids{i}(j);
                Ecell_fileij = getEcellFileName(mat_dir,layer_set_num,Efield_name,...
                                              nrn_pop_name,nrn_model_ver,i,...
                                              cell_id,interp_method,...
                                              0,'sample_method',sample_method_struct.method,...
                                              'sample_method_struct',sample_method_struct);
                vcell_fileij = getEcellFileName(mat_dir,layer_set_num,Efield_name,...
                                              nrn_pop_name,nrn_model_ver,i,...
                                              cell_id,interp_method,...
                                              1,'sample_method',sample_method_struct.method,...
                                              'sample_method_struct',sample_method_struct);                  
                if ~exist(Ecell_fileij,'file') || (~exist(vcell_fileij,'file') && include_v) || over_write
                    %  get azimuthal rotations for jth cell in layer
                    phisij = phis{i}{j}; % numPositions x 1 vector of azimuthal rotations
                    Call = cell(num_positions,1); % sampling coordinates in global coord system
                    Ecell_global = cell(num_positions,1); % E vecs in global coord system
                    Ecell = cell(num_positions,1); % E vecs in local coord system
                    if include_v
                        vcell = cell(num_positions,1);
                    end
                    fprintf('Looping through all %g positions...\n',num_positions);
                    tic;
                    parfor k = 1:num_positions
                        %             fprintf('Position %g\n',j)
                        % Extract coordinates of jth position
                        %             Cij = Call(1+(j-1)*numComp:j*numComp,:); % ith cell, jth position coordinates
                        Cijk = placeCell(cell_originsi(k,:),cell_normalsi(k,:),C,phisij(k));
%                         Cijk = samplePts(cell_originsi(k,:),cell_normalsi(k,:),phisij(k),sample_method_struct);
%                         Edata = get_fields_at_coordinates(m,Cij,out_fill);
%                         fprintf('position %g...\n',k)
                        Edata = get_fields_at_coordinates(Efield_ROI_file,Cijk,out_fill);
                        Edata_str_arr = [Edata{:}];
                        names = {Edata_str_arr.name};
                        Eint = Edata_str_arr(strcmp(names,'E')).data;
                        % Add to full array of Efield vectors (global
                        % coordinate system)
                        Ecell_global{k} = Eint;
                        % rotate E-field to standard cell orientation
                        Eint = reorientEfield(cell_normalsi(k,:),phisij(k),Eint);
                        % Add to full array of Efield vectors (local
                        % coordinate system)
                        Ecell{k} = Eint;
                        Call{k} = Cijk;
                        if include_v
                           vcell{k} = Edata_str_arr(strcmp(names,'v')).data*1e3; % convert from V to mV
                        else
%                             fprintf('WARNING: scalar potential data not availabile\n')
                        end
                    end
                    fprintf('Done\n');
                    toc
                    % check nans
                    if strcmp(in.out_fill,'nan')
                        fprintf('Number of nan (positions with undefined comp values): %g\n',...
                              sum(cellfun(@(x) any(isnan(x(:,1))),Ecell,'UniformOutput',1)));
                    end
                    % Save cell arrays in .mat files                    
                    Ecell_dataij = struct(); 
                    Ecell_dataij.Ecell = Ecell; 
                    Ecell_dataij.Ecell_global = Ecell_global; 
                    Ecell_dataij.Call = Call; 
                    Ecell_dataij.sample_method_struct = sample_method_struct; 
                    Ecell_dataij.C = C;
                    save(Ecell_fileij,'-STRUCT','Ecell_dataij');
                    if include_v
                        vcell_dataij = struct(); 
                        vcell_dataij.vcell = vcell; 
                        vcell_dataij.Call = Call;
                        vcell_dataij.sample_method_struct;
                        vcell_dataij.C = C; 
                        save(vcell_fileij,'-STRUCT','vcell_dataij');
                        fprintf('L%g of %g: Saved E and v for cell %g of %g (%g total cells)\n',i,num_layers,j,num_cells_in_layer,num_cells);
                    else
                        fprintf('L%g of %g: Saved E for cell %g of %g (%g total cells)\n',i,num_layers,j,num_cells_in_layer,num_cells);
                    end
                else
                    fprintf('file exists for cell %g, over_write = %g\n',cell_id,over_write);
                end
            end
        end
        delete(poolobj)
        if exist(pc_storage_dir,'dir')
            rmdir(pc_storage_dir,'s'); % delete parfor temp files if necessary
        end
    otherwise
        error('%s not a valid interp_method',interp_method);
end
end

