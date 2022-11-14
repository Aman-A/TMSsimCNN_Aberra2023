function file_name = getEcellFileName(data_dir,layer_set_num,Efield_name,...
                                            nrn_pop_name,nrn_model_ver,...
                                            cell_layer,cell_id,interp_method,...
                                            use_scalar_potentials,varargin)
%GETECELLFILENAME Output full file path to specific Ecell file
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
in.sample_method = 'all';
in.sample_method_struct = []; % should be input if sample_method ~= 'all'
in.comp_type = 'all';
in.sectype = 'all';
in.reposition_mode = 'off'; % default no repositioning
in = sl.in.processVarargin(in,varargin);
% full path name
Edata_dir = fullfile(data_dir,'output_data','nrn_efields',...
    sprintf('layer_set_%g',layer_set_num),Efield_name,...
    getNrnPopFileName(nrn_pop_name,nrn_model_ver,in.reposition_mode));
% Efield file name depends on interpolation method and E vectors vs. potentials
if strncmp(interp_method,'simnibs_mesh_interp',length('simnibs_mesh_interp'))
    file_name = sprintf('L%g_cell%g_s',cell_layer,cell_id);
elseif strcmp(interp_method,'scattered_interp') % scattered interp
    file_name = sprintf('L%g_cell%g',cell_layer,cell_id);
else
    error('%s interp_method not supported\n',interp_method);
end
if use_scalar_potentials % load scalar potentials (tES only)
    file_name = [file_name 'v'];
end
if ~strcmp(in.sample_method,'all')
    file_name = sprintf('%s_%s',file_name,in.sample_method);
    if strcmp(in.sample_method,'box')
        file_name = sprintf('%s_N%g',file_name,in.sample_method_struct.Nx);
        if in.sample_method_struct.lz ~= 1 % add length of box if not 1 mm
            file_name = sprintf('%s_l%g',file_name,in.sample_method_struct.lz);
        end
        if isfield(in.sample_method_struct,'rshift') % if field exists and not empty or [0 0 0]
           if ~isempty(in.sample_method_struct.rshift) && ~all(in.sample_method_struct.rshift == 0)
               file_name = sprintf('%s_rz%.2f',file_name,in.sample_method_struct.rshift(3)); % add z-component of rshift to name
           end
        end
    elseif strcmp(in.sample_method,'cylinder')
        file_name = sprintf('%s_cyl_Nr%g_Nphi%g_Nz%g_lr%g_lz%g',file_name,in.sample_method_struct.Nr,...
                                                                in.sample_method_struct.Nphi,...
                                                                in.sample_method_struct.Nz,...
                                                                in.sample_method_struct.lr,...
                                                                in.sample_method_struct.lz);
    end
else
    if ~(strcmp(in.comp_type,'all') && strcmp(in.sectype,'all')) % not all for both
        file_name = sprintf('%s_c%s_s%s',...
            file_name,in.comp_type([1:2,regexp(in.comp_type,'_')]),...
            in.sectype(1));
    end
end
file_name = fullfile(Edata_dir,[file_name '.mat']); 
end