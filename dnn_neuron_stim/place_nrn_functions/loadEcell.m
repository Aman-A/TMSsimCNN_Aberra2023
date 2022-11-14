function [Eout_all,file_exists,varargout] = loadEcell(cell_layer,cell_id,...
                                            nrn_pop_name,nrn_model_ver,layer_set_num,Efield_name,...
                                            interp_method,use_scalar_potentials,...
                                            data_dir,varargin)

in.comp_type = 'all';
in.sectype = 'all';
in.load_Ecell = 1;
in.sample_method = 'all'; % if loading interpEfieldSample (sample_method_struct.method)
in.sample_method_struct = []; 
in.reverse_E = []; % flip E-field vectors, i.e. 180� rotation or polarity reversal
in.reposition_mode = 'off'; % default no repositioning
in = sl.in.processVarargin(in,varargin);
% Get file name
if isempty(data_dir)
   data_dir = addPaths_dnn_neuron_stim;
end
if strcmp(Efield_name(end-1:end),'_r') % reversed E-field dist
    Efield_name = Efield_name(1:end-2);    
    in.reverse_E = 1;
elseif isempty(in.reverse_E) % unset, assume no reversal
    in.reverse_E = 0; 
end
settings.sample_method = in.sample_method;
settings.sample_method_struct = in.sample_method_struct;
settings.comp_type = in.comp_type; 
settings.sectype = in.sectype; 
settings.reposition_mode = in.reposition_mode; 
Ecell_file = getEcellFileName(data_dir,layer_set_num,Efield_name,nrn_pop_name,...
                              nrn_model_ver,cell_layer,cell_id,interp_method,...
                              use_scalar_potentials,settings); 
if exist(Ecell_file,'file')
    file_exists = 1;
else
    file_exists = 0;
    fprintf('WARNING: %s does not exist\n',Ecell_file);
end
if in.load_Ecell && file_exists
    cell_data = load(Ecell_file);
    if use_scalar_potentials
        fprintf('Loaded scalar potential distribution from %s\n',Ecell_file);
        vcell = cell_data.vcell;
        Eout_all = vcell;
    else
        fprintf('Loaded E-field distribution from %s\n',Ecell_file);
        Ecell = cell_data.Ecell;
        % Extract E-field vectors for this position
        if in.reverse_E
            if strcmp(in.sample_method,'all')
                % flip E-field vectors but not cell_normal (1st row)
                % TODO: REMOVE CELL_NORMAL FROM ECELL, KEEP IN SEPARATE
                % VARIABLE
                Ecell = cellfun(@(x) [1;-1*ones(length(x)-1,1)].*x,...
                                Ecell,'UniformOutput',0); 
            else
                Ecell = cellfun(@(x) -1*x,Ecell,'UniformOutput',0); 
            end
            fprintf('Flipping E-field data from base distribution: %s\n',Efield_name);
        end
        Eout_all = Ecell;
    end
    if nargout == 3
       if isfield(cell_data,'sample_method_struct')
            varargout = {cell_data.sample_method_struct};  
       else
           error('sample_method_struct not included in Ecell file'); 
       end
    end
else
    Eout_all = [];
end