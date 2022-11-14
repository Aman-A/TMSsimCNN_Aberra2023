function [data_file_name, model_prefix_pre_pop] = getNrnSimDataFileName(...
                                                    model_prefix_pre_name,varargin)
if nargin == 0
    model_prefix_pre_name = 'tdcs_r4_1.8'; 
end
in.model_prefix_pre = []; 
in.model_prefix_params = []; % defaults defined below
in.nrn_pop_names = {'nrn_pop1','nrn_pop12'}; 
in.func = []; % 
in.comp_type = [];
in.sectype = [];
in.func_method = []; 
in.name_style = 2; % 2 - new, 1 - old
in = sl.in.processVarargin(in,varargin);
% create model_prefix_pre (contains sim params up to P_<nrn_pop_name> part of model_prefix)
if ~isempty(in.model_prefix_pre)
   model_prefix_pre = in.model_prefix_pre; % input full name 
else
   p = get_model_prefix_params(in.model_prefix_params); % checks fields and assigns defaults for missing fields
   if p.load_potentials == 1 % TMS/tES in head model
       if strcmp(p.sim_type,'polarization') && strcmp(p.waveform_type,'rect')
           if isempty(p.dur) || isempty(p.amp) || isempty(p.Efield_name)
              error('Make sure to specify dur, amp, and Efield_name');  
           end
           if in.name_style == 1
               model_prefix_pre = sprintf('%s_%gms_%.1fmA_ls_%g_E_%s',p.nrn_model_ver,...
                   p.dur,p.amp,p.layer_set_num,p.Efield_name);               
           elseif in.name_style == 2
               model_prefix_pre = sprintf('%s_%gms_ls_%g_E_%s',p.nrn_model_ver,...
                   p.dur,p.layer_set_num,p.Efield_name);
           end
       elseif strcmp(p.sim_type,'threshold') && strcmp(p.waveform_type,'tms')
           % ignore name_style for now
            model_prefix_pre = sprintf('%s_w%g_ls_%g_E_%s',p.nrn_model_ver,...
                                        p.mode,p.layer_set_num,p.Efield_name); 
       end
   end
end
% make part of model_prefix before nrn_pop name (P_)
model_prefix_pre_pop = [model_prefix_pre_name '_' model_prefix_pre]; 
% Get nrn_pop_names suffix
nrn_pop_names = in.nrn_pop_names; 
if ischar(nrn_pop_names)
   nrn_pop_names = {nrn_pop_names}; % but in cell array format 
end
if length(nrn_pop_names) > 1
    nrn_pop_str = [nrn_pop_names{1} '-' nrn_pop_names{end}];
    comb_pops = 1;
else
    nrn_pop_str = nrn_pop_names{1}; 
    comb_pops = 0;
end
% Processed data file vs. raw data final suffix
if ~isempty(in.func) &&  ~isempty(in.comp_type) && ~isempty(in.sectype) && ~isempty(in.func_method)
    if isnumeric(in.func)
        func_str = num2str(in.func,'q%g'); % quantile function
    else
        func_str = in.func;
    end
    file_suff = sprintf('_%s_%s_c%s_s%s',in.func_method,func_str,...
        in.comp_type([1:2,regexp(in.comp_type,'_')]),...
        in.sectype(1)); % abbreviate comp_type and secytpe
else
    if comb_pops
        file_suff = '_all';
    else
        file_suff = ''; 
    end
end

data_file_name = [model_prefix_pre_pop '_P_' nrn_pop_str file_suff];

end

function in = get_model_prefix_params(varargin)
% default fields/params
in.layer_set_num = 3;
in.nrn_model_ver = 'maxH';
in.mode = 1; % only relevant for waveform_type = 'tms'
in.dur = []; % user needs to specify
in.amp = []; % user needs to specify 
in.Efield_name  = '';
in.sim_type = 'polarization';
in.load_potentials = 1;   
in.waveform_type = 'rect'; 
in = sl.in.processVarargin(in,varargin);

end