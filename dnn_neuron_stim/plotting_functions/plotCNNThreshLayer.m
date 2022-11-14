function [nrn_data,dnn_data,plot_errs] = plotCNNThreshLayer(model_prefix_pre_base,layer_set_num,layer_num,...
                            cell_ids,nrn_pop_inds,Efield_name,method,Nsamples,...
                            lz,E_mode,weights_files,varargin)
if nargin == 0        
    model_prefix_pre_base = 'tms'; 
    layer_set_num = 1;
    layer_num = 4;
    cell_ids = 16;  
    nrn_pop_inds = 1:6;
    Efield_name = 'M1_PA_MagVenture_MC_B70';
    method = 'cell';
    Nsamples = 9;
    lz = 2; 
    E_mode = '3D_cart';    
end
if nargin < 11
    weights_files = 'tms_thresh_N9_3D_sph_cell_E_center_cell%g_3dconvnet1e3'; 
%     weights_files = arrayfun(@(x) sprintf('tms_thresh_N9_3D_sph_cell_E_center_cell%g_3dconvnet1e3',...
%                         x),cell_ids,'UniformOutput',0); 
end
if ischar(weights_files) && any(~cellfun(@isempty,regexp(weights_files,{'%g','%d'})))
   % assume weights_files is string to format
   weights_file_basestring = weights_files; 
   weights_files = arrayfun(@(x) sprintf(weights_file_basestring,...
                        x),cell_ids,'UniformOutput',0); 
    fprintf('Converted weights_file_basestring %s to %g weights_files\n',...
            weights_file_basestring,length(cell_ids)); 
end
dnn_dir = addPaths_dnn_neuron_stim;
sim_data_dir = dnn_dir; 
in.mode = 1; % waveform mode
in.plot_output_var = 'threshold';
in.nrn_pop_basename = 'nrn_pop';
in.Eopts.drop_phi = 0;
in.Eopts.theta_mode = 'cos';
in.Eopts.phi_mode = 'cos';
in.amp_mode = 'E_center';
in.cax_lims = [80 400];
in.err_cax_lims = [-50 15]; 
in.rshift = [0,0,0]; % default no spatial shift of sampling grid
in.z_lims = [22 52.4057];
in.ax_view = [-89.2 45];
in.plot_vertices = 1;
in.plot_cell_func = []; % gets set inside plotActivationLayer if empty
in.err_plot_cell_func = []; % gets set inside plotActivationLayer if empty
in.plot_errs = 'per_err'; % percent error
in.plot_region = 'all';
in.sim_data_dir = sim_data_dir;
in.save_fig = 0;
in.fig_fold = fullfile(dnn_dir,'figs');
in = sl.in.processVarargin(in,varargin); 
Eopts = Eopts_handler(in.Eopts); 
%% plotActivationLayer options
output_var = 'threshold';
nrn_pop_names = arrayfun(@(x) sprintf('%s%g',in.nrn_pop_basename,x),...
                            nrn_pop_inds,'UniformOutput',0); 
opts = struct; 
opts.mode = in.mode; 
opts.method = method; 
opts.output_var = output_var; 
opts.plot_output_var = in.plot_output_var; 
opts.cax_lims = in.cax_lims;
opts.amp_mode = in.amp_mode;
opts.save_fig = in.save_fig; 
opts.plot_errs = in.plot_errs; 
opts.dnn_nrn_pop_names = nrn_pop_names; 
opts.sim_data_dir = in.sim_data_dir; 
opts.err_cax_lims = in.err_cax_lims; 
opts.fig_fold = in.fig_fold;
opts.ax_view = in.ax_view; 
opts.z_lims = in.z_lims;
opts.rshift = in.rshift; 
opts.plot_cell_func = in.plot_cell_func; 
opts.err_plot_cell_func = in.err_plot_cell_func; 
opts.plot_vertices = in.plot_vertices;
opts.plot_region = in.plot_region; 
%% Plot DNN estimation maps and error
[nrn_data,dnn_data,plot_errs] = plotActivationLayer(model_prefix_pre_base,layer_set_num,...
                    layer_num,cell_ids,nrn_pop_names,Efield_name,...
                    Nsamples,lz,E_mode,Eopts,weights_files,opts);