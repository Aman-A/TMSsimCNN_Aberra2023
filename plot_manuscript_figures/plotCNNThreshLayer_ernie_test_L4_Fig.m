function plotCNNThreshLayer_ernie_test_L4_Fig(save_fig,cbar_mode,save_data)
if nargin == 0
   save_fig = 1;  
   cbar_mode = 1; % 1 horz 2 vert
   save_data = 1;
end
model_prefix_pre_base = 'tms_dt5'; 
layer_set_num = 1;
layer_num = 4;
cell_ids = 16:20;  
nrn_pop_inds = 1:12;
Efield_name = 'M1_PA_MagVenture_MC_B70_ernie';
method = 'cell';
Nsamples = 9;
lz = 1.5; 
E_mode = '3D_cart';    
weights_file_basestring = sprintf('tms_thresh2_l%g_N9_3D_cart_cell_E_center_cell%%g_train_3dcnn_rs_adr_e2k_seed1',lz);
% weights_file_basestring = 'tms_thresh_N9_3D_sph_cell_E_center_cell%g_3dconvnet1e3_bn0'; 
% weights_file_basestring = 'tms_thresh2_N9_3D_sph_cell_E_center_cell%g_randsearch_tms2g';
opts = struct(); 
opts.mode = 1; 
opts.plot_output_var = 'threshold';
opts.plot_errs = 'per_err'; % 'per_err' or 'per_err_med'
opts.nrn_pop_basename = 'nrn_pop';
opts.Eopts.drop_phi = 0;
opts.Eopts.theta_mode = 'cos';
opts.Eopts.phi_mode = 'cos';
opts.amp_mode = 'E_center';
opts.plot_cell_func = 'median';
opts.err_plot_cell_func = 'median';
opts.cax_lims = [90 400];
opts.plot_vertices = 0; 
if cbar_mode == 1 % horz sep colorbars for cnn/unif
%     opts.err_cax_lims = [-0.67 2.15]; 
    if opts.plot_vertices == 1
        opts.err_cax_lims = [-0.82 1.66]; % l=1.5
        opts.fig_fold = 'figs';
    else        
        if strcmp(opts.plot_errs,'per_err') % median percent error 
            opts.err_cax_lims = [-1.945464 3.511911]; 
            opts.fig_fold = 'figs2';
        elseif strcmp(opts.plot_errs,'per_err_med') % percent error of median threshold
            opts.err_cax_lims = [-4.57 5.13]; 
            opts.fig_fold = 'figs3';
        end
    end
else % vert same colorbar for both
    opts.err_cax_lims = [-18 26]; % match unifE 
    opts.fig_fold = 'figs';
end
% opts.err_cax_lims = [-20 100];  % match unifE for cell16/nrn_pop1
% opts.err_cax_lims = [-12.92 12.24];  
opts.z_lims = []; 
opts.save_fig = save_fig;
[nrn_data,dnn_data,plot_errs] = plotCNNThreshLayer(model_prefix_pre_base,layer_set_num,layer_num,...
                    cell_ids,nrn_pop_inds,Efield_name,method,Nsamples,...
                    lz,E_mode,weights_file_basestring,opts);
if save_data
    data_fold = 'summary_data';
    data_filename = 'L4_test';
    save(fullfile(data_fold,[data_filename '.mat']),'nrn_data','dnn_data','plot_errs');
    fprintf('Saved layer data to %s\n',data_filename);
end
%% Plot colorbars
if cbar_mode == 1 && save_fig
    plot_horz_cbars_figs(layer_num,opts.cax_lims,opts.err_cax_lims,...
                        opts.save_fig,[1 2],opts.fig_fold); 
elseif cbar_mode == 2 && save_fig
    plot_vert_cbars_figs(layer_num,opts.cax_lims,opts.err_cax_lims,opts.save_fig); 
end
end