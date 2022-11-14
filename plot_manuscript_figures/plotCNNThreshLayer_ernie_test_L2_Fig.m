function plotCNNThreshLayer_ernie_test_L2_Fig(save_fig,cbar_mode,save_data)
if nargin == 0
   save_fig = 0;  
   cbar_mode = 1; % 1 horz 2 vert
   save_data = 1;
end
model_prefix_pre_base = 'tms_dt5'; 
layer_set_num = 1;
layer_num = 2;
cell_ids = 6:10;  
nrn_pop_inds = 1:12;
Efield_name = 'M1_PA_MagVenture_MC_B70_ernie';
% Efield_name = 'M1_PA_MagVenture_MC_B70';
method = 'cell';
Nsamples = 9;
lz = 2; 
rshift = [0 0 -0.4868]; 
E_mode = '3D_cart';    
weights_file_basestring = 'tms_thresh2_N9_3D_cart_cell_E_center_cell%g_train_3dcnn_rs_adr_e2k_seed1';
% weights_file_basestring = 'tms_thresh2_N9_3D_sph_cell_E_center_cell%g_train_3dcnn_rs_adr_e2k_s5';
% weights_file_basestring = 'tms_thresh_N9_3D_sph_cell_E_center_cell%g_3dconvnet1e3_bn0'; 
% weights_file_basestring = 'tms_thresh2_N9_3D_sph_cell_E_center_cell%g_randsearch_tms2g';
opts = struct(); 
opts.mode = 1; 
opts.plot_output_var = 'threshold';
opts.plot_errs = 'per_err'; % 'per_err' 
% opts.plot_errs = 'rsq'; % R^2 at each position across clones/rotations
opts.nrn_pop_basename = 'nrn_pop';
opts.Eopts.drop_phi = 0;
opts.Eopts.theta_mode = 'cos';
opts.Eopts.phi_mode = 'cos';
opts.Eopts.reorder_E = 0;
opts.amp_mode = 'E_center';
opts.rshift = rshift; 
opts.plot_cell_func = 'median';
opts.err_plot_cell_func = 'median';
opts.cax_lims = [110 400];
opts.plot_vertices = 0; 
opts.plot_region = 'all'; % 'G_precentral','S_central','G_postcentral'
if cbar_mode == 1
    if opts.plot_vertices == 1
        opts.err_cax_lims = [-7.57 24.66]; % plot_vertices = 1
        opts.fig_fold = 'figs';
    else
        if strcmp(opts.plot_errs,'per_err') % median percent error 
            opts.err_cax_lims  = [-23.743633, 37.063662]; % plot_vertices = 0 
            opts.fig_fold = 'figs2';
        elseif strcmp(opts.plot_errs,'per_err_med') % percent error of median threshold
            opts.err_cax_lims = [-19.26 45.54];
            opts.fig_fold = 'figs3';
        end
%         opts.err_cax_lims = [0.2057,0.9998]; % R^2
    end
    % opts.err_cax_lims = [-12.92 12.24];  
else
    opts.err_cax_lims = [-23 67]; % matches unifE 
end
opts.z_lims = []; 
opts.ax_view = [-89.2 45]; % default view
% opts.ax_view = [-51.9 15.2]; % ventral view
opts.save_fig = save_fig;

[nrn_data,dnn_data,plot_errs] = plotCNNThreshLayer(model_prefix_pre_base,layer_set_num,layer_num,...
                    cell_ids,nrn_pop_inds,Efield_name,method,Nsamples,...
                    lz,E_mode,weights_file_basestring,opts);
if save_data
    data_fold = 'summary_data';
    data_filename = 'L2_test';
    if ~exist(data_fold,'dir')
        mkdir(data_fold);
    end
    save(fullfile(data_fold,[data_filename '.mat']),'nrn_data','dnn_data','plot_errs');
    fprintf('Saved layer data to %s\n',data_filename);
end
%% Plot colorbars
if cbar_mode == 1 && save_fig
    plot_horz_cbars_figs(layer_num,opts.cax_lims,opts.err_cax_lims,...
                        opts.save_fig,[1 2],opts.fig_fold); 
elseif cbar_mode == 2 && savefig
    plot_vert_cbars_figs(layer_num,opts.cax_lims,opts.err_cax_lims,opts.save_fig,2); 
end
end