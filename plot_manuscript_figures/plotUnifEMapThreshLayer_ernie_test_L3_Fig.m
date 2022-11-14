%% Plot unifE map threshold on l5 surface
function plotUnifEMapThreshLayer_ernie_test_L3_Fig(save_fig,cbar_mode,save_data)
if nargin == 0
   save_fig = 0; 
   cbar_mode = 1; % 1 horz 2 vert
   save_data = 0; 
end
model_prefix_pre_base = 'tms_dt5'; 
layer_set_num = 1;        
nrn_pop_inds = 1:12;
nrn_pop_names = arrayfun(@(x) sprintf('nrn_pop%g',x),nrn_pop_inds,'UniformOutput',0); 
Efield_name = 'M1_PA_MagVenture_MC_B70_ernie';
opts = struct(); 
opts.mode = 1; 
opts.cell_ids = {[];[];[11:15];[];[]}; 
opts.nrn_model_ver = 'maxH';
opts.plot_func = 'median';
opts.plot_err_func = 'median';
opts.plot_errs = 'per_err'; % per_err or per_err_med
opts.cax_lims = [110 400]; % 
opts.plot_vertices = 0; 
if opts.plot_vertices == 1
    opts.err_cax_lims = [-17 20]; % -30 40 (plot_vertices = 0)
    opts.fig_fold = 'figs';
else    
    if strcmp(opts.plot_errs,'per_err') % median percent error
        opts.err_cax_lims = [-29.53939 38.783412]; % plot_vertices = 0
        opts.fig_fold = 'figs2';
    elseif strcmp(opts.plot_errs,'per_err_med') % percent error of median threshold
        opts.err_cax_lims = [-35.1298 37.7762]; % plot_vertices = 0
        opts.fig_fold = 'figs3';
    end
end
opts.z_lims = []; 
opts.save_fig = save_fig; 
opts.fig_name_suffix = 'L3'; 

[threshEs,threshEsM,plot_errs] = plotUnifEMapThreshLayer(model_prefix_pre_base,layer_set_num,nrn_pop_names,...
                        Efield_name,opts);
layer_num = find(~cellfun(@isempty,opts.cell_ids,'UniformOutput',1));   
if save_data
    data_fold = 'summary_data';
    data_filename = sprintf('L%g_test_unif',layer_num);
    save(fullfile(data_fold,[data_filename '.mat']),'threshEs','threshEsM','plot_errs');
    fprintf('Saved layer data to %s\n',data_filename);
end
if cbar_mode == 1 && save_fig
    plot_horz_cbars_figs(layer_num,opts.cax_lims,opts.err_cax_lims,...
                        opts.save_fig,2,opts.fig_fold); 
elseif cbar_mode == 2 && save_fig
    plot_vert_cbars_figs(layer_num,opts.cax_lims,opts.err_cax_lims,opts.save_fig,2); 
end                     
end