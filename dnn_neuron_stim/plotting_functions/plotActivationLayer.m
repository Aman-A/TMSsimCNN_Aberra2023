function [nrn_data,dnn_data,plot_errs] = plotActivationLayer(model_prefix_pre_base,layer_set_num,...
                            layer_num,cell_ids,nrn_pop_names,Efield_name,...
                            Nsamples,lz,E_mode,Eopts,weights_files,varargin)
if nargin == 0
    model_prefix_pre_base = 'tms';
    layer_set_num = 1;
    layer_num = 4;
    cell_ids = 16;
    nrn_pop_names = arrayfun(@(x) sprintf('nrn_pop%g',x),1,'UniformOutput',0);
    Efield_name = 'M1_PA_MagVenture_MC_B70_ernie';
    Nsamples = 9;
    lz = 1.5;
    E_mode = '3D_cart';
    Eopts = Eopts_handler(); % assign defaults
    weights_files = {'tms_thresh2_N9_3D_cart_cell_E_center_cell16_train_3dcnn_rs_adr_e2k_seed1'};    
end
[dnn_dir,dnn_data_dir] = addPaths_dnn_neuron_stim;
sim_data_dir = dnn_data_dir; 
in.stim_amp = 100; % A/us
in.output_var = 'probability';
in.plot_output_var = [];
in.plot_errs = [];
in.per_err_cutoff = 5; % how many positions below this percent error (if plot_errs is per_err)
in.method = 'stat';
in.nrn_model_ver = 'maxH';
in.dnn_nrn_pop_names = 'nrn_pop1';
in.mode = 1; % waveform mode
in.rshift = [0,0,0]; 
in.interp_method = 'scattered_interp';
in.waveform_type = 'tms';
in.use_scalar_potentials = 0;
in.sample_method = 'box';
in.n_sample_pts = 20;
in.amp_mode = 'E_center';
in.sim_data_dir = sim_data_dir;
in.est_dir = fullfile(dnn_data_dir,'cnn_data','est_data');
in.plot_vertices = 0;
in.plot_region = 'all';
in.plot_cell_func = []; % set below if still empty
in.err_plot_cell_func = []; % set below if still empty
% in.cax_lims = [0.1 0.6];
in.cax_lims = [0 1];
in.err_cax_lims = [];
in.z_lims = [22 52.4057];
in.ax_view = [-89.2 45];
in.cmap = flipud(fake_parula(1000));
% in.cmap = [0.8*ones(1,3);fake_parula(1000)];
in.plot_actual_v_pred = 0;
in.fig_fold = fullfile(dnn_dir,'figs');
in.save_fig = 1;
in = sl.in.processVarargin(in,varargin);
% More argument handling
Eopts = Eopts_handler(Eopts);
if isempty(in.plot_output_var)
    in.plot_output_var = in.output_var;
    fprintf('Plotting default for this dataset: %s\n',in.plot_output_var)
else
    if ~strcmp(in.output_var,in.plot_output_var)
        fprintf('Plotting %s for DNN that outputs: %s\n',in.plot_output_var,in.output_var)
    end
end
if isempty(in.plot_errs)
    if strcmp(in.plot_output_var,'probability') % default err to plot is based on default set above
        in.plot_errs = 'err';
    elseif strcmp(in.plot_output_var,'threshold')
        in.plot_errs = 'per_err';
    end
end
if ischar(weights_files)
   weights_files = {weights_files}; % single weights_file
end
cell_ids_all = cell_ids;
if strcmp(in.method,'cell')
   assert(length(weights_files) == length(cell_ids),...
        'For cell method, length(weights_files) should equal length(cell_ids)');
end
if strcmp(in.method,'stat') || strcmp(in.plot_output_var,'probability')
    cell_ids = cell_ids(1);
    fprintf('Using first cell_id input for stat method: %g\n',cell_ids);
end
%% Load ground truth data
model_params.sim_type = 'threshold';
model_params.waveform_type = in.waveform_type;
model_params.nrn_model_ver = in.nrn_model_ver;
model_params.mode = in.mode;
model_params.layer_set_num = layer_set_num;
model_params.Efield_name = Efield_name;
model_prefix = getNrnSimDataFileName(model_prefix_pre_base,...
                                    'nrn_pop_names',nrn_pop_names,...
                                    'model_prefix_params',model_params);

thresh_data = load(fullfile(in.sim_data_dir,'nrn_sim_data',[model_prefix '.mat']));
fprintf('Loaded ground truth data: %s\n',model_prefix);
[~,cell_inds_all] = intersect(thresh_data.cell_model_names,cellModelNames(cell_ids_all));
% [~,cell_inds] = intersect(thresh_data.cell_model_names,cellModelNames(cell_ids));
threshEs = thresh_data.threshEs(cell_inds_all); % double array of thresh in A/�s
if strcmp(in.plot_output_var,'probability')
    % convert to firing prob at each position
    fp_nrn_sim_data_file = fullfile(in.est_dir,sprintf('fp%g_%s.mat',in.stim_amp,model_prefix));
    if exist(fp_nrn_sim_data_file,'file')
        fp_data = load(fp_nrn_sim_data_file);
        nrn_data = fp_data.nrn_data;
        fprintf('Loaded firing probs of nrn sim data from %s\n',fp_nrn_sim_data_file);
    else
        fprintf('Calculating firing probs for %g cells on nrn sim data\n',length(threshEs));
        nrn_data = {thresh_to_firing_prob_kernel(cell2mat(threshEs),'stim_amp',in.stim_amp)};
        save(fp_nrn_sim_data_file,'nrn_data');
        fprintf('Saved firing probs of nrn sim data to %s\n',fp_nrn_sim_data_file);
    end
%     nrn_data = firing_probs;
    if isempty(in.plot_cell_func)
        in.plot_cell_func = 'none';
    end    
elseif strcmp(in.output_var,'threshold')
    nrn_data = threshEs;
    if isempty(in.plot_cell_func)
        in.plot_cell_func = 'median';
    end
end
if isempty(in.err_plot_cell_func)
    in.err_plot_cell_func = in.plot_cell_func;
end
% Load layer for plotting
layers = loadLayers(layer_set_num);

% Figure file names
if strcmp(in.plot_output_var,'probability')
    fig_name1 = sprintf('surfFP_%s_l%g_%s_%.1famp',Efield_name,layer_num,model_prefix,in.stim_amp);
    fig_name2 = sprintf('surfFP_%s_%s_%.1famp',Efield_name,weights_files{1},in.stim_amp);
    fig_name3 = sprintf('surfFP_err_%s_%s_%.1famp',Efield_name,weights_files{1},in.stim_amp);
    fig_name4 = sprintf('actual_v_pred_dnn_FP_%s_%s_%.1famp',Efield_name,weights_files{1},in.stim_amp);
elseif strcmp(in.plot_output_var,'threshold')
    fig_name1 = sprintf('surfTh_%s_w%g_l%g_%s',Efield_name,in.mode,layer_num,model_prefix);
    fig_name2 = sprintf('surfTh_%s_w%g_%s',Efield_name,in.mode,weights_files{1});
    if strcmp(in.plot_errs,'per_err')
        fig_name3 = sprintf('surfTh_pererr_%s_w%g_%s',Efield_name,in.mode,weights_files{1});        
    elseif strcmp(in.plot_errs,'rsq')
        fig_name3 = sprintf('surfTh_rsq_%s_w%g_%s',Efield_name,in.mode,weights_files{1});        
    else
        fig_name3 = sprintf('surfTh_%s_%s_w%g_%s',in.plot_errs,Efield_name,in.mode,weights_files{1});
    end
    if strcmp(in.amp_mode,'E_center')
        fig_name4 = sprintf('actual_vs_pred_dnn_Eth_%s_w%g_%s',Efield_name,in.mode,weights_files{1});
    else
        fig_name4 = sprintf('actual_vs_pred_dnn_dIdt_%s_w%g_%s',Efield_name,in.mode,weights_files{1});
    end
end

%% Plot simulated data
cell_ids_cell = cell(1,length(layers)); cell_ids_cell{layer_num} = cell_ids;
plot_cell_model_names = cellModelNames(cell_ids);
args_all = struct(); % args for all plots
args_all.cell_ids = cell_ids_cell;
args_all.cell_func = in.plot_cell_func;
args_all.clims = in.cax_lims;
args_all.cmap = in.cmap;
args_all.z_lims = in.z_lims;
args_all.ax_view = in.ax_view; 
args_all.fig_fold = in.fig_fold;
args_all.cell_model_names = plot_cell_model_names;
args_all.cbar_on = 0;
args_all.plot_vertices = in.plot_vertices;
args_all.save_fig = in.save_fig;
args_all.plot_region = in.plot_region;
figure;
ax = gca;
args1 = args_all;
args1.ax = ax;
args1.fig_name = fig_name1;
args1.threshEs = nrn_data;
plotThreshLayersFig(layer_set_num,[],args1);
fprintf('Plotted NEURON simulation data\n');
%% Load Ecell for running DNN inference
fprintf('Estimating response at %g positions\n',length(threshEs{1}));
sample_method_struct = struct('Nx',Nsamples,'Ny',Nsamples,'Nz',Nsamples,...
                               'lx',lz,'ly',lz,'lz',lz,'rshift',in.rshift,...
                               'method',in.sample_method); 
[dnn_data, scale_factors] = dnnEstLayer(cell_ids,weights_files,nrn_pop_names,...
                                        E_mode,layer_num,layer_set_num,Efield_name,...
                                        Eopts,in.nrn_model_ver,in.mode,in.interp_method,...
                                        sample_method_struct,in.output_var,...
                                        'est_dir',in.est_dir,...
                                        'sim_data_dir',in.sim_data_dir,...
                                        'stim_amp',in.stim_amp); 
%% Convert DNN thresholds to probability
if strcmp(in.output_var,'threshold') && strcmp(in.plot_output_var,'probability')
    threshEs_dnn = cell2mat(dnn_data);
    fprintf('Converting DNN estimated thresholds to firing probabilites for %g positions, %g clones+rotations\n',...
            size(threshEs_dnn,1),size(threshEs_dnn,2));
    dnn_data = {thresh_to_firing_prob_kernel(threshEs_dnn,'stim_amp',in.stim_amp)};
end
%% Plot
figure;
ax2 = gca;
args2 = args_all;
args2.ax = ax2;
args2.fig_name = fig_name2;
args2.threshEs = dnn_data;
plotThreshLayersFig(layer_set_num,[],args2);
%% Get errs
dnn_errs = cellfun(@(x,y) calc_errs(x,y,'all'),nrn_data,dnn_data,...
                    'UniformOutput',0); % predicted - actual
fprintf('Error of DNN at estimating %s across %g clones and %g rotations\n',...
        in.plot_output_var,length(cell_ids_all),length(nrn_pop_names));
plot_errs = cell(1,length(dnn_errs));
for i = 1:length(dnn_errs)
    err_namesi = fieldnames(dnn_errs{i});
    for j = 1:length(err_namesi)
        if ~strcmp(err_namesi{j},'errs')
            fprintf('%s: %f\n',err_namesi{j},mean(dnn_errs{i}.(err_namesi{j})));
        end
    end
    if strcmp(in.plot_errs,'err')
        plot_errsi = dnn_errs{i}.errs;
    elseif strcmp(in.plot_errs,'per_err')
        plot_errsi = 100*dnn_errs{i}.errs./nrn_data{i}; % percent error    
    elseif strcmp(in.plot_errs,'rsq')        
        plot_errsi = zeros(size(dnn_data{i},1),1);
        for k = 1:size(dnn_data{i},1)
            [Rk,~] = corrcoef(dnn_data{i}(k,:),nrn_data{i}(k,:));
            plot_errsi(k) = Rk(1,2)^2; 
        end
    else
        plot_errsi = []; % placeholder, set below if per_err_med
    end
    plot_errs{i} = plot_errsi;
end
%% Plot errs
fig = figure;
ax3 = gca;
args3 = args_all;
args3.ax = ax3;
args3.clims = in.err_cax_lims;
args3.cmap = 'bwr';
% args3.cmap = fake_parula(1000);
args3.fig_name = fig_name3;
if strcmp(in.plot_errs,'per_err_med')
    med_nrn_data = median(cell2mat(nrn_data),2);
    med_dnn_data = median(cell2mat(dnn_data),2);
    plot_errs = {100*(med_dnn_data-med_nrn_data)./med_nrn_data};    
    args3.cell_func = 'none';
    args3.cell_ids{layer_num} = args3.cell_ids{layer_num}(1); 
    args3.cell_model_names = args3.cell_model_names(1);
else
    args3.cell_func = in.err_plot_cell_func; 
end
args3.threshEs = plot_errs;
med_errs = plotThreshLayersFig(layer_set_num,[],args3);
if isempty(in.err_cax_lims)
    cax_lims = caxis;
    fprintf('Color axes: %f - %f\n',cax_lims(1),cax_lims(2));
end
if strcmp(in.plot_errs,'per_err')
    for i = 1:length(plot_errs)
        fprintf('Cell %g\n',i);
        fprintf('%1f %% below %.1f %% error\n',...
            100*sum(plot_errs{i}<in.per_err_cutoff,'all')/numel(plot_errs{i}),in.per_err_cutoff);
    end
    plot_errs_mat = cell2mat(plot_errs);
    fprintf('All: %1f %% below %.1f %% error\n',...
            100*sum(plot_errs_mat<in.per_err_cutoff,'all')/numel(plot_errs_mat),in.per_err_cutoff);
end
% colormap(ax3,bluewhitered(1000));
if in.save_fig
   savefig(fig,fullfile(in.fig_fold,[fig_name3 '.fig']));
   % this export_fig call causes matlab on M1 Pro to hang/crash for some reason...
%    export_fig(fig,fullfile(in.fig_fold,[fig_name3 '.png']),'-png','-cmyk','-r250');
    exportgraphics(fig,fullfile(in.fig_fold,[fig_name3 '.png']),...
                    'Colorspace','rgb','Resolution',250);
   fprintf('Saved %s\n',fig_name3);
end
%% Plot actual vs predicted
if in.plot_actual_v_pred
    fig = figure;
    fig.Units = 'inches';
    fig.Position(3:4) = [2.5 2];
    % fig.Position(3:4) = [6.08 4.13];
    Rsqs = zeros(length(nrn_data),1);
    min_th_actual = 1e3; % for setting axis limits
    max_th_actual = 0;
    for i = 1:length(nrn_data)
        if strcmp(in.output_var,'threshold')
            if strcmp(in.amp_mode,'E_center')
                th_nrni = nrn_data{i}.*scale_factors{i}; % put back in V/m
                th_dnni = dnn_data{i}.*scale_factors{i}; % put back in V/m
                th_nrni = th_nrni(:);
                th_dnni = th_dnni(:);
                unit_str = 'V/m';
            else
                th_nrni = nrn_data{i}(:);
                th_dnni = dnn_data{i}(:);
                unit_str = 'A/\mu s';
            end
        elseif strcmp(in.output_var,'probability')
           th_nrni = nrn_data{i}(:);
           th_dnni = dnn_data{i}(:);
           unit_str = 'prob.';
        end
        R = corrcoef(th_nrni,th_dnni);
        Rsqs(i) = R(1,2)^2;
        plot(th_nrni,th_dnni,'.'); hold on;
        xlabel(sprintf('Actual %s (%s)',in.output_var,unit_str));
        ylabel(sprintf('Predicted %s (%s)',in.output_var,unit_str));
        box off;
        min_th_actual = min([min_th_actual;th_nrni]);
        max_th_actual = max([max_th_actual;th_nrni]);
    end
    % legend(arrayfun(@(x) sprintf('L5 PC %g',x),1:length(cell_ids_all),'UniformOutput',0),'Location','Best')
    if strcmp(in.plot_output_var,'threshold')
    %     title(sprintf('Median absolute percent error: %.4f %%',dnn_errs{1}.medape));
    %     axis([100 500 100 500]);
        plot([min_th_actual max_th_actual],[min_th_actual max_th_actual],'--k');
        xlim([min_th_actual max_th_actual]);
        ylim(xlim);
    else
    %     title(sprintf('Median absolute error: %.4f %%',dnn_errs{1}.mae));
        plot([0 1],[0 1],'--k');
        axis([0 1 0 1]);
    end
    box off;
    % grid on;
    set(gca,'FontName','Arial','FontSize',10);
    fprintf('R squared vals:\n')
    disp(Rsqs);
    if in.save_fig
        savefig(fig,fullfile(in.fig_fold,[fig_name4 '.fig']));
%         export_fig(fig,fullfile(in.fig_fold,[fig_name4 '.png']),'-png','-r300');
        exportgraphics(fig,fullfile(in.fig_fold,[fig_name4 '.png']),...
                    'Resolution',300);
        fprintf('Saved %s\n',fig_name4);
    end
end
end
