function [threshEs,threshEsM,plot_errs] = plotUnifEMapThreshLayer(model_prefix_pre_base,...
                                                                 layer_set_num,...
                                                                 nrn_pop_names,...
                                                                 Efield_name,...
                                                                 varargin)
if nargin == 0
    model_prefix_pre_base = 'tms';
    layer_set_num = 1;
    nrn_pop_names = arrayfun(@(x) sprintf('nrn_pop%g',x),1:12,'UniformOutput',0);
    Efield_name = 'M1_PA_MagVenture_MC_B70';
end
[dnn_dir,dnn_data_dir] = addPaths_dnn_neuron_stim; % make sure DNN_neuron_stim is on path
in.cell_ids = {[];[];[];16:20;[]};
in.nrn_model_ver = 'maxH';
in.mode = 1; % waveform mode
in.waveform_type = 'tms';
in.sim_data_dir = dnn_data_dir;
in.plot_vertices = 0;
in.plot_func = 'median';
in.plot_err_func = []; 
in.plot_errs = 'per_err';
in.per_err_cutoff = 5; % how many positions below this percent error
in.cax_lims = [80 400];
in.err_cax_lims = [];
in.z_lims = [22 52.4057];
in.cmap = flipud(fake_parula(1000));
in.save_fig = 0;
in.amp_mode = 'E_center';
in.fig_fold = fullfile(dnn_dir,'figs');
in.fig_name_suffix = '';
in.plot_actual_v_pred = 0; 
in = sl.in.processVarargin(in,varargin);
fig_name1 = sprintf('surfThresh_map_%s_w%g',Efield_name,in.mode);
fig_name2 = sprintf('surfTh_map_%s_%s_w%g',in.plot_errs,Efield_name,in.mode);
if strcmp(in.amp_mode,'E_center')
    fig_name3 = sprintf('actual_vs_pred_map_Eth_%s_w%g',Efield_name,in.mode);
else
    fig_name3 = sprintf('actual_vs_pred_map_dIdt_%s_w%g',Efield_name,in.mode);
end

if ~isempty(in.fig_name_suffix)
   fig_name1 = [fig_name1 '_' in.fig_name_suffix];  
   fig_name2 = [fig_name2 '_' in.fig_name_suffix]; 
   fig_name3 = [fig_name3 '_' in.fig_name_suffix];  
end
layer_num = ~cellfun(@isempty,in.cell_ids);
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

cell_ids_lin = [in.cell_ids{:}];
cell_model_names = cellModelNames(cell_ids_lin);

thresh_data = load(fullfile(in.sim_data_dir,'nrn_sim_data',[model_prefix '.mat']));
[~,cell_inds] = intersect(thresh_data.cell_model_names,cell_model_names);
threshEs = thresh_data.threshEs(cell_inds); % double array of thresh in A/�s
fprintf('Loaded neuron simulation data\n');
%% Load unifE map
map_layer_set_num = layer_set_num;
map_model_prefix = 'utms_dt5_maxH_w1_dth5_dph5';
nrn_model_ver = 'maxH';
% map_Efield_name = 'M1_PA_MagVenture_MC_B70'; % was used to generate existing maps
map_Efield_name = Efield_name;
map_nrn_pop_names = arrayfun(@(x) sprintf('nrn_pop%g',x),1:12,'UniformOutput',0);
map_cell_ids = {1:5,6:10,11:15,16:20,21:25};
% map_nrn_pop_names = nrn_pop_names;
 % getMapsSomaE from simTBSnrn/mat
Maps = getMapsSomaE(map_model_prefix,map_layer_set_num,map_Efield_name,map_nrn_pop_names,...
                    nrn_model_ver,'threshold',map_cell_ids);
fprintf('Loaded threshEs using soma E lookup\n');
% Load layersE for E-field
layersE = loadLayers(layer_set_num,'opt','layersE','Efield_name',Efield_name);
if in.mode < 0 % reversed current, flip layersE Efield vecs
    reverse_E = 1;
else
    reverse_E = 0;
end
threshEsM = mapEstLayers(Maps,layersE,in.cell_ids,nrn_pop_names,'reverse_E',reverse_E);
% threshEsM_layer = calcDataLayers(threshEsM,cell_model_names,in.cell_ids,...
%                                  'func',in.func,'norm_mode','none');

%% Plot unifE approximation data
cell_ids1_lin = [in.cell_ids{:}];
if length(cell_ids1_lin) == 1
   plot_cell_model_names = {cellModelNames(cell_ids1_lin)};
else
   plot_cell_model_names = cellModelNames(cell_ids1_lin);
end
if isempty(in.plot_err_func)
   in.plot_err_func = in.plot_func;  
end
args_all = struct();
args_all.cell_ids = in.cell_ids;
args_all.cell_func = in.plot_func;
args_all.clims = in.cax_lims;
args_all.cmap = in.cmap;
args_all.z_lims = in.z_lims;
args_all.fig_fold = in.fig_fold;
args_all.cell_model_names = plot_cell_model_names;
args_all.cbar_on = 0;
args_all.plot_vertices = in.plot_vertices;
args_all.save_fig = in.save_fig;
figure;
ax = gca;
args1 = args_all;
args1.ax = ax;
args1.fig_name = fig_name1;
args1.threshEs = threshEsM;
plotThreshLayersFig(layer_set_num,[],args1);
fprintf('Plotted threshold with func %s across %g clones and %g rotations\n',...
        in.plot_func,length(threshEsM),length(nrn_pop_names));
%% Get errs
unif_errs = cellfun(@(x,y) calc_errs(x,y,'all'),threshEs,threshEsM,'UniformOutput',0);
fprintf('Error of unif E approx at estimating firing probs across clones and %g rotations\n',...
        length(nrn_pop_names));
for i = 1:length(unif_errs)
    unif_errsi = unif_errs{i};
    err_names = fieldnames(unif_errsi);
    fprintf('%s:\n',plot_cell_model_names{i});
    for j = 1:length(err_names)
        if ~strcmp(err_names{j},'errs')
            fprintf('%s: %f\n',err_names{j},mean(unif_errsi.(err_names{j})));
        end
    end
end
if strcmp(in.plot_errs,'err')
    plot_errs = cellfun(@(x) x.errs,unif_errs,'UniformOutput',0);
elseif strcmp(in.plot_errs,'per_err')
    % percent error at each position/rotation for each clone
    plot_errs = cellfun(@(x,y) 100*x.errs./y,unif_errs,threshEs,'UniformOutput',0);
    for i = 1:length(plot_errs)
        fprintf('%s:\n',plot_cell_model_names{i});
        fprintf('%1f %% below %.1f %% error\n',...
            100*sum(plot_errs{i}<in.per_err_cutoff,'all')/numel(plot_errs{i}),in.per_err_cutoff);
    end
    plot_errs_mat = cell2mat(plot_errs);
    fprintf('All: %1f %% below %.1f %% error\n',...
            100*sum(plot_errs_mat<in.per_err_cutoff,'all')/numel(plot_errs_mat),in.per_err_cutoff);
end
%% Plot errs
args2 = args_all;
args2.clims = in.err_cax_lims;
args2.cmap = 'bwr';
args2.fig_name = fig_name2;
if strcmp(in.plot_errs,'per_err_med')
    med_thresh = median(cell2mat(threshEs),2);
    med_threshM = median(cell2mat(threshEsM),2);
    plot_errs = {100*(med_threshM - med_thresh)./med_thresh};
    args2.cell_func = 'none';
    args2.cell_ids{layer_num} = args2.cell_ids{layer_num}(1); 
    args2.cell_model_names = args2.cell_model_names(1);
else
    args2.cell_func = in.plot_err_func; 
end
args2.threshEs = plot_errs;

figure;
ax = gca;
args2.ax = ax;
plotThreshLayersFig(layer_set_num,[],args2);
if isempty(in.err_cax_lims)
    cax_lims = caxis;
    fprintf('Color axes: %f - %f\n',cax_lims(1),cax_lims(2));
end
fprintf('Plotted %s with func %s\n',in.plot_errs,in.plot_err_func);
%% Plot actual vs predicted
if in.plot_actual_v_pred    
    if size(layersE(layer_num).Efield,2) == 3
        Emag_layer = vmag(layersE(layer_num).Efield);
    else
       Emag_layer = vmag(layersE(layer_num).Efield(:,4:6));
    end
    fig = figure;
    fig.Units = 'inches';
    % fig.Position(3:4) = [6.08 4.13];
    fig.Position(3:4) = [2.5 2];
    Rsqs = zeros(length(threshEs),1);
    min_th_actual = 1e3; % for setting axis limits
    max_th_actual = 0;
    for i = 1:length(threshEs)
        if strcmp(in.amp_mode,'E_center')
            thi = threshEs{i}.*Emag_layer; % put back in V/m
            thMi = threshEsM{i}.*Emag_layer; % put back in V/m
            thi = thi(:);
            thMi = thMi(:);
            unit_str = 'V/m';
        else
            thi = threshEs{i}(:);
            thMi = threshEsM{i}(:);
            unit_str = 'A/\mu s';
        end
        R = corrcoef(thi,thMi);
        Rsqs(i) = R(1,2)^2;
        plot(thi,thMi,'.'); hold on;
        xlabel(sprintf('Actual threshold (%s)',unit_str));
        ylabel(sprintf('Predicted threshold (%s)',unit_str));
        box off;
        min_th_actual = min([min_th_actual;thi]);
        max_th_actual = max([max_th_actual;thi]);
    end
    plot([min_th_actual max_th_actual],[min_th_actual max_th_actual],'--k'); 
    % legend(arrayfun(@(x) sprintf('L5 PC %g',x),1:length(cell_ids_lin),'UniformOutput',0),'Location','Best')
    % title(sprintf('Median absolute percent error: %.4f %%',unif_errs{1}.medape));
    % xlim([100 500]); ylim([100 500]);
    ax = gca;
    ax.XLim = [min_th_actual max_th_actual];
    ax.YLim = ax.XLim;
    box off;
    % grid on;
    set(ax,'FontName','Arial','FontSize',10);
    fprintf('R squared vals:\n')
    disp(Rsqs);
    if in.save_fig
        savefig(fig,fullfile(in.fig_fold,[fig_name3 '.fig']));
        export_fig(fig,fullfile(in.fig_fold,[fig_name3 '.png']),'-png','-r300');
    end
end
end