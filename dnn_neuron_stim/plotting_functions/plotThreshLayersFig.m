function [med_thresh,varargout] = plotThreshLayersFig(layer_set_num,model_prefix,varargin)
% Plots median threshold across clones and rotations at each location on
% layers 1-5, adjacent to each other, as in Fig. 3c of Aberra 2020
mat_dir = addPaths_dnn_neuron_stim;
if nargin==0
   model_pre = 'tms';
   layer_set_num = 1;
   nrn_model_ver = 'maxH';
   mode = 1; % monophasic MagProX100 pulse   
   Efield_name = 'M1_PA_MagVenture_MC_B70_ernie';
   nrn_pop = 'nrn_pop1';
   model_prefix = sprintf('%s_%s_w%g_ls_%g_E_%s_P_%s',model_pre,nrn_model_ver,mode,...
           layer_set_num,Efield_name,nrn_pop);   
end
% Optional settings
in.data_fold = fullfile(mat_dir,'nrn_sim_data');
in.cell_ids = {[],[6:10],[11:15],[16:20],[]}; 
in.shift_dir = [0,-35,0]; % vector of direction of shift
in.plot_vertices = 0;
in.plot_region_name = ''; % FDI_rep_inds
% Data processing settings
in.cell_func = 'median'; % function to apply across threshold of cells in each position
in.cutoff = 0; % binary color map, all thresholds below this quantile are one color
in.scale_data = 1; % scale all values by scalar
in.norm_mode = 'none';
% Plot settings
in.ax_view = [-89.2 45]; % [-89.2 70.8]
in.cmap = [flipud(fake_parula(1000));0.8 0.8 0.8]; % add gray for values above cutoff
in.clims = [80 220]; % A/us
in.cmap_scale = 'linear'; % 'linear' or 'log'
in.z_lims = [22 52.4057]; % or []
in.lt_pos = [-411 -807 836]; % [170.3353 52.0246 1.2195e3]
in.ax = []; % axis handle
in.fig_units = 'normalized';
in.fig_size = [1 1];
in.save_fig = 0;
in.fig_name = sprintf('thresh_layers_%s',model_prefix);
in.fig_fold = fullfile(mat_dir,'figures');
in.fig_format = 'png';
in.cbar_on = 1;
in.plot_region = 'all'; % plot sub-region based on label data (region_finds/vinds)
% Input data rather than loading
in.threshEs = []; 
in.cell_model_names = []; 
in = sl.in.processVarargin(in,varargin);
%% Load data
layers = loadLayers(layer_set_num);
if isempty(in.threshEs)
    data_struct = load(fullfile(in.data_fold,[model_prefix '.mat']));
    threshEs = data_struct.threshEs;
    cell_model_names = data_struct.cell_model_names;    
else
    threshEs = in.threshEs;
    cell_model_names = in.cell_model_names; 
    if ischar(cell_model_names)
       cell_model_names = {cell_model_names};  
    end
end
plot_layers = find(cellfun(@(x) ~isempty(x),in.cell_ids,'UniformOutput',1));
%% Plot
opts.func = in.cell_func;
opts.shift_dir = in.shift_dir;
opts.plot_vertices = in.plot_vertices;
opts.scale_data = in.scale_data;
opts.cmap_scale = in.cmap_scale; 
opts.norm_mode = in.norm_mode;
opts.ax = in.ax;
opts.plot_region = in.plot_region;
med_thresh = plotDataLayers(layers,threshEs,cell_model_names,...
                           in.cell_ids,opts);
if isempty(in.ax)
    ax = gca;
else
    ax = in.ax;
end
fig = ax.Parent;
fig.Units = in.fig_units;
fig.Position(3:4) = in.fig_size;
if ~isempty(in.clims)
    if strcmp(in.cmap_scale,'log')
        caxis(ax,log10(in.clims));
    else
        caxis(ax,in.clims);
    end
end
if ischar(in.cmap)
    if strcmp(in.cmap,'bwr')
        colormap(ax,bluewhitered(1000));
    elseif strcmp(in.cmap,'rwb')
        colormap(ax,redwhiteblue(1000));
    end
else
    colormap(in.cmap);
end
if ~isempty(in.ax_view)
    view(ax,in.ax_view);
end
% change light
ax.Children(1).Position = in.lt_pos;
ax.Children(1).Style = 'local';
% cut off sulcus
if ~isempty(in.z_lims)
   ax.ZLim = in.z_lims;
end
% add colorbar
if in.cbar_on
    colorbar(ax,'FontSize',16);
end
%%
if in.cutoff % plot binary plot of thresholds below quant_cutoff quantile
    p_all = cell(length(plot_layers),1);
    if length(p_all) > 1
        for i = 1:length(p_all); p_all{i} = ax.Children(end-i+1); end
    else
        p_all = ax.Children(2); % skip light
    end
%     elems_below = cellfun(@(x) find(x <= min(x)*cutoff),med_thresh,'UniformOutput',0);
    if in.cutoff < 1
        if any(strcmp(in.norm_mode,{'max_exc','exc'}))
%             elems = cellfun(@(x) find(x > quantile(x,in.cutoff)),med_thresh,'UniformOutput',0);
            elems = cellfun(@(x) find(x > in.cutoff),med_thresh,'UniformOutput',0);
        else
            elems = cellfun(@(x) find(x <= quantile(x,in.cutoff)),med_thresh,'UniformOutput',0);
        end
        for i = 1:length(plot_layers)
            p_all{i}.FaceVertexCData = zeros(length(p_all{i}.FaceVertexCData),1);
            p_all{i}.FaceVertexCData(elems{i}) = ones(length(elems{i}),1);
        end
        area_activated = cellfun(@(x,y) sum(getPatchAreas(x.Vertices, x.Faces,y)),p_all,elems,'UniformOutput',1)';
        display(area_activated);
        colormap([.8 .8 .8; 1 0 0]);
        caxis([0 1]);
        varargout{1} = area_activated;
    else
        caxis([1 in.cutoff]);
    end
end
%% Save figs
if in.save_fig    
    printFig(fig,in.fig_fold,in.fig_name,'formats',{'fig',in.fig_format},...
            'resolutions',{[],'-r250'});
end
