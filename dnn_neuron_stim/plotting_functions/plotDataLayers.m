function data_layer = plotDataLayers(layers,data,cell_model_names,cell_ids,varargin)
% Plots median value of simulation data (thresholds,time constants, etc.)
% on layer surfaces, shifts each layer by shift_dir   
in.shift_dir = [0,-35,0]; % default shift
in.plot_vertices = 1; % interpolate data onto vertices (cells at element centers)
in.func = 'median';
in.norm_mode = 'none'; % no normalization
in.scale_data = 1;
in.cmap_scale = 'linear'; % 'linear' or 'log'
in.sim_layerROI_name = '';
in.ROIi = []; 
in.ax = []; 
in.plot_region = 'all';
in = sl.in.processVarargin(in,varargin);
%% Convert to vertex data
% if in.plot_vertices
%     data = getAllVertexCData(data,layers,cell_ids,cell_model_names);
% end
%% Plot single figure with all layer surfaces
data_layer = calcDataLayers(data,cell_model_names,cell_ids,...
                                'func',in.func,'norm_mode',in.norm_mode);
plot_layers = find(cellfun(@(x) ~isempty(x),cell_ids,'UniformOutput',1)); 
data_layer = data_layer(plot_layers); 
layers = layers(plot_layers); 
%% Apply scaling factor
if in.scale_data ~= 1
   data_layer = cellfun(@(x) x*in.scale_data,data_layer,'UniformOutput',0);
   fprintf('Scaled data_layer by %.1f\n',in.scale_data); 
end
%% Extract ROI surface if data was simulated for just this ROI
if ~isempty(in.ROIi)
    ROIi = in.ROIi;
    extract_roi = 1;
elseif ~isempty(in.sim_layerROI_name)
    layerROIdata = loadLayerROI(layers(1).layer_set_num,in.sim_layerROI_name); 
    ROIi = layerROIdata.ROIi; 
    extract_roi = 1;
else
    extract_roi = 0;
end
% num_layers = length(cell_ids); 
if extract_roi
   layers2 = layers; 
   ROIi = ROIi(plot_layers); 
   for i = 1:length(plot_layers)
       [lv,lf] = removeMeshFaces(layers(i).surface.vertices,layers(i).surface.faces,~ROIi{i});
       layers2(i).surface.vertices = lv;
       layers2(i).surface.faces = lf;
       if length(data_layer{i}) ~= size(layers2(i).surface.faces,1)
          data_layer{i} = data_layer{i}(ROIi{i});  
       end
   end
   layers = layers2; 
end
%% Get data points on vertices using weighted average
if in.plot_vertices    
    for i = 1:length(plot_layers)
%         if ~isempty(cell_ids{i})
        if isempty(data_layer{i}) % plot transparent placeholder layer if no data
            data_layer{i} = nan(size(layers(i).surface.vertices,1),1); 
        else
            data_layer{i} = getVertexCData(data_layer{i},layers(i).surface,layers(i).cell_origins); 
        end
%         end
    end
end
%% Check if masking elements based on parcellation data and get region indices
if isfield(layers(1),'region_labels') % only check if layer set has label data
    if ischar(in.plot_region)
        if strcmp(in.plot_region,'all')
            mask_elements = 0;
        else % input as names
            include_region_inds = find(strcmp(in.plot_region,layers(1).region_labels),1);
            if isempty(include_region_inds)
                mask_elements = 0;
            else
                mask_elements = 1;
            end
        end
    elseif iscell(in.plot_region)
        [~,~,include_region_inds] = intersect(in.plot_region,layers(1).region_labels,...
                                                'stable');
        if isempty(include_region_inds)            
            mask_elements = 0;             
            fprintf('%s is not a valid region\n',in.plot_region{:});
        else            
            mask_elements = 1;
        end
    elseif isnumeric(in.plot_region) % input as indices
        mask_elements = 1;   
        include_region_inds = in.plot_region; 
    end
else
    mask_elements = 0;
end
%% Plot
if isempty(in.ax)
    fig = figure;
    ax = gca;
else
    ax = in.ax;
end
for i = 1:length(plot_layers)
%     num_cells_in_layer = length(cell_ids{i});    
%     if num_cells_in_layer >= 1                                    
    layer_surf = layers(i).surface;
    layer_surf.vertices = layer_surf.vertices + (i-1)*repmat(in.shift_dir,size(layer_surf.vertices,1),1);
    p = patch(ax,layer_surf);
    if strcmp(in.cmap_scale,'log')
        p.FaceVertexCData = log10(data_layer{i});
    else
        p.FaceVertexCData = data_layer{i};
    end
    if in.plot_vertices
        p.FaceColor = 'interp';
    else
        p.FaceColor = 'flat';
    end
    p.CDataMapping = 'scaled';
    p.EdgeColor = 'none';
    hold on;     
    if mask_elements
        if in.plot_vertices
            alphas = zeros(size(layer_surf.vertices,1),1);
            for j = 1:length(include_region_inds)
                alphas(layers(i).region_vinds == include_region_inds(j)) = 1;
            end
            p.FaceAlpha = 'interp';
        else
            alphas = zeros(size(layer_surf.faces,1),1); 
            for j = 1:length(include_region_inds)
                alphas(layers(i).region_finds == include_region_inds(j)) = 1;
            end
            p.FaceAlpha = 'flat';
        end
        p.FaceVertexAlphaData = alphas;      
        p.AlphaDataMapping = 'none';
    end
%     else
%         fprintf('No cells in layer %g\n',i);
%     end
end            
camlight(ax); lighting(ax,'gouraud');
axis(ax,'equal','tight','off'); 
end