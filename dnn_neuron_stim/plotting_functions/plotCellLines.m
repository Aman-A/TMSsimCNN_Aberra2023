function [C,cell_surf] = plotCellLines(varargin)
%plotCellLines plot cell morphology using lines (edges of surface)
%
%   [C,cell_surf] = plotCellLines(varargin);
%
%   Usage Notes
%   ------
%   Input cell_id and nrn_model_ver to plot neuron with default coordinates
%   at origin (in µm) or to plot population of same cell, load cell data
%   externally and input with input origin, normal vector, and azimuthal
%   rotation
%
%   Inputs
%   ------
%
%   Optional Inputs
%   ---------------
%   cell_id : integer 
%             cell_id corresponding to specific cell_model_name      
%   nrn_model_ver : string
%                   name corresponding to set of neuron model parameters,
%                   specified in outputMorphParams
%   inputting cell_id and/or nrn_model_ver alone are sufficient for plotting 
%   
%   cell_data : struct with following fields (also output by
%               loadCellData.m)
%       C : num_comp x 3 array
%           coordinates of compartments for cell (default is to load from
%           cell_data/<nrn_model_ver>)
%       comp_types : num_comp x 1 vector
%                   compartment type of each compartment (default is to load):
%                   0 - soma, 1 - axon, 2 - node, 3 - myelin, 4 - unmyelin, 
%                   5 - basal, 6 - apic
%                   used for coloring if vals is empty
%       parent_inds : num_comp-1 x 1 vector
%                   index of parent compartment for each compartment,
%                   starting from 2nd compartment (1st coordinate is always
%                   root soma)
%   vals : num_comp x 1 vector
%          values for coloring neurons. (default is to use comp_types).
%          Otherwise can be used to plot potential/E-field distribution on
%          compartments, or any other set of values. 
%   cell_origin : 3 x 1 vector
%                 position to shift cell to in mm (default is no shift)
%   cell_normal : 3 x 1 vector
%                 unit vector corresponding to somatodendritic axis/element
%                 normal in layer mesh for rotating cell (default is no
%                 rotation)
%   phi : scalar
%         azimuthal rotation of cell model (default is 0)
%   lw : scalar
%        line width of lines used to plot morphology (default is 0.5)
%   Examples
%   ---------- 
%   1) **return coordinates and plot default cell position/orientation**
%   cell_id = 6; % 'L23_PC_cADpyr229_1'
%   nrn_model_ver = 'maxH'; 
%   C = plotCellLines(cell_id,nrn_model_ver); 
%   
%   2) ** pre-load cell data and plot at default cell position/orientation**
%   cell_id = 6; % 'L23_PC_cADpyr229_1'
%   nrn_model_ver = 'maxH'; 
%   cell_data = loadCellData(cell_id,nrn_model_ver); % loads C, comp_types, and parent_inds
%   figure; plotCellLines(cell_id,nrn_model_ver,'cell_data',cell_data); 
%   
%   3) ** pre-load cell_data and input values for coloring compartments
%   cell_data = loadCellData(cell_id,nrn_model_ver); % loads C, comp_types, and parent_inds
%   vals = cell_data.C(:,3); % use z coordinate
%   figure;
%   plotCellLines(cell_id,nrn_model_ver,'cell_data',cell_data,'vals',vals);
%   
%   4) ** plot specific cell within NeuronPop ** 
%   cell_id = 6; 
%   nrn_model_ver = 'maxH';
%   layer_set_num = 1; 
%   l_ind = 2; % layer index
%   c_ind = 1; % cell index within layer
%   pos_ind = 1500; % position index
%   layers = loadLayers(layer_set_num);
%   cell_data = loadCellData(cell_id,nrn_model_ver); % loads C, comp_types, and parent_inds
%   cell_data.C = cell_data.C*1e-3; % convert to mm
%   cell_origin = layers(l_ind).cell_origins(pos_ind,:);
%   cell_normal = layers(l_ind).cell_normals(pos_ind,:); 
%   nrn_pop_name = 'nrn_pop1';
%   NeuronPop = loadNeuronPop(nrn_pop_name,nrn_model_ver); 
%   phi = NeuronPop.phis{l_ind}{c_ind}(pos_ind);
%   figure;
%   plotCellLines(cell_id,nrn_model_ver,'cell_data',cell_data,'cell_origin',...
%                   cell_origin,'cell_normal',cell_normal,'phi',phi);
%   axis equal; view([0 0]); axis tight; 
in.cell_id = 6;
in.nrn_model_ver = 'maxH';
in.cell_data = struct(); % default load
in.comp_types = [];
in.parent_inds = [];
in.vals = []; % place holder, replace with comp_types later if not changed by user 
in.axon_color_mode = 'red'; 
in.cell_origin = []; % no shift by default
in.cell_normal = []; % no rotation by default
in.phi = []; % no rotation by default
in.lw = 0.5;
in.plot_soma_sphere = 0;
in.soma_diam = []; 
in = sl.in.processVarargin(in,varargin);
% cell_data_files = fieldnames(in.cell_data);
cell_data_files = {'C','comp_types','parent_inds'}; 
% Load non-input morphology data
% load_data_files = cell_data_files(cellfun(@(x) isempty(in.cell_data.(x)),cell_data_files)); 
load_data_files = cell_data_files(cellfun(@(x) ~isfield(in.cell_data,x),cell_data_files)); 
if ~isempty(load_data_files)
    cell_data = loadCellData(in.cell_id,in.nrn_model_ver,load_data_files);
    for i = 1:length(load_data_files)
       in.cell_data.(load_data_files{i}) = cell_data.(load_data_files{i});  
    end
end
if isempty(in.vals) % Use compartment types for coloring
   in.vals = in.cell_data.comp_types; % use compartment types 
   if strcmp(in.axon_color_mode,'red')
       in.vals(in.vals==2|in.vals==3|in.vals==4) = 1; % set all axonal to 1       
   else
       in.vals(in.vals==3) = 0; % myelin black
       in.vals(in.vals==2) = 1; % node red   
       in.vals(in.vals==4) = 4; % unmyelin orange
   end   
   in.vals(in.vals==5) = 2; % apical dend blue
   in.vals(in.vals==6) = 3; % basal dendrite green
   use_rgb_cmap = 1;
   use_sing_color = 0;
elseif length(in.vals) == 1 && isnumeric(in.vals) % single value to color morphology
    in.vals = in.vals*ones(length(in.cell_data.C),1);
    use_rgb_cmap = 0;
    use_sing_color = 0;
elseif length(in.vals) == 3 || ischar(in.vals) % color as [r,g,b] or 'color name'
    use_sing_color = 1;
    use_rgb_cmap = 0; 
else % use input color values
   use_rgb_cmap = 0; 
   use_sing_color = 0;
end
if ~isempty(in.cell_origin) && ~isempty(in.cell_normal) % shift/rotate cell accordingly    
    if isempty(in.phi)
        fprintf('WARNING: No phi input, using randomized azimuthal rotation\n'); 
        C = placeCell(in.cell_origin,in.cell_normal,in.cell_data.C); 
    else
        C = placeCell(in.cell_origin,in.cell_normal,in.cell_data.C,in.phi); 
    end
else
    C = in.cell_data.C; % use loaded/input C
    if ~isempty(in.phi) % just input phi, no cell_origin/cell_normal
        C = placeCell(C(1,:),[0 0 1],C,in.phi); 
    end
end
% colors = ['r','r','r','r','g','b'];
ax = gca;
hold(ax,'on'); 
% compartment sec types 
% 0 - soma, 1 - axon, 2 - node, 3 - myelin, 4 - unmyelin, 5 - basal, 6 - apic

if in.plot_soma_sphere % Plot sphere at center of soma
    if isempty(in.soma_diam)
        soma_diam = mean(in.cell_data.diams(in.cell_data.comp_types == 0));
    else
        soma_diam = in.soma_diam;
    end
    soma_points = C(in.cell_data.comp_types==0,:);
    soma_point = mean(soma_points,1);
%     soma_point = soma_points(end,:);
    [Xs,Ys,Zs] = sphere(50); 
    soma = surf(ax,Xs*soma_diam+soma_point(1),Ys*soma_diam+soma_point(2),...
                Zs*soma_diam+soma_point(3),in.vals(1)*ones(size(Zs)),...
                'Edgecolor','none','FaceColor','flat','FaceLighting','gouraud');
end
num_soma_comp = sum(in.cell_data.comp_types == 0);
num_comp = size(C,1);
Cp = zeros(3,(num_comp-num_soma_comp)*3); % 2 coordinates per compartment (not incl soma)
% make coordinate array for plotting
Cp(:,1:3:end) = C(in.cell_data.parent_inds(num_soma_comp:end),:)';
Cp(:,2:3:end) = C((1+num_soma_comp):end,:)';
Cp(:,3:3:end) = nan(3,num_comp-num_soma_comp);
if use_sing_color
    surface(ax,[Cp(1,:);Cp(1,:)],[Cp(2,:);Cp(2,:)],[Cp(3,:);Cp(3,:)],...
      'FaceColor','flat','EdgeColor',in.vals,'LineWidth',in.lw);
else
    valsp =  zeros(1,(num_comp-num_soma_comp)*3);
    valsp(1:3:end) = in.vals((1+num_soma_comp):end); 
%     valsp(1:3:end) = in.vals(in.cell_data.parent_inds);
    valsp(2:3:end) = in.vals((1+num_soma_comp):end);
    valsp(3:3:end) = nan(1,num_comp-num_soma_comp);
    cell_surf = surface(ax,[Cp(1,:);Cp(1,:)],[Cp(2,:);Cp(2,:)],[Cp(3,:);Cp(3,:)],[valsp;valsp],...
      'FaceColor','flat','EdgeColor','flat','LineWidth',in.lw); 
end
if use_rgb_cmap % use r,g,b for coloring axon, basal, apical dendrs
    if strcmp(in.axon_color_mode,'red')
        colormap(gca,[1 0 0;0 1 0;0 0 1]);
        caxis(gca,[0 4]);
    else
        colormap(gca,[0 0 0; 1 0 0;0 1 0;0 0 1;1 0.65 0]);
        caxis(gca,[0 5])
    end
end
if ~all(ax.DataAspectRatio==[1,1,1]) % axis equal not already set
    axis off; axis equal; view([0 0]); axis tight; % preferred viewing options
end
end
