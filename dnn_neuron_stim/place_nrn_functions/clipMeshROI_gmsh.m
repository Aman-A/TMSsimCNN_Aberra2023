function [mesh_roi,include_v] = clipMeshROI_gmsh(mesh,ROI,include_data)
if nargin == 2
   include_data = 1; % set to 0 for meshes with no node/element data 
end
% CLIPMESHROI_GMSH clips gmsh mesh within region of interest, defined by 
% ROI = [xmin xmax ymin ymax zmin zmax]
% only includes tetrahedron data and data within GM and WM tissues
% (excludes CSF, skull, and skin)
% applies to potential data if included in mesh
% Extract all tet center coordinates            
tet_centers = mesh_get_tetrahedron_centers(mesh);
% Get inds of points in GM/WM
mesh_inds = mesh.tetrahedron_regions==2 | mesh.tetrahedron_regions==1;
% get inds of points in GM/WM and in ROI
tet_ROIinds = (tet_centers(:,1)>ROI(1)) & (tet_centers(:,1)<ROI(2)) & ...
(tet_centers(:,2)>ROI(3)) & (tet_centers(:,2)<ROI(4)) & ...
(tet_centers(:,3)>ROI(5)) & (tet_centers(:,3)<ROI(6)) & mesh_inds;
% get mesh data in ROI
tetrahedra_roi = mesh.tetrahedra(tet_ROIinds,:); % get tetrahedra in GM/WM and ROI
% get E in GM/WM and ROI
if include_data
    field_idx = get_field_idx(mesh,'E','element');
    E = mesh.element_data{field_idx}.tetdata(tet_ROIinds,:); % [Ex Ey Ez] at tet centers
end
% get renumbered mesh
[node_inds,~,tetrahedra] = unique(tetrahedra_roi(:));
nodes = mesh.nodes(node_inds,:);
tetrahedra = reshape(tetrahedra,size(tetrahedra_roi));
if isempty(mesh.node_data)
    include_v = 0;
else
    include_v = 1;
    v = mesh.node_data{1}.data(node_inds);
end
mesh_roi = mesh_empty;
if include_data
    mesh_roi.element_data = {struct('name','E','tetdata',E,'tridata',[])}; % give placeholder data for tri
    if include_v
        mesh_roi.node_data = {struct('name','v','data',v)};
    end
end
mesh_roi.nodes = nodes;
mesh_roi.tetrahedra = int32(tetrahedra);
mesh_roi.tetrahedron_regions = mesh.tetrahedron_regions(tet_ROIinds);
mesh_roi.triangles = [];
mesh_roi.triangle_regions = [];
end