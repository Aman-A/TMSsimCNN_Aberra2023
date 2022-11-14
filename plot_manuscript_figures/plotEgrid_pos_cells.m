% Plots sampling grid points for multiple cells on the same axis
function plotEgrid_pos_cells(cell_ids,nrn_model_ver,cell_origin,cell_normal,...
                            phi,sample_method_struct,shift_vec,save_fig)
if nargin == 0
    cell_ids = 16:20;
    nrn_model_ver = 'maxH';
    cell_origin = [0 0 0];
    cell_normal = [0 0 1];
    phi = 0; 
    sample_method_struct.method = 'box';
    N = 9; 
    sample_method_struct.Nx = N;
    sample_method_struct.Ny = N;
    sample_method_struct.Nz = N;
    l = 2e3; % um
    sample_method_struct.lx = l;
    sample_method_struct.ly = l;
    sample_method_struct.lz = l;
    sample_method_struct.rshift = [0 0 0]; 
%     sample_method_struct.rshift = [0 0 0.2024]*1e3; % [0 0 -0.1041] 
    shift_vec = [sample_method_struct.lx*1.5 0 0]; 
    save_fig = 0;
end
Cs = samplePts(sample_method_struct,'print_level',1); 
min_max_Cz = zeros(length(cell_ids),2);
fig = figure; 
fig.Position(3:4) = [1320 590]; 
for i = 1:length(cell_ids)
    cell_origini = cell_origin + shift_vec*(i-1);
    Ci = plotCellLines('cell_id',cell_ids(i),'nrn_model_ver',nrn_model_ver,...
                  'cell_origin',cell_origini,'cell_normal',cell_normal,...
                  'phi',phi,'lw',1); 
    min_max_Cz(i,:) = [min(Ci(:,3)),max(Ci(:,3))]; 
    Csi = placeCell(cell_origini,cell_normal,Cs,phi); 
    plot3(Csi(:,1),Csi(:,2),Csi(:,3),'k.'); 
    axis equal; 
end
% zlim([min(Csi(:,3)) max(Csi(:,3))])
zlim([min(min_max_Cz(:,1)),max(min_max_Cz(:,2))]); 
if isfield(sample_method_struct,'lx')
   xlim([-sample_method_struct.lx shift_vec(1)*length(cell_ids)])
else
   xlim([-sample_method_struct.lr shift_vec(1)*length(cell_ids)])    
end
axis off;
view([3.8 22.05]); 
if save_fig 
    if isfield(sample_method_struct,'lx')
        fig_name = sprintf('samplepts_cells%g-%g_%s_Nz%g_lx%g_ly%g_lz%g',...
            cell_ids(1),cell_ids(end),sample_method_struct.method(1:3),...
            sample_method_struct.Nz,...
            sample_method_struct.lx*1e-3,... % convert to mm
            sample_method_struct.ly*1e-3,sample_method_struct.lz*1e-3);
    else
        fig_name = sprintf('samplepts_cells%g-%g_%s_Nz%g_lr%g_lz%g',...
            cell_ids(1),cell_ids(end),sample_method_struct.method(1:3),...
            sample_method_struct.Nz,...
            sample_method_struct.lr*1e-3,... % convert to mm
            sample_method_struct.lz*1e-3);
    end    
    printFig(fig,'figs',fig_name,'formats',{'fig','png'},...
            'resolutions',{[],'-r300'});    
end
end