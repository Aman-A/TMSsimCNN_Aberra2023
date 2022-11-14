%%
function plot_vert_cbars_figs(layer_num,cax_lims,err_cax_lims,save_fig,...
                              plot_inds)
% addpath('..');
% cmap = flipud(fake_parula(1000));
% % threshold 
% plot_horz_cbar(1,cmap,'fkparF',[90 400],'linear',[20 1],...
%                 'centimeters',1,'Arial',18,'figs')
% % error 
% plot_horz_cbar(1,[],'bwr',[-12.92 12.24],'linear',[10 0.8],... % [20 1]
%                 'centimeters',0,'Arial',18,'figs')            
%%             
% plotVertCbar(1,[-0.1 0.4],[],'bwr',[5 15],[0.2 4.5])
if nargin < 5
   plot_inds = [1 2]; % plot both 
end
%% Plot colorbars
cbar_opts.font_name = 'Arial';
cbar_opts.font_size = 10; 
cbar_opts.cbar_location = 'EastOutside';
cbar_opts.orientation = 'Vertical';
cbar_opts.fig_fold = 'figs';
% Threshold plot
if any(plot_inds == 1)
    cbar1_size = [0.1 1.5]; 
    cax_lims1 = cax_lims;
    cmap1 = flipud(fake_parula(1000));
    cmap_name1 = 'fkparF'; 
    cbar_opts.cbar_name= sprintf('cnn_cbar_l%g',layer_num);
    plot_cbar(save_fig,cmap1,cmap_name1,cax_lims1,'linear',cbar1_size,...
                'inches',1,cbar_opts)
end
% Error plot
if any(plot_inds == 2)
    cbar2_size = [0.1 1.5]; 
    cax_lims2 = err_cax_lims; 
    cbar_opts.cbar_location = 'EastOutside';
    cbar_opts.cbar_name = sprintf('err_cnn_cbar_l%g',layer_num); 
    plot_cbar(save_fig,'','bwr',cax_lims2,'linear',cbar2_size,...
                'inches',0,cbar_opts)
end
end