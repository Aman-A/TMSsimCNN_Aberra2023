function plot_cbar(save_fig,cmap,cmap_name,cax_lims,ax_scale,cb_size,...
                        cb_size_units,saturate_lims,varargin)
if nargin == 0
   save_fig = 0;  
   % cmap settings
   % cmap = jet(1000);
   cmap = flipud(fake_parula(1000)); % or inferno, parula,etc.
   cmap_name = 'fkparF';
   % cmap_name = 'jet';
   cax_lims = [200 1000];
   % cax_lims = [0.9 1.7];
   % cax_lims = [-5 120];
   ax_scale = 'log'; % or 'linear'
   cb_size = [0.7750    0.0199]; 
   cb_size_units = 'normalized';
   saturate_lims = 1; 
end
% plot settings
in.font_name = 'Arial';
in.font_size = 12;
in.cbar_name = ''; 
in.cbar_location = 'SouthOutside';
in.orientation = 'Horizontal'; % or 'Vertical'
in.fig_fold = 'figures';
in.format_cbar = 0; 
in.units = 'V/m';
in = sl.in.processVarargin(in,varargin); 
%% Plot
fig = figure;
fig.Units = 'centimeters';
fig.Position(3:4) = [20 8]; 
axis off; 
cb = colorbar('Location',in.cbar_location,'Orientation',in.orientation);
cb.FontName = in.font_name; 
cb.FontSize = in.font_size;
cb.Units = cb_size_units;
cb.Position(3:4) = cb_size; 
cb.Color = 'k';
if strcmp(ax_scale,'log')
    caxis(log10(cax_lims))
%     cb.TickLabels = num2str(10.^(cb.Ticks),'%.1f\n');     
    format_cbar_ticks(cb,ax_scale,in); 
else
    caxis(cax_lims)
    if saturate_lims
        if ~any(cb.Ticks == cax_lims(end))
           cb.Ticks = [cb.Ticks cax_lims(end)];           
        end
        cb.TickLabels{end} = sprintf('\\geq%.0f',cb.Ticks(end)); 
    end
end
if strcmp(cmap_name,'rwb')
    colormap(redwhiteblue(1000));
elseif strcmp(cmap_name,'bwr')
    colormap(bluewhitered(1000)); 
else
    colormap(cmap); 
end 
if save_fig
    if strcmp(in.orientation,'Vertical')
        orient_name = 'vert';
    else
        orient_name = 'horz';
    end
    if isempty(in.cbar_name)
        if strcmp(ax_scale,'log')
            cbar_name = sprintf('%s_cb_%s_log_%.3f-%.3f',orient_name,cmap_name,cax_lims(1),cax_lims(2));
        else
            cbar_name = sprintf('%s_cb_%s_%.3f-%.3f',orient_name,cmap_name,cax_lims(1),cax_lims(2));
        end   
    else
        cbar_name = in.cbar_name; 
    end
    savefig(fig,fullfile(in.fig_fold,[cbar_name '.fig']));
    exportgraphics(fig,fullfile(in.fig_fold,[cbar_name '.eps']),...
        'Colorspace','cmyk','BackgroundColor','none');
    fprintf('Saved %s cbar %s\n',in.orientation,cbar_name);
end
end
function format_cbar_ticks(hcb,ax_scale,in,include_units)    
    if nargin < 4
       include_units = 0; 
    end
    set(hcb,'TickDirection','out');
    if include_units
        set(hcb,'TickLabelInterpreter','latex');
    end
    set(hcb,'FontSize',in.font_size);
    ticklabels = get(hcb,'TickLabels');
    for ii=1:length(ticklabels)
        tickvalue=str2double(ticklabels{ii});
        if strcmp(ax_scale,'log')
            if include_units
                ticklabels{ii}=sprintf('$$%4.1f  \\:\\rm{%s}$$',10^tickvalue,in.units);
            else
                 ticklabels{ii}=sprintf('%g',round(10^tickvalue,0));
            end
        else % linear
            if abs(tickvalue)<0.01 % use for small deltaVms
                ticklabels{ii}=sprintf('$$%4.3f  \\:\\rm{%s}$$',tickvalue,in.units);
            else
                ticklabels{ii}=sprintf('$$%4.0f  \\:\\rm{%s}$$',tickvalue,in.units);
            end
        end
    end
    set(hcb,'TickLabels',ticklabels);
end