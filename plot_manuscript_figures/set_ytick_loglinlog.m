function set_ytick_loglinlog(y_lims,ymax,linlimit,ax)
if nargin < 4
   ax = gca; 
end
% ymax = 100;
ylim(ax, loglinlogtransform(y_lims, linlimit));
decneg = -1:0.1:-0.2;
decpos = 0.2:0.1:1;
n = ceil(log10(ymax));
ytick = [];
yticklabel_vec = [];
while n > log10(linlimit)
    ytick = [ytick,decneg*10^n];
    yticklabel_vec = [yticklabel_vec,-10^n];
    n = n-1;
end
ytick = [ytick, -linlimit:linlimit/2:linlimit];
n = log10(linlimit);
while n < ceil(log10(ymax))
    ytick = [ytick, decpos*10^(n+1)];
    yticklabel_vec = [yticklabel_vec,10^(n+1)];
    n = n+1;
end
yticklabel_vec = [yticklabel_vec, -linlimit:linlimit/2:linlimit];
%     yticklabel_vec = [yticklabel_vec,0];

yticklabels = cell(size(ytick));
[yticklabels{:}]=deal('');
for ii = 1 : length(ytick)
    switch ytick(ii)
        case num2cell(yticklabel_vec)
            yticklabels{ii} = sprintf('$$%g\\%%$$',ytick(ii));
    end
end
set(ax, 'YTick',  loglinlogtransform(ytick, linlimit), ...
    'YTickLabel', yticklabels, 'TickDir', 'out',...
    'TickLabelInterpreter','latex');