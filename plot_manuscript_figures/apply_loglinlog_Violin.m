function apply_loglinlog_Violin(violin_handle,linlimit)

handle_fields = {'ScatterPlot','ScatterPlotOut','ViolinPlot','BoxPlot','WhiskerPlot',...
           'MedianPlot','NotchPlots','MeanPlot'};
for i = 1:length(handle_fields)        
    handle_fieldi = handle_fields{i};
    if isprop(violin_handle,handle_fieldi) && all(ishandle(violin_handle.(handle_fieldi)))
        if length(violin_handle.(handle_fieldi)) == 1
            violin_handle.(handle_fieldi).YData = ...
                loglinlogtransform(violin_handle.(handle_fieldi).YData,linlimit);     
        else
            for j = 1:length(violin_handle.(handle_fieldi))
                violin_handle.(handle_fieldi)(j).YData = ...
                    loglinlogtransform(violin_handle.(handle_fieldi)(j).YData,linlimit);  
            end
        end
    end
end
end