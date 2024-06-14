function [v] = invisibleViolin(v)
    % make violin plot invisible to legend    
    % v is a handle to a violinplot
   
v.ScatterPlot.HandleVisibility = 'off';
v.ViolinPlot.HandleVisibility = 'off';
v.BoxPlot.HandleVisibility = 'off';
v.WhiskerPlot.HandleVisibility = 'off';
v.MedianPlot.HandleVisibility = 'off';
v.MeanPlot.HandleVisibility = 'off';
v.NotchPlots(1).HandleVisibility = 'off';
v.NotchPlots(2).HandleVisibility = 'off';

end

