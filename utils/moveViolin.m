function v=moveViolin(v,dx)
    % v is a handle to a violinplot
    % move it dx on the x axis
    
v.ScatterPlot.XData = v.ScatterPlot.XData+dx;
v.ViolinPlot.XData = v.ViolinPlot.XData+dx;
v.BoxPlot.XData = v.BoxPlot.XData+dx;
try
    v.WhiskerPlot.XData = v.WhiskerPlot.XData+dx;
catch
    a=1;
end
v.MedianPlot.XData = v.MedianPlot.XData+dx;
v.MeanPlot.XData = v.MeanPlot.XData+dx;


end