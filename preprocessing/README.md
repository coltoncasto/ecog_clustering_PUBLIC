Contained in this subdirectory are the two main scripts used to run preprocessing. They are intended to processing signal and extract trial information from the **raw** files as we recieve them from the collection sites. Although the raw data are not availible publicly, we have provided the code for transparency. For those interested in tracing the code, a sample command to run preprocessing has been provided below.

```matlab
crunch('AMC026','MITSWJNTask',...
       'fromScratch',true,...
       'doneVisualInspection',false,...
       'order','DefaultECoG',...
       'isPlotVisible',true);
```