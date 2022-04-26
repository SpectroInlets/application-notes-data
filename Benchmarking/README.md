## This script enables reproducing the figures from the Spectro Inlets Technical Notes 1 and 2 (Gas Exchange and System Benchmarking)


### How to work with this script?

In the top part of the script select which data to load, which part of the data treatment to execute and how and where to save the resulting plots. 

To use test data included with the script, select following settings: 
```python
DATA_SOURCE = "ixdat" 
data_directory = THIS_DIR / "data"
```

Select the other settings according to how to plot the data.

### Plot your own data
If you want to plot your benchmarking data, first the data needs to be imported using `DATA_SOURCE="raw"`. 
Select the `DATA_DIRECTORY` where your data is saved. If you provide the full path using Windows it can be necessary to give the path as a raw string, i.e. use following syntax: `r"C:\Users\MyUser\MyFolder"`

The first time you run the analysis, you will need to find the time spans for the relevant parts of the data (i.e. CVs, the HER sequence, the gas exchange sequence). This can for example be done by running the script (which might result in error messages that can be ignored for now), and then calling `full_data = main()` in the Python console. This will save the imported dataset in a variable that can then be called directly from the console. 

Now, use `full_data.plot_measurement(tspan=[t1, t2])` to plot the experiment over different time spans. If you use an interactive Python plot, you can leftclick on the plot to get the time printed to the console. You can also left click at one point and right click at another point and then get the time span printed to the console. Double clicking will remove the marker in the plot. 

Once you have found the relevant time spans, go through the script and replace the time span whereever `tspan` is given accordingly. 

The script can automatically detect the onset of a decay. However, the algorithm is not very robust with regards to noise in the data. This means that it can be necessary to adjust the "gradient_cutoff" manually. To help with this, when executed, the script will plot (but not save), the gradient of the signal of interest vs. time. For defining the onset of the decay, choose an absolute cutoff value that is high enough to distinguish from noise. Note that for noisy data, selecting the approprioate cutoff value can be challenging.

Happy analysing!


