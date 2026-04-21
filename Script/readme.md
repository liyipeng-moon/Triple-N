# Code to analysis and generating figure in triple-N paper.

Run CombinedCode.m to go through all figures related to triple-N paper. please extract others/FMRI.zip and others/ModelFeature.zip before running human- and model- related analysis

### Dependency
[GSN toolbox](https://github.com/cvnlab/GSN)  
[colormap_matplotlib](https://github.com/kjamison/colormap_matplotlib_matlab)  
[shadedErrorBar](https://github.com/raacampbell/shadedErrorBar)  
[npy-matlab](https://github.com/kwikteam/npy-matlab)  
[violinplot](https://github.com/bastibe/Violinplot-Matlab/)  
[BombCell](https://github.com/Julie-Fabre/bombcell/)  

### Python version
Here is the python code to reproduce main result from our data, written by Kesheng Wang:
[TripleNpy](https://github.com/JustMuteAll/TripleNpy/tree/main)

### related to plot_EF2_cd.ipynb
Since `npy-matlab` cannot load structured `.npy` files, we provide a Jupyter notebook for visualizing probe drift.
You can install Kilosort4 following the instructions in [here](https://github.com/liyipeng-moon/Triple-N/tree/main/Preprocess).  
### related to plot_F4_PY_plot_surface_data.ipynb
This is notebook to plot any stats into human brain surface. You need to download the subjects' surface file before plotting stats on it.