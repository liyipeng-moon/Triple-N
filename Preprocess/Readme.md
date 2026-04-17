# Preprocessing code for NNN
This repository contains the preprocessing code for NNN data, as well as other datasets collected in BaoLab (PKU) using the [Macaque Fixation task](https://github.com/liyipeng-moon/PassiveViewing_in_ML))

The preprocessing consists of two main steps: 
Python step – Run **Analysis.ipynb** within each session folder to perform spike-sorting analysis without behavioral information. This step runs in Python. 
Matlab step – Run **NSD_Process_pipeline.m** to align the electrophysiology system with the behavioral system. 


### Python Environment Setup
```bash
conda create --name si_kilosort python=3.10
conda activate si_kilosort
pip install spikeinterface[full]
python -m pip install kilosort[gui]
pip uninstall torch
conda install pytorch pytorch-cuda=11.8 -c pytorch -c nvidia
```

### Matlab Dependency:
Add the following tools to the `/util` directory:  
[BombCell](https://github.com/Julie-Fabre/bombcell/)    
[prettify_matlab](https://github.com/Julie-Fabre/prettify_matlab)  
[ShadedErrorBar](https://www.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar)  
[npy-matlab](https://github.com/kwikteam/npy-matlab)  
[MonkeyLogic](https://monkeylogic.nimh.nih.gov/download.html)

### Run `NSD_Process_pipeline.m` to go through 4 steps:
1) Load_Data_function.m  
This function extract basic information from the monkeylogic `.bhv` file, and the meta data from spikeglx `.meta`, extract event code from NI data, and check which trial is valid based on eye signal.  
Output: `META_YYMMDD_Subject_NSD1000_LOC.mat`  
2) PostProcess_function_raw.m  
This function load the output of Kilosort4, and run BombCell to classify unit into noise, single unit, mua, and non-somatic unit.  
Output: `GoodUnitRaw*.mat`  
3) PostProcess_function.m  
This function align the spike into image onset time, generating `GoodUnit` structure consisting raster data and PSTH data  
Output: `GoodUnit*.mat`   
1) PostProcess_function_LFP.m  
This function load LFP signal and align to image onset time.  
Output: `GoodLFP*.mat`  

### Run `S0_ConvertMatTOh5.m` to generate H5 file from GoodUnit, if you don't have it in your downloaded data.