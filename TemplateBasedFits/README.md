# Fit 

The environment with necessary packages can be accessed by 
```
conda activate /afs/cern.ch/user/m/mcolonna/public/.conda/redist_env
```
For plotting, mplhep pkg is necessary but it can be installed in your base enviroment or can be found in lb-conda environments.



## Configuring Templates
- ```config_bin.json``` files are where we define the configuration of templates with histogram name, binning scheme, FF parameterisation name. These need to be created before running the fit. You can use utils/create_config_files.py script.

- Fit templates already include b-kin-occ and TrackCalib2 corrections for mmiss2 fits. There are 3 different templates/fits : WSPi removed, WSPi included (david_temps_v2), and full-q2 (for post-fit).
- Templates are already inside this directory and created from [https://gitlab.cern.ch/lhcb-slb/b02dstmunu_q2_shape/-/tree/dvicoben/Templates/histos_wbins] 


## Pre-limary steps before running fit
```utils/``` directory contains some scripts used for 
- if base templates change, new .dat files used in the fit can be created using -> ```json_converter_dat.py``` (Please note that we do not change base templates it is just the way of reading them in fit)
- Pre-fit templates can be plotted and their yields can be saved to be used in fit starting points -> ```plotTemplates-and-saveYields.py``` 
- creating config files -> ```create_config_files.py```
- grouping Dstst bkg templates -> ```sum_dstst_temps.py```

# Post-fit
- ```initial-reweight-comparison-full-w.py``` is for comparison of data/mc with b-kin-occ + track calib corrections
- ```3d_reweighting-full-w.py``` is for comaprison of data/mc with additional reweightings (see script for each step of sequential 3d reweighting) 

## Pre-limary steps before running post-fit
- You should have mmiss2 fit first and your fitted yields. For quick comparison, we look at full w fit.
 

