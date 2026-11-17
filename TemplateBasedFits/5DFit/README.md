# 5D fit

Scripts for the different steps to run the 5D (3 angles + mmiss + w) fit.

Paths to tuples/histograms can be change in ```fit_setup.py```

## Creating histograms

To build the histograms needed for the fits run:

```./make_templates.py <sample> <year> <pol> <dim>```

with 

* ```sample```: 
    * The file suffix for the sample of interest.
    * Values: [```'nDstMu'```, ```'nDstTau'```, ```'nbDstst'```, ```'nbsDstst'```, ```'nDstst_higher'```, ```'nDstDX'```, ```'nWSMu'```, ```'data'```]
* ```year```: 
    * The data-taking year of the sample.
    * Values: [2016,2017,2018]
* ```pol```: 
    * Magnet polarity of the sample.
    * Values: [```'up'```, ```'down'```]
* ```dim```: 
    * Dimensionality of the templates.
    * Values: ```'1D'``` for mmiss templates in w bins, ```'5D'``` for mmiss templates in angles + w bins.

Similarly the templates for the control region fit can be produced with :

```./make_templates.py <sample> <year> <pol>```


## Merging histograms

The histograms created with the previous steps are split by magnet polarity and data-taking year. To create merged histograms for the different data-taking years and histograms for the full Run2 dataset, run:


```./merge_up_down.py <year> <dim> <control>```

with

* ```year```: 
    * The data-taking year of the sample. The ```Run2``` option takes as input the histograms produced by the script running for the three data-taking years independently.
    * Values: [2016,2017,2018,Run2]
* ```dim```: 
    * Dimensionality of the templates.
    * Values: ```'1D'``` for mmiss templates in w bins, ```'5D'``` for mmiss templates in angles + w bins.
* ```control```: 
    * Whether to merge the histograms for the control region fit or for the signal region fit.
    * Values: ```'control'``` for the control region histograms. Anything else for the signal region. 


## Running 1D fits

The final (5D) fit takes inputs from simultaneous 1D fits to the signal and control regions to derived constraints on several components. To run the 1D fit do:

```./fit_simultaneous.py -d '1D' -year <year> --constrain_feeddown```

with the following options as baseline:

* ```-year```: 
    * The data-taking year of the sample.
    * Values: [2016,2017,2018,Run2]
* ```-d '1D'```: 
    * For a 1D fit.
* ```--constrain_feeddown```: 
    * To constrain the relative contributions of light B meson feeddown and Bs feeddown backgrounds.

Additinal options include:

* ```-i```: 
    * Which contributions to include in the fit, by default all contributions are included. If ```-i``` is used, only the specified contributions will be included
    * Values: [```'nDstMu'```, ```'nDstTau'```, ```'nbDstst'```, ```'nbsDstst'```, ```'nDstst_higher'```, ```'nDstDX'```, ```'nWSMu'```]
* ```--fix_tau```: 
    * Wheter to fix the B->D*taunu contribution. Default is False.
* ```--constrain_feeddown```: 
    * To constrain the relative contributions of light B meson feeddown and Bs feeddown backgrounds.
* ```--smear```: 
    * Whether to float the vertex smearing parameter in the fit. Default is False.
* ```--toy```: 
    * Wheter to run a toy fit. Default is False.


## Running 5D fits

To run the final independent fit in the different 5D bins, do:

```./fit_indep.py -b <bin> -d '5D' -year <year> --constrain_feeddown```

with the following options as baseline:

* ```-b```: 
    * The phasespace bin number for which to run the fit.
    * Values: [0:224]
* ```-year```: 
    * The data-taking year of the sample.
    * Values: [2016,2017,2018,Run2]
* ```-d '5D'```: 
    * For a fit in the 5D bins.
* ```--constrain_feeddown```: 
    * To constrain the relative contributions of light B meson feeddown and Bs feeddown backgrounds.

Additinal options include:

* ```-i```: 
    * Which contributions to include in the fit, by default all contributions are included. If ```-i``` is used, only the specified contributions will be included
    * Values: [```'nDstMu'```, ```'nDstTau'```, ```'nbDstst'```, ```'nbsDstst'```, ```'nDstst_higher'```, ```'nDstDX'```, ```'nWSMu'```]
* ```--fix_tau```: 
    * Wheter to fix the B->D*taunu contribution. Default is False.
* ```--constrain_feeddown```: 
    * To constrain the relative contributions of light B meson feeddown and Bs feeddown backgrounds.
* ```--smear```: 
    * Whether to float the vertex smearing parameter in the fit. Default is False.
* ```--toy```: 
    * Wheter to run a toy fit. Default is False.
