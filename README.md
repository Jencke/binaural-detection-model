This repository contains the code for the binaural-detection model used in the publication arXiv:2111.04637

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5643429.svg)](https://doi.org/10.5281/zenodo.5643429)

# Dependencies
The model depends on the following python packages:
 * numpy
 * scipy
 * pandas
 * tables
 * matplotlib

All of witch should be included in popular python distributions such as [Anaconda](https://www.anaconda.com/). If you are using virtual environments with pip, you can install all requirements by running:

`pip install -r requirements`

# Repository Structure
The repository contains all scripts to run and plot the experiments discussed in the manuscript. It also contains a `data.h5` file which contains pre-calculated results in HDF5 format. There are also two script `calc_all.py` and `plot_all.py` which will run an plot all experiments respectively.

## How to run individual experiments.
The experiment subfolder acts as a python package. Experiments are best loaded individually To calculate and plot the results for the experiment of Lanford & Jeffress 1964 one would thus run:

``` python
from experiments import langford1964

langford1964.calc() # run the experiment store results in data.h5
langford1964.plot() # plot the results which are loaded from the data.h5 file
```

### The `calc` function
Calling the calc function without parameters runs the model with the parameters as stated in the manuscript. Model parameters can, however, be easily changed by setting the parameters

 * rho_hat
 * bin_noise
 * mon_noise

for example:
``` python
langford1964.calc(rho_hat=0.95, bin_noise=0.33, mon_noise=0.70)
```

Be aware that the calc function overwrites previous results that might be stored in `data.h5` to prevent this, set the save parameter to False:

``` python
langford1964.calc(rho_hat=0.95, bin_noise=0.33, mon_noise=0.70, save=False)
```

Alternatively, one can also provide the filename for a new buffer file:

``` python
langford1964.calc(rho_hat=0.95, bin_noise=0.33, mon_noise=0.70, save='newdata.h5')
```

### The `plot` function
As the name suggests, the plot function plots the model results. By default, the function plots pre-calculated values as stored in the `data.h5` file. One can provide the `file` paramter to load data from another file:

``` python
langford1964.plot(file='newdata.h5')
```

## Model Structure
All model code is contained within the `experiments` folder. The actual model is implemented in `model.py`.

Individual experiments are split into subfolders named following the structure `authorYEAR`. The folder `langford1964` thus contains scripts for the experiment of Langford & Jeffress 1964. Each of these folders contains a `calc.py` file which includes the code for running the calculations and saving the results in a buffer file called `data.h5`. The `plot.py` file in the subfolder then contains the code for plotting the results from the buffer file as well as the experimental results.

Please be aware that these files only provide the functions for calculating and plotting the results and can not be called directly.
