#### ****************************************************************
####    Create: 2020-07-01 
####    Last Modified: 2021-03-24 14:36:23
#### ****************************************************************


# diag_feedback_E3SM
- a developing version to diagnose feedbacks for E3SMv1 and developing versions

#### description:
- This package is used to run some cloud feedback analysis for E3SM developing versions.

#### Environments:
- python (anaconda)
- latest CDAT
- cartopy [conda install -c conda-forge cartopy]
- NCO [conda install -c conda-forge nco]
- sklearn [conda install -c conda-forge scikit-learn]
- statsmodels [conda install -c conda-forge statsmodels]
- global-land-mask [pip install global-land-mask]

#### Pre-needed data
- please download radiative kernel data (Huang et al., 2017) from: https://www.dropbox.com/sh/ngfb7bxhwcbxwu8/AAC6AIha5rLjsl3lUZPiLO6Ua/toa?dl=0&subfolder_nav_tracking=1 first. Then, set RadKernel_dir = *datadir* in main.py.

#### How to use it?
- main.py is the control script. Please modify it following all related settings in that script.
- main-plot.py is used to get several plots for diagnosis. 
- additionally, tropo.f90 is a fortran file. so, use f2py to convert it into one .so file: f2py -c tropo.f90 -m tropo [If some errors occurs, please try to 'module load gcc' first and then execute the command again.]

- other cal_xxx.py files are defined functions used by main.py. Please don't modify them.

#### outputs?
- the directory "figure" includes some basic output figures.
- the directory "data" includes some NC files and CSV files generated by main.py and they will be used by main-plot.py for further analysis.


#### history:
- Aug 6, 2020: added plotting scripts as v1
- Sep 30, 2020: use different color for LW,SW,NET in CldRadKernel plot; include amipFuture as a comparison
- Nov 10, 2020: update input data for other CMIP models
- Feb 18, 2021: add plot_tas_latlon
- Mar 21, 2021: add plot_CldRadKernel_latlon 
- Jul 20, 2021: re-organize  
- Oct 10, 2021: clear-up and re-organize main.py and main-plot.py 
- Dec 22, 2021: generate a web page to better organize and visualize all figures.

