

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

#### How to use it?
- main.py is the control script. Please modify it following all related settings in that script.
- other xxx.py files are defined functions used by main.py. Please don't modify them.
- main-plot.py is used to get several plots for diagnosis. 
- additionally, tropo.f90 is a fortrain file. so, use f2py to convert it into one .so file: f2py -c tropo.f90 -m tropo

#### outputs?
- the directory "figure" includes some basic output figures.
- the directory "data" includes some NC files and CSV files for further detailed analysis.


#### history:
- Aug 6, 2020: added plotting scripts as v1
- Sep 30, 2020: use different color for LW,SW,NET in CldRadKernel plot; include amipFuture as a comparison
- Nov 10, 2020: update input data for other CMIP models
