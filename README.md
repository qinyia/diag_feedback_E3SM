

# diag_feedback_E3SM
<<<<<<< HEAD
- a developing version to diagnose feedbacks for E3SMv1 and developing versions

#### description:
- This package is used to run some cloud feedback analysis for E3SM developing versions.
=======
a developing version to diagnose feedbacks for E3SMv1 and developing versions

###########################################################################################################
#### description for v0:
This package is used to run some cloud feedback analysis for E3SM developing versions.

>>>>>>> 798c54ec24f0c0e42f689fdd330245e3f9b4c0df

#### Environments:
- python (anaconda)
- latest CDAT
- cartopy

#### How to use it?
- main.py is the control script. Please modify it following all related settings in that script.
<<<<<<< HEAD
- other xxx.py files are defined functions used by main.py. Please don't modify them.
- main-plot.py is used to get several plots for diagnosis. 

#### outputs?
- the directory "figure" includes some basic output figures.
- the directory "data" includes some NC files and CSV files for further detailed analysis.


#### history:
- Aug 6, 2020: added plotting scripts as v1
- Sep 30, 2020: use different color for LW,SW,NET in CldRadKernel plot; include amipFuture as a comparison
- Nov 10, 2020: update input data for other CMIP models
=======
- other .py files are defined functions used by main.py. Please don't modify them.

#### outputs?
- the directory "figure" includes some basic output figures.
- the directory "data" includes some NC files and CSV files for further detailed analysis. Please help to tar it and send it to me for further plotting. Thanks ^-^


###########################################################################################################
#### Aug 6, 2020: added plotting scripts as v1

after completing main.py run for all cases, pls modify and run main-plot.py to get several plots for diagnosis.

>>>>>>> 798c54ec24f0c0e42f689fdd330245e3f9b4c0df
