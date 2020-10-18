# diag_feedback_E3SM
a developing version to diagnose feedbacks for E3SMv1 and developing versions

###########################################################################################################
#### description for v0:
This package is used to run some cloud feedback analysis for E3SM developing versions.


#### Environments:
- python (anaconda)
- latest CDAT
- cartopy

#### How to use it?
- main.py is the control script. Please modify it following all related settings in that script.
- other .py files are defined functions used by main.py. Please don't modify them.

#### outputs?
- the directory "figure" includes some basic output figures.
- the directory "data" includes some NC files and CSV files for further detailed analysis. Please help to tar it and send it to me for further plotting. Thanks ^-^


###########################################################################################################
#### Aug 6, 2020: added plotting scripts as v1

after completing main.py run for all cases, pls modify and run main-plot.py to get several plots for diagnosis.

