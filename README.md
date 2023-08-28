#### ****************************************************************
####    Create: 2020-07-01 
####    Last Modified: 2023-08-28
#### ****************************************************************


# diag_feedback_E3SM

#### Description:
- This package is used to run some cloud feedback analysis for E3SM developing versions and CMIP models.

#### Environments:
- python (anaconda)
- xarray 
- cartopy [conda install -c conda-forge cartopy]
- NCO [conda install -c conda-forge nco]
- sklearn [conda install -c conda-forge scikit-learn]
- statsmodels [conda install -c conda-forge statsmodels]
- global-land-mask [pip install global-land-mask]

Installation:

- conda create -n diagfbk -c conda-forge -c xarray dask netCDF4 bottleneck cartopy nco scikit-learn statsmodels
- conda activate diagfbk
- pip install global_land_mask
- cd ./codes 
- f2py -c tropo4d.f90 -m tropo4d [a *.so file will be generated]

#### Pre-needed data
- please download radiative kernel data (Huang et al., 2017) from: https://www.dropbox.com/sh/ngfb7bxhwcbxwu8/AAC6AIha5rLjsl3lUZPiLO6Ua/toa?dl=0&subfolder_nav_tracking=1 first. Then, set RadKernel_dir = *datadir* in main-test.py (for E3SM) or main-test-mip.py (for other CMIP models).

#### How to use it?
##### E3SM raw output
- Modify cases_lookup.py: it is used to map the shortname of each pair of CTL and P4K cases with their long casename. Please add your defined shortname of case as a new dict key with real CTL and P4K experiment names as its value. After that, you can use the shortname in main-test.py to refer the cases.

- main-test.py is the control script. Please modify it following all related settings in that script.
-- Four main options: PreProcess, hasCOSP, RunDiag and GetFigure. Turning on all these will lead to a final webpage on your server's webpage directory. 

##### CMIP model
- main-test-mip.py is the control script. 
-- Three main options: hasCOSP, RunDiag and GetFigure.


#### Outputs?
- If you don't have web sever on your machine, you would check the directory "[shortname-of-the-case]/figure/" to check all generated figures. 
- The directory "data" includes NC files from diagnostics. They can be used to plot figures you want.  
- The directory "diagfigure" includes some figures while running *some* diagnostic codes.  

#### NOTES
- tropo.f90 is a fortran file, so f2py is used to convert it into one .so file as shown above. If some errors occurs, please try to 'module load gcc' first and then execute the command again.


