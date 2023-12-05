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
- cartopy
- NCO
- sklearn
- statsmodels
- global-land-mask

Installation:

- `conda create -n diagfbk -c conda-forge python=3.11.0 xarray dask netCDF4 bottleneck cartopy nco scikit-learn statsmodels`
- `conda activate diagfbk`
- `pip install global_land_mask`
- `cd ./codes` 
- `f2py -c tropo4d.f90 -m tropo4d` [a *.so file will be generated]

#### Pre-needed data
- please download radiative kernel data (Huang et al., 2017) from [their website](https://www.dropbox.com/sh/ngfb7bxhwcbxwu8/AAC6AIha5rLjsl3lUZPiLO6Ua/toa?dl=0&subfolder_nav_tracking=1) or [this backup link](https://portal.nersc.gov/project/mp193/qinyi/DATA/Huang_kernel_data/) first. Then, set RadKernel_dir = *datadir* in main-test.py (for E3SM) or main-test-mip.py (for other CMIP models).

#### How to use it?
##### For E3SM raw output
1. **Modify cases_lookup.py**: it is used to map the shortname of each pair of CTL and P4K cases with their long casename. Please add your defined shortname of case following the example. After that, you can use the shortname in main-test.py to refer the cases.
2. **Modify main-test.py**: it is the control script.
-- Four main options: PreProcess, hasCOSP, RunDiag and GetFigure. Turning on all these will lead to a final webpage on your server's webpage directory. 

##### For CMIP model
- **Modify main-test-mip.py**: it is the control script for processing CMIP models. 
-- Three main options: hasCOSP, RunDiag and GetFigure.


#### Outputs?
##### Have webserver:
-   Go to the webdir and open the webpage.
##### No webserver: 
-   Download the "**[shortname-of-the-case]/figure/view.tar**" for E3SM or "**[model_name-variant]/figure/view.tar**" for CMIP models to your local computer.
-   `untar view.tar`
-   `cd view`
-   `cd viewer` and open the index.html. The generated webpage will be opened in your browser. 


#### NOTES: 
- tropo.f90 is a fortran file, so f2py is used to convert it into one .so file as shown above. If some errors occurs, please try to 'module load gcc' first and then execute the command again.
- The directory "data" includes NC files from RunDiag. They can be used to plot figures that you want.  
- The directory "diagfigure" includes some figures while running *some* diagnostic codes.  

