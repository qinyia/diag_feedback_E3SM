#!/bin/env python
# -*- coding: utf-8 -*-

"""

Stephen Po-Chedley 19 March 2021

xarray wrangle (xrw): I/O and query utilities for accessing CMIPx data.

@author: pochedls
"""

import sqlite3
import xarray as xr
import numpy as np
#<qinyi 2021-10-25 #------------------
import copy
from xarray.coding.variables import SerializationWarning
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=SerializationWarning)
#>qinyi 2021-10-25 #------------------


def get_cmip_paths(**kwargs):
    """
    get_cmip_paths(**kwargs).

    Function returns a list of data paths based on user-defined search
    criteria (at least one search constraint must be provided). The
    optional arguments, include:

        db : path to sqlite database file
        mip_era : mip_era for CMIP data
        activity : activity for CMIP data
        experiment : experiment for CMIP data
        realm : realm for CMIP data
        frequency : frequency for CMIP data
        variable : variable for CMIP data
        model : model for CMIP data
        member : realization for CMIP data
        gridLabel : grid label for CMIP data
        trim : Boolean to trim off duplicate files (default True)
        verbose : Boolean if information about search should be printed

    Example Usage:
    -------------
        dpaths = get_cmip_paths(model='CCSM4',
                                experiment='historical',
                                variable='tas',
                                mip_era = 'CMIP5',
                                activity = 'CMIP',
                                realm = 'atmos',
                                frequency = 'mon',
                                member = 'r1i1p1')

    Returns:
    -------
        ['/p/css03/esgf_publish/cmip5/output1/NCAR/CCSM4/historical/mon/atmos/Amon/r1i1p1/v20160829/tas/']
    """
    #  Define default dictionary for search parameters
    searchDict = {'mip_era': '*',
                  'activity': '*',
                  'experiment': '*',
                  'realm': '*',
                  'frequency': '*',
                  'variable': '*',
                  'model': '*',
                  'member': '*',
                  'gridLabel': '*',
                  'cmipTable': '*'} #<qinyi 2021-10-31 #------------------

    #  Ensure search arguments were provided
    if len(kwargs.keys()) == 0:
        print('No search criteria provided. Provide search constraints, such '
              'as: \n\n'
              'mip_era, activity, experiment, realm, frequency, '
              'variable, model, member')
        return

    #  Replace default arguments with user-provided arguments
    for key in kwargs:
        if key in searchDict.keys():
            searchDict[key] = kwargs[key]
        if key == 'db':
            db = kwargs['db']
        else:
            db = '/p/user_pub/xclim/persist/xml.db'
        if key == 'verbose':
            verbose = kwargs['verbose']
        else:
            verbose = True
        if key == 'trim':
            trim = kwargs['trim']
        else:
            trim = True
        if key == 'criteria':
            criteria = kwargs['criteria']
        else:
            criteria = ['ver', 'cdate', 'tpoints']

    # Construct sqlite search
    constraints = [key + "='" + searchDict[key] + "'" for key in searchDict.keys() if searchDict[key] != '*']
    # ensure we don't use ignored / retired data
    constraints = " and ".join(constraints) + " and retired = 0 and ignored = 0;"
    query = "select * from paths where " + constraints

    # perform sqlite search
    if verbose:
        print('Searching for data...')
    con = sqlite3.connect(db)
    cur = con.cursor()
    cur.execute(query)
    headers = list(map(lambda x: x[0], cur.description))
    result = cur.fetchall()
    con.close()

    # create a dictionary of paths
    # each path is a key linking to other metadata
    pathDict = {}
    for l in result:
        dpath = l[0]
        pathDict[dpath] = {}
        for i, h in enumerate(headers[1:]):
            pathDict[dpath][h] = l[i+1]

    dpaths = list(pathDict.keys())

    # Deduplicate paths (if specified)
    if trim:
        if verbose:
            print('De-duplicating list...')
        dpaths = trimModelList(pathDict, criteria=criteria, verbose=verbose)

    return dpaths


def trimModelList(pathDict,
                  criteria=['cdate', 'ver', 'tpoints'],
                  verbose=False):
    """
    dpathsOut = trimModelList(pathDict).

    trimModelList takes in a list of file paths and returns a list of
    paths such that there is one path per model and realization.

    The returned dpaths are priorized by a cascading criteria, which can be
    optionally specified. The default order is:
        ver:        prioritizes dpaths based on version id
        cdate:      prioritizes dpaths that were created more recently
        tpoints:    prioritizes dpaths with more time steps
        publish:    prioritizes dpaths that have 'publish' in their path

    The cascading criteria can be altered by specifying an optional
    argument, criteria, with a list of the strings above (e.g.,
    criteria=['publish', 'tpoints', 'ver', 'cdate']).

    An additional optional argument is verbose (boolean), which will output
    diagnostic information during execution. By default verbose is False.
    """
    keyMap = {}
    models = []
    rips = []
    # loop over dpaths and store metadata
    for dpath in pathDict.keys():
        # get initial metadata for sorting
        model = pathDict[dpath]['model']
        rip = pathDict[dpath]['member']
        ver = versionWeight(pathDict[dpath]['version'])
        # collect all models and rips
        models.append(model)
        rips.append(rip)
        # store data in dictionary
        keyMap[dpath] = {'model': model, 'rip': rip, 'ver': ver}
    # get unique models / rips
    rips = list(set(rips))
    models = list(set(models))
    models.sort()
    rips.sort()

    # loop over each model + realization and filter by /criteria/
    dpathsOut = []
    for model in models:
        if verbose:
            print('  ' + model)
        for rip in rips:
            # subset dpaths for each model / realization
            subdpaths = [fn for fn in keyMap.keys() if
                         (keyMap[fn]['model'] == model
                          and keyMap[fn]['rip'] == rip)]
            # continue whittling down path list until only one is left
            # by iteratively using each criteria
            for crit in criteria:
                subdpaths, keyMap = filter_dpaths(subdpaths, keyMap, crit)
            # if more than one path is left after criteria is applied,
            # choose first one
            if len(subdpaths) > 0:
                dpathsOut.append(subdpaths[0])

    # if verbose mode, print off dpaths and selection (*)
    if verbose:
        # loop over models / rips and see if dpaths are in output list
        # print dpaths used with asterisk and unchosen dpaths next
        for model in models:
            for rip in rips:
                subdpaths = [fn for fn in keyMap.keys()
                            if (keyMap[fn]['model'] == model
                                and keyMap[fn]['rip'] == rip)]
                if len(subdpaths) > 1:
                    lowdpaths = []
                    for fn in subdpaths:
                        if fn in dpathsOut:
                            fn1 = fn
                        else:
                            lowdpaths.append(fn)
                    print('* ' + fn1)
                    for fn in lowdpaths:
                        print(fn)
                elif len(subdpaths) == 1:
                    print('* ' + subdpaths[0])

    return dpathsOut


def get_dataset_metadata(dpath):
    """
    cdate, publish, tpoints = getFileMeta(dpath).

    getFileMeta takes a path (dpath) for a CMIP xml and returns:
        cdate:      creation date
        publish:    boolean if the underlying data is in the LLNL publish
                    directories (if it has been locally republished)
        tpoints:    the number of timesteps in the dataset
    """
    try:
        #<qinyi 2021-11-01 #------------------
        if 'HadGEM2-ES' in dpath:
            ds = xr.open_mfdataset(dpath+'*1[8-9]????-??????*')
        else:
            ds = xr.open_mfdataset(dpath + '*', combine='by_coords')
        #>qinyi 2021-11-01 #------------------

    except:
        print('problem opening, path skipped:\n', dpath)
        cdate, publish, tpoints = [0 for _ in range(3)]
        return cdate, publish, tpoints
    if 'creation_date' in ds.attrs:
        cdate = ds.creation_date
    else:
        cdate = '1989-03-06T17:00:00Z'  # Assume PCMDI dawn of time
    # If missing, default to PCMDI dawn of time 9am Monday 6th March 1989
    # most dates are of form: 2012-02-13T00:40:33Z
    # some are: Thu Aug 11 22:49:09 EST 2011 - just make 20110101
    if cdate[0].isalpha():
        cdate = int(cdate.split(' ')[-1] + '0101')
    else:
        cdate = int(cdate.split('T')[0].replace('-', ''))
    # check if republished
    if bool(dpath.find('publish')):
        publish = True
    else:
        publish = False
    axisList = ds.coords
    if 'time' in axisList:
        tpoints = len(ds['time'])
    else:
        tpoints = 0
    ds.close()
    return cdate, publish, tpoints


def filter_dpaths(dpaths, keyMap, crit):
    """
    OutList = filter_dpaths(files, keyMap, crit).

    Function to narrow down the number of paths in a list based on
    file metadata and selection criteria.

    Inputs:
        dpaths:     list of file paths
        keyMap:     dictionary (xml filenames are keys) which includes
                    the metadata associated with the xml (e.g., ver, cdate,
                    publish, tpoints)
        crit:       string of criteria (e.g., 'tpoints') in which to filter
                    the input list of files

    filter_dpaths will take the top value(s) from the list. For example, if
    the filter criteria is 'publish' and a list with two republished and two
    unpublished files are passed to the function, the two published files
    will be returned. The criteria can be boolean or int (will return True
    booleans and the max integer). Alternatively, if the criteria is creation
    date it will return files with the most recent creation date

    If only one file is in the list, the original list is returned (with
    one file).
    """
    if len(dpaths) < 2:
        return dpaths, keyMap

    # only update dictionary if need be
    dpaths_new = copy.deepcopy(dpaths)
    for dpath in dpaths:
        if crit not in keyMap[dpath].keys():
            cdate, publish, tpoints = get_dataset_metadata(dpath)
            # if the dataset cannot be loaded, remove it
            if 0 in (cdate, publish, tpoints):
                dpaths_new.remove(dpath)
            # update dictionary
            keyMap[dpath]['cdate'] = cdate
            keyMap[dpath]['publish'] = publish
            keyMap[dpath]['tpoints'] = tpoints

    #<qinyi 2021-10-25 #------------------
    if len(dpaths_new) == 0:
        return dpaths_new, keyMap
    #>qinyi 2021-10-25 #------------------

    # get values
    values = []
    for dpath in dpaths_new:
        v = keyMap[dpath][crit]
        values.append(v)

    # determine optimal value
    if type(v) == int:
        vmax = np.max(values)
    else:
        vmax = True
    # create output list
    OutList = []
    for dpath in dpaths_new:
        if keyMap[dpath][crit] == vmax:
            OutList.append(dpath)
    return OutList, keyMap


def versionWeight(V):
    """
    V = versionWeight(ver).

    versionWeight takes a version string for a CMIP xml and returns a
    numeric in which the larger the number, the more recent the version.
    Typically an int corresponding to the date (e.g., 20190829), but will
    give precedence to version numbers (e.g., v1 is returned as 100000000).
    """
    if V == 'latest':
        V = 0
    else:
        V = int(V.replace('v', ''))
        if V < 10:
            V = V*100000000
    V = int(V)
    return V

