
def get_lutable(version,exp):
    '''
    ***** Note: please follow this example to add your case info. ******

    'case_shortname': [
    CTL_casename, CTL_startyear, CTL_endyear,
    P4K_casename, P4K_startyear, P4K_endyear,
    ]

    '''

    lu_table = {
    'v1_coupled':[
    'None',None,None,
    'None',None,None,
    ],
    'v2_coupled':[
    'None',None,None,
    'None',None,None,
    ], 
    'v1_amip4K':[
    'None',None,None,
    'None',None,None,
    ], 
    'v1':[
    '20211208.F2010C5-CMIP6-LR.IC.ne30_oECv3.compy.1080',2,3,
    '20211210.F2010C5-CMIP6-LR.IC.p4Ka.ne30_oECv3.compy.1080',2,3,
    ],
    'v2test':[
    '20221013.v2.F2010-CICE.IC.back.clubb.g1.1.ne30pg2_EC.compy',2,3,
    '20221013.v2.F2010-CICE.IC.back.clubb.g1.1.p4Ka.ne30pg2_EC.compy',2,3,
    ],
    'v2test2':[
    '20221013.v2.F2010-CICE.IC.back.clubb.g1.1.ne30pg2_EC.compy2',3,4,
    '20221013.v2.F2010-CICE.IC.back.clubb.g1.1.p4Ka.ne30pg2_EC.compy2',3,4,
    ],
    }

    if exp == 'amip':
        return lu_table[version][0],lu_table[version][1],lu_table[version][2]
    else:
        return lu_table[version][3],lu_table[version][4],lu_table[version][5]


