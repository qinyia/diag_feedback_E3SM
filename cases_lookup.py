
def get_lutable(version,exp):
    lu_table = {
    'v1_coupled':[
    'None',
    'None',
    ],
    'v2_coupled':[
    'None',
    'None',
    ],
    'v1':[
    '20211208.F2010C5-CMIP6-LR.IC.ne30_oECv3.compy.1080',
    '20211210.F2010C5-CMIP6-LR.IC.p4Ka.ne30_oECv3.compy.1080',
    ],
    'v2':[
    '20211208.v2.F2010-CICE.IC.ne30pg2_EC.compy',
    '20211210.v2.F2010-CICE.IC.p4Ka.ne30pg2_EC.compy',
    ],
    }

    if exp == 'amip':
        return lu_table[version][0]
    else:
        return lu_table[version][1]
