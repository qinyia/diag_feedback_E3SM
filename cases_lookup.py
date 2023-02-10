
def get_lutable(version,exp):
    lu_table = {
    'v1':[
    '20211208.F2010C5-CMIP6-LR.IC.ne30_oECv3.compy.1080',
    '20211210.F2010C5-CMIP6-LR.IC.p4Ka.ne30_oECv3.compy.1080',
    ],
    'v2test':[
    '20221013.v2.F2010-CICE.IC.back.clubb.g1.1.ne30pg2_EC.compy',
    '20221013.v2.F2010-CICE.IC.back.clubb.g1.1.p4Ka.ne30pg2_EC.compy',
    ],
    }

    if exp == 'amip':
        return lu_table[version][0]
    else:
        return lu_table[version][1]
