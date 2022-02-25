#****************************************************************
#
#    Filename: contributed.py
#
#    Author: Yi Qin - qin4@llnl.gov
#    Description: 
#               pivot (flip) the contents of array "x" about some arbitrary
#               user specified longitude. The effect is similar to "lonFlip"
#               However, lonFlip will pivot about the mid point [whatever that is]
#               while thus function allows the user to specify what lon to pivot about.
#               
#               grid must be "global" [no cyclic point] and it assumes that the
#               rightmost dimension of "x" is a coordinate variable corresponding
#               to longitude.
#               usage    xNew = lonPivot (x, 20.)    ; pivot about 20E
#                           x = lonPivot (x, 20.)    ; can overwrite
#               Referred from the contributed.ncl from NCL
#    Input: 
#    Output: 
#    Create: 2021-09-27 18:51:17
#    Last Modified: 2021-09-27 18:51:17
#****************************************************************

import numpy as np
import copy

def find_element_in_list(element, list_element):
    try:
        index_element = list(list_element).index(element)
        return index_element
    except ValueError:
        return None

def lonPivot(x, Lons, pivotLon):
    '''
    flip longitudes
    '''
    sz = x.shape
    nDim = len(sz)

    if nDim > 5:
        print("lonflip: too many dims: nDim=",nDim)
        return x

    temp = copy.deepcopy(x)
    indP = find_element_in_list(pivotLon, Lons)
    if indP == None:
        print('lonPivot: bad pivot value')
        return None

    mlon = len(Lons) # of longitudes 
    indL = mlon - 1  # last index
    n = indL - indP
    #print('indL=',indL,'indP=',indP,'n =',n)

    if nDim == 1:
        temp[0:n+1] = x[indP:indL+1]
        temp[n+1:]  = x[0:indP]
    if nDim == 2:
        temp[:,0:n+1] = x[:,indP:indL+1]
        temp[:,n+1:]  = x[:,0:indP]
    if nDim == 3:
        temp[:,:,0:n+1] = x[:,:,indP:indL+1]
        temp[:,:,n+1:]  = x[:,:,0:indP]
    if nDim == 4:
        temp[:,:,:,0:n+1] = x[:,:,:,indP:indL+1]
        temp[:,:,:,n+1:]  = x[:,:,:,0:indP]
    if nDim == 5:
        temp[:,:,:,:,0:n+1] = x[:,:,:,:,indP:indL+1]
        temp[:,:,:,:,n+1:]  = x[:,:,:,:,0:indP]

    tlon = np.empty((mlon))
    tlon[0:n+1] = Lons[indP:indL+1]
    tlon[n+1:]  = Lons[0:indP]

    if tlon[0] > 0: # (say) 20, 25, ..., 350,355,0,5,...
        indt = [ n for n,m in enumerate(tlon) if m<tlon[0] ] 
        if len(indt) > 0: 
            tlon[indt] = tlon[indt] + 360.

    if tlon[0] >= 180. or tlon[0] == 360:
        tlon = tlon - 360.

    return (temp,tlon)

if __name__ == '__main__':

    x = np.array([[1,2,3,4,5,6,7],[1,2,3,4,5,6,7]])
    Lons = np.array([2.5, 22.5, 117.5, 177.5, 182.5, 347.5, 357.5])
    pivotLon = 182.5
    print('x=',x.shape,x)
    print('Lons=',Lons.shape,Lons)
    print('pivotLon=',pivotLon)

    xx,tlon = lonPivot(x, Lons, pivotLon)
    print('xx = ',xx)
    print('tlon = ',tlon)

