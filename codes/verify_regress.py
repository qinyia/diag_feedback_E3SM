
import numpy as np
import genutil
import cdms2 as cdms
from PlotDefinedFunction import linearregression_nd

y = np.array([[1,1,1],
[2,2,2],
[3,3,3],
[4,4,4],
[5,5,5],
[6,6,6],
[7,7,8],
[8,8,7],])

x = np.array([[1],[2],[3],[4],[5],[6],[7],[8]])
print(x.shape, y.shape)

error = 2
probability=1

V, E, Pt1, Pt2, Pf1, Pf2 = linearregression_nd(y,x,error=error, probability=probability,
                       noslope=None, nointercept=None)

print('result1=', V, E, Pt1, Pt2, Pf1, Pf2)

y1 = cdms.asVariable(y)
x1 = cdms.asVariable(np.reshape(x,(x.shape[0])))

newtime = cdms.createAxis(range(x.shape[0]))
newtime.id = "time" # name of dimension
newtime.designateTime()  # tell cdms to add attributes that make it time
newtime.units = "months since 2017"
y1.setAxis(0,newtime)
x1.setAxis(0,newtime)

V_new, E_new, Pt1_new, Pt2_new, Pf1_new, Pf2_new = genutil.statistics.linearregression(y1,x=x1,error=error, probability=probability, noslope=None, nointercept=None)
print('result2=', V, E, Pt1, Pt2, Pf1, Pf2)

print(np.array(V)-np.array(V_new), np.array(E)-np.array(E_new), np.array(Pt1)-np.array(Pt1_new), np.array(Pt2)-np.array(Pt2_new), np.array(Pf1)-np.array(Pf1_new), np.array(Pf2)-np.array(Pf2_new))


