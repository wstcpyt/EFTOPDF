from scipy.interpolate import interp1d
import numpy as np
#extended refrative index from database
Ag_k = np.genfromtxt('oc_si/asi-n.txt', delimiter='\t',invalid_raise=False)
ry = np.array([])
for object in Ag_k:
    ry = np.append(ry,object[1])
rx = np.array([])
for object in Ag_k:
    rx = np.append(rx,object[0])
print(ry)
fr = interp1d(rx, ry)
rxx = np.linspace(300,1100, 81)
ryy = fr(rxx)
np.savetxt('newx.txt',rxx)
np.savetxt('newy.txt',ryy)
print(ryy)
