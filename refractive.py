from scipy.interpolate import interp1d
import numpy as np
#extended refrative index from database
Ag_k = np.genfromtxt('oc_si/Ag-k.txt', delimiter='\t',invalid_raise=False)
ry = np.array([])
for object in Ag_k:
    ry = np.append(ry,object[1])
rx = np.array([])
for object in Ag_k:
    rx = np.append(rx,object[0])
fr = interp1d(rx, ry,kind='cubic')
rxx = np.linspace(300,1100, 81)
ryy = fr(rxx)
print(ryy)
