__author__ = 'yutongpang'
import numpy as np
txtarray = np.genfromtxt('oc_si/asi-n.txt',delimiter='\t',invalid_raise=False)
print(txtarray)
reflectionarray = np.array([])
for object in txtarray:
    refractiveindex = object[1]
    reflectionobject  = ((1-refractiveindex)/(1+refractiveindex))**2
    reflectionarray = np.append(reflectionarray,reflectionobject)

fdtdref = np.genfromtxt('oc_si/reflection.txt',delimiter = '\t', invalid_raise=False)
fdtdrefarray = np.array([])
for object in fdtdref:
    reflectionobject = object[1]
    fdtdrefarray = np.append(fdtdrefarray,reflectionobject)

qearray = np.genfromtxt('oc_si/si-qe.txt',delimiter='\t',invalid_raise=False)
fdtdrefarray = np.genfromtxt('si_florida/RF.txt',delimiter='\t',invalid_raise=False)
iqedataarray = np.array([])
i=0
for object in qearray:
    eqe = object[1]
    iqe = object[1]/(1-fdtdrefarray[i]/100.0)
    iqedataarray = np.append(iqedataarray,iqe)
    i=i+1
rxx = np.arange(300,1110,10)
ryy = iqedataarray
np.savetxt('x.txt',rxx)
np.savetxt('y.txt',ryy)

from pylab import *
x = np.arange(300,1110,10)
ax1 = subplot(111)
#ax1.set_yscale('log')
#ax.set_xscale('log')
ax1.scatter (x,iqedataarray,marker='o',label='TSVD Regularization',color='black')

show()