__author__ = 'yutongpang'
import numpy as np
txtarray = np.genfromtxt('oc_si/Si-n.txt',delimiter='\t',invalid_raise=False)
print(txtarray)
reflectionarray = np.array([])
for object in txtarray:
    refractiveindex = object[1]
    reflectionobject  = ((1-refractiveindex)/(1+refractiveindex))**2
    reflectionarray = np.append(reflectionarray,reflectionobject)

qearray = np.genfromtxt('oc_si/si-qe.txt',delimiter='\t',invalid_raise=False)
iqedataarray = np.array([])
i=0
for object in qearray:
    eqe = object[1]
    iqe = object[1]/(1-reflectionarray[i])
    iqedataarray = np.append(iqedataarray,iqe)
    i=i+1
x = np.arange(300,1110,10)

from pylab import *
x = np.arange(300,1110,10)
ax1 = subplot(111)
#ax1.set_yscale('log')
#ax.set_xscale('log')
ax1.scatter (x,iqedataarray,marker='o',label='TSVD Regularization',color='black')

show()