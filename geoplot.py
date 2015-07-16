__author__ = 'yutongpang'
import numpy as np
from SVD import SVD
size = 20
d = 0.25

optimumarray = np.array([])

for loop in range(0, 1):
    s = (size, size)
    #derror = d + d/loop

    #########################################
    AM = np.zeros(s)

    for i in range(0, size):
        for j in range(0, size):
            s = 1.0/size * (j + 0.5)
            t = 1.0/size * (i + 0.5)
            AM[i][j] = 1.0/size * d/pow(pow(d, 2.0) + pow(s-t, 2.0), 1.5)
    #noise = np.random.normal(0,1,(20,20))
    #AM = AM +noise/5
    ###########################################
    fexact = np.array([])
    for i in range(0, size):
        if i < 8:
            fexact = np.append(fexact, 2.0)
        else:
            fexact = np.append(fexact, 1.0)

    b = np.dot(AM, fexact)
    ##############################################
    for i in range(0, size):
        for j in range(0, size):
            s = 1.0/size * (j + 0.5)
            t = 1.0/size * (i + 0.5)
            AM[i][j] = 1.0/size * 0.27/pow(pow(0.25, 2.0) + pow(s-t, 2.0), 1.5)

    from SVD import SVD
    svd = SVD(size)
    f_tikhonov = svd.f_tikhonov(0.00001, AM, b.T)
    fvalue = 0.0
    for f_tikobj in f_tikhonov:
        fvalue += f_tikobj * f_tikobj

    fvalue1 = 0.0
    bre = np.dot(AM, f_tikhonov)
    bdif = bre - b
    for bdifobj in bdif:
        fvalue1 += bdifobj*bdifobj

    optimumarray = np.append(optimumarray, fvalue)



print(optimumarray)

from pylab import *
x= np.arange(0,20)

ax1 = subplot(111)
#ax.set_xscale('log')
ax1.scatter (x/20.0,f_tikhonov,marker='o',label='TSVD Regularization',color='black')
show()
