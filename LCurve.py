__author__ = 'yutongpang'
import numpy as np
from SVD import SVD
size = 40
d = 0.25
class Geo:
    def __init__(self, d1, d2, lam):
        self.lam = lam
        s = (size, size)
        #########################################
        self.AM = np.zeros(s)

        for i in range(0, size):
            for j in range(0, size):
                s = 1.0/size * (j + 0.5)
                t = 1.0/size * (i + 0.5)
                self.AM[i][j] = 1.0/size * d/pow(pow(d, 2.0) + pow(s-t, 2.0), 1.5)
        #noise = np.random.normal(0,1,(20,20))
        #AM = AM +noise/5
        ###########################################
        fexact = np.array([])
        for i in range(0, size):
            if i < 8:
                fexact = np.append(fexact, 2.0)
            else:
                fexact = np.append(fexact, 1.0)

        self.b = np.dot(self.AM, fexact)
        ##############################################
        for i in range(0, size):
            for j in range(0, size):
                s = 1.0/size * (j + 0.5)
                t = 1.0/size * (i + 0.5)
                self.AM[i][j] = 1.0/size * d1/pow(pow(d2, 2.0) + pow(s-t, 2.0), 1.5)

    def gettikhonov_x(self):
        from SVD import SVD
        svd = SVD(size)
        f_tikhonov = svd.f_tikhonov(self.lam, self.AM, self.b.T)
        return f_tikhonov


if __name__ == '__main__':
    ro = np.array([])
    yita = np.array([])
    lambsweep = np.arange(0.00001, 0.1, 0.0001)
    for lambobj in lambsweep:
        print(lambobj)
        geo = Geo(0.25, 0.2504, lambobj)
        f_tikhonov = geo.gettikhonov_x()
        from numpy import linalg as LA
        difference = np.dot(geo.AM, f_tikhonov) - geo.b
        rovalue = LA.norm(difference)
        yitavalue = LA.norm(f_tikhonov)
        yita = np.append(yita, yitavalue)
        ro = np.append(ro, rovalue)
    from pylab import *


    ax1 = subplot(111)
    #ax1.set_yscale('log')
    #ax.set_xscale('log')
    ax1.scatter (lambsweep,yita,marker='o',label='TSVD Regularization',color='black')
    #ax1.scatter(x,s,marker=(3, 1),label=r"${\sigma _i}$",color='black')
    show()