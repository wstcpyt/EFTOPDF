__author__ = 'yutongpang'
import numpy as np
from SVD import SVD
size = 20
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

    def getfitness(self):
        from SVD import SVD
        svd = SVD(size)
        f_tikhonov = svd.f_tikhonov(self.lam, self.AM, self.b.T)
        fvalue = 0.0
        for f_tikobj in f_tikhonov:
            fvalue += f_tikobj * f_tikobj
        return np.std(f_tikhonov, dtype=np.float64)/np.mean(f_tikhonov)



def get1dplot():
    from pylab import *
    ax1 = subplot(111)
    ax1.set_yscale('log')
    lambsweep = np.array([10 ,1 , 0.1, 0.0001, 0.0000001])
    for lambobj in lambsweep:
        sweep = np.arange(50, 300, 1)/1000.0
        #print(sweep)
        plotsubarray = np.array([])
        for obj in sweep:
            geo = Geo(obj, obj, lambobj)
            plotsubarray = np.append(plotsubarray, geo.getfitness())
        #print(plotsubarray)
        ax1.plot (sweep,plotsubarray,label=str(lambobj))

    legend = ax1.legend(loc='upper left', shadow=True,prop={'size':15})
    show()

def get2dplot():
    d1 = np.arange(245, 255, 1)/1000.0
    d2 = np.arange(245, 255, 1)/1000.0
    colormap = np.zeros((len(d1), len(d1)))
    for i in range(0, len(d1)):
        for j in range(0, len(d2)):
            geo = Geo(d1[i], d2[j], 0.000000001)
            colormap[i][j] = geo.getfitness()
    print(colormap)
    import matplotlib.pyplot as plt
    plt.imshow(colormap)
    plt.colorbar(orientation='vertical')
    plt.show()