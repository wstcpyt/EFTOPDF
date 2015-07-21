__author__ = 'yutongpang'
import numpy as np
from numpy.linalg import inv
import time
from numpy import linalg as LA
from math import sqrt
size = 40
d = 0.25
class TikhonovLeastSquare:
    def __init__(self, A, b, lambdavalue):
        self.A = A
        self.b = b
        self.lambdavalue = lambdavalue

    def calculateX(self):
        matrixsize = 40
        identity = np.identity(matrixsize)
        firstterm = inv(np.dot(self.A.T, self.A) + self.lambdavalue**2*identity)
        secondterm = np.dot(self.A.T, self.b)
        xarray = np.dot(firstterm, secondterm)
        return xarray


class Geo:
    def __init__(self, d1, d2, lam):
        self.lam = lam
        s = (size, size)
        #########################################
        self.AM = np.zeros(s)
        self.AMexact = np.zeros(s)

        for i in range(0, size):
            for j in range(0, size):
                s = 1.0/size * (j + 0.5)
                t = 1.0/size * (i + 0.5)
                self.AM[i][j] = 1.0/size * (d1)/pow(pow(d2, 2.0) + pow(s-t, 2.0), 1.5)
                self.AMexact[i][j] = 1.0/size * d/pow(pow(d, 2.0) + pow(s-t, 2.0), 1.5)
        #noise = np.random.normal(0,1,(20,20))
        #AM = AM +noise/5
        ###########################################
        fexact = np.array([])
        for i in range(0, size):
            if i < 8:
                fexact = np.append(fexact, 2.0)
            else:
                fexact = np.append(fexact, 1.0)

        self.b = np.dot(self.AMexact, fexact)
        ##############################################

    def getleastsquarex(self):
        tikhonovLeastSquare = TikhonovLeastSquare(self.AM, self.b, self.lam)
        xarray = tikhonovLeastSquare.calculateX()
        return xarray

    def gettikhonov_x(self):
            from SVD import SVD
            svd = SVD(size)
            f_tikhonov = svd.f_tikhonov(self.lam, self.AM, self.b.T)
            self.phi = svd.phi
            return f_tikhonov



class Curvature:
    def __init__(self,A, x, z, b, lambdavalue):
        self.A = A
        self.x = x
        self.z = z
        self.b = b
        self.lambdavalue = lambdavalue

    def calculatecValue(self):
        self.calculateyita()
        self.calculatero()
        self.calculateyitaprime()

        firstterm = 2 * self.yita * self.ro / self.yitaprime
        dividend = pow(self.lambdavalue, 2.0) * self.yitaprime * self.ro + 2 * self.lambdavalue * self.yita * self.ro + pow(self.lambdavalue, 4.0) * self.yita * self.yitaprime;
        divisor = (self.lambdavalue**2 * self.yita**2 + self.ro**2)**1.5
        self.cValue = firstterm * dividend / divisor

    def calculateyita(self):
        self.yita = LA.norm(self.x)**2

    def calculatero(self):
        difference = np.dot(self.A, self.x) - self.b
        self.ro = LA.norm(difference)**2

    def calculateyitaprime(self):
        vv = np.dot(self.x.T, self.z)
        self.yitaprime = 4.0/self.lambdavalue * vv


class GSV:
    def __init__(self, A, x, b, phi):
        self.A = A
        self.x = x
        self.b = b
        self.phi = phi

    def calculateGvalue(self):
        arraysize = len(self.b)
        dividend = LA.norm(np.dot(self.A, self.x) - self.b)**2
        sumphi = 0.0
        for i in range(0, arraysize):
            sumphi += self.phi[i]

        divisor = (arraysize - sumphi)**2
        self.Gvalue = dividend/divisor




if __name__ == '__main__':
        # ro = np.array([])
        # yita = np.array([])
        # cvalue = np.array([])
        # lambsweep = np.arange(0.0001, 0.1, 0.0001)
        # for lambobj in lambsweep:
        #     print(lambobj)
        #     geo = Geo(0.25, 0.265, lambobj)
        #     f_tikhonov_X = geo.gettikhonov_x()
        #     bdifference = np.dot(geo.AM, f_tikhonov_X) - geo.b
        #     from SVD import SVD
        #     svd1 = SVD(size)
        #     f_tikhonov_Z = svd1.f_tikhonov(lambobj, geo.AM, bdifference)
        #     curvature = Curvature(geo.AM, f_tikhonov_X, f_tikhonov_Z, geo.b, lambobj)
        #     curvature.calculatecValue()
        #     cvalue = np.append(cvalue, curvature.cValue+1E-5)
        #     yita = np.append(yita, curvature.yita)
        #     ro = np.append(ro, curvature.ro)
        lambsweep = np.arange(0.000001, 0.001, 0.00001)
        bsweep = np.arange(-0.2, 0.2, 0.05)
        bsweepsize = len(bsweep)
        minimumarray = np.zeros((bsweepsize, bsweepsize))
        for i in range(0, bsweepsize):
            print(i)
            for j in range(0, bsweepsize):
                print(j)
                Garray = np.array([])
                for lambsweepobj in lambsweep:
                    geo = Geo(0.25, 0.25, lambsweepobj)
                    geo.b[0] = geo.b[0] + bsweep[i]
                    geo.b[1] = geo.b[1] + bsweep[j]
                    tikhonov_x = geo.gettikhonov_x()
                    leastsquare_x = geo.getleastsquarex()
                    gsv = GSV(geo.AM, leastsquare_x, geo.b, geo.phi)
                    gsv.calculateGvalue()
                    Garray = np.append(Garray, gsv.Gvalue)
                minumumvalue = np.amin(Garray)
                minimumarray[i][j] = minumumvalue

        print(minimumarray)
        from pylab import *
        ax1 = subplot(111)
        x= np.arange(40)
        ax1.set_yscale('log')
        ax1.set_xscale('log')
        #ax1.scatter (lambsweep,Garray,marker='o',label='TSVD Regularization',color='black')
        #ax1.scatter(x,s,marker=(3, 1),label=r"${\sigma _i}$",color='black')
        #show()
