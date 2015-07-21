__author__ = 'yutongpang'
import numpy as np
from numpy.linalg import inv
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

    def getleastsquarex(self):
        tikhonovLeastSquare = TikhonovLeastSquare(self.AM, self.b, self.lam)
        xarray = tikhonovLeastSquare.calculateX()
        return xarray


if __name__ == '__main__':
        geo = Geo(0.25, 0.25, 0.0000001)
        f_tikhonov = geo.getleastsquarex()
        from numpy import linalg as LA
        difference = np.dot(geo.AM, f_tikhonov) - geo.b
        rovalue = LA.norm(difference)
        yitavalue = LA.norm(f_tikhonov)
        print(rovalue)
        print(yitavalue)


