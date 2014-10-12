__author__ = 'yutongpang'
import numpy as np
class SVD():
    def __init__(self,n):
        self.n=n
    def kernel(self,s,t):  #expression of the kernel for the first type Fredholm integral equation of the first kind
        d=0.5
        value = d*(d**2+(s-t)**2)**(-3/2)
        return value
    def am(self):  # calculate the matrix expression for the operater
        n=self.n
        j = np.arange(0,n,dtype=np.float)
        sj = (j)/n
        tj = (j)/n
        s = (n,n)
        AM = np.zeros(s)
        for i in range(0,n):
            for j in range(0,n):
                AM[i,j] = self.kernel(sj[i],tj[j])/n
        return AM
    def g(self,AM):  # measured function
        n= self.n
        f = np.zeros((n,))
        for i in range(0,8):
            f[i] = 2
        for i in range (8,n):
            f[i] = 1
        gvalue = np.dot(AM,f)
        return gvalue,f

    def svdmatrix(self,AM):   #svd decomposition
        U, s, V = np.linalg.svd(AM,full_matrices=True)
        return U, s, V
#########################################################
    def picardparameter(self,AM,g): #picadrd parameter
        n = self.n
        utb = np.zeros((n,))
        utbs = np.zeros((n,))
        U, s, V = self.svdmatrix(AM)
        ut = U.T
        for i in range(0,n):
            utb[i] = (np.dot(ut[i],g))
            utbs[i] = (np.dot(ut[i],g))/s[i]

        return utb, utbs
    def lcurve(self,AM,g):    #lcurve
        n=self.n
        utb, utbs = self.picardparameter(AM,g)
        U,s,V=self.svdmatrix(AM)
        lam = np.arange(0.01,10,0.5)
        lcurvex = np.zeros((len(lam),))
        lcurvey = np.zeros((len(lam),))
        for i in range(0,len(lam)):
            lamobject =lam[i]
            lcurveyobject = 0
            for j in range(0,n):
                philam=s[j]**2/(s[j]**2+lamobject**2)
                lcurveyobject=lcurveyobject+(philam*utbs[j])**2
            lcurvey[i]=lcurveyobject
            lcurvexobject = 0
            for j in range(0,n):
                philam=s[j]**2/(s[j]**2+lamobject**2)
                lcurvexobject=lcurvexobject+((1-philam)*utb[j])**2
            identitymatrix =np.matrix(np.identity(n))
            epsi= np.dot((identitymatrix-U*U.T),g)
            epsinorm = np.linalg.norm(epsi)
            lcurvex[i]=np.sqrt(lcurvexobject+epsinorm**2)
        return lcurvex,lcurvey


#######################################################
    def f_tsvd(self,filtern,AM,g):  #f_tsvd regularization
        f_tsvd = 0
        utb, utbs = self.picardparameter(AM,g)
        U,s,V=self.svdmatrix(AM)
        for i in range(0,filtern):
            f_tsvd = f_tsvd + utbs[i]*V[i]
        return f_tsvd

    def f_tikhonov(self,lam,AM,g): #f_tikhonov regularization
        n=self.n
        f_tikhonov=0
        utb, utbs =self.picardparameter(AM,g)
        U,s,V=self.svdmatrix(AM)
        for i in range(0,n):
            f_tikhonov=f_tikhonov+s[i]**2/(s[i]**2+lam**2)*utbs[i]*V[i]
        return f_tikhonov




