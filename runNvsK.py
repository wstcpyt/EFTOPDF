__author__ = 'yutongpang'
import numpy as np
from scipy.interpolate import interp1d
from scipy import signal
from scipy import interpolate
from SVD import SVD
iqescr=1
zmo=2300
leff=920.0
smooverd=4.3e-4
def collectionfunction(z):
    numerator= iqescr*((np.cosh((z-zmo)/leff))/leff-smooverd*np.sinh((z-zmo)/leff))
    denominator= np.cosh((zmo-350.0)/leff)/leff+smooverd*np.sinh((zmo-350.0)/leff)
    iqevalue=numerator/denominator
    return iqevalue
#n=81 for pdf
cc =  np.array([])
for z in range(0,2300):
    if z < 351:
        cc = np.append(cc,1)
    else:
        cc = np.append(cc,collectionfunction(z))

# calculate iqe vs waveklength
wl = np.arange(300,1110,120)
iqewl = np.array([])
for wlobject in wl:
    file = open('pdf/%d_NORMAL_PDF.dat'%(wlobject))
    gf = np.array([])
    gffoward = np.array([])
    for line in file:
        gffoward = np.append(gffoward,float(line.split()[1]))
    gfreversed = gffoward[::-1]
    i=0
    for g in gfreversed:
        if i>315 and i<2616:
            gf = np.append(gf,g)
        i = i+1
    iqewl = np.append(iqewl,np.dot(gf,cc)/2300)


n = len(wl)
print (n)
discreteno = (2616.0-315.0)/n
s = (n,n)
#calculate the matrix
AM = np.zeros(s)
for i in range(0,n):
    search = open('pdf/%d_NORMAL_PDF.dat'%(wl[i]))
    rawAMarray = np.array([])
    AMarray = np.array([])
    AMfoward = np.array([])
    for line in search:
        AMfoward = np.append(AMfoward,float(line.split()[1]))
    AMreversed = AMfoward[::-1]
    modcount = 0
    for object in AMreversed:
        if modcount > 315 and modcount<2616:
            rawAMarray = np.append(rawAMarray,object)
        modcount = modcount + 1
    for k in range(0,n):
        rawcount = 0
        sumAM = 0.0
        for rawobject in rawAMarray:
            if rawcount > discreteno*k and rawcount < discreteno*(k+1):
                sumAM = sumAM + rawobject
            rawcount = rawcount+1
        sumAM = sumAM/discreteno
        AMarray = np.append(AMarray,sumAM)
    for j in range(0,n):
        #AMarray = AMarray[::-1]
        AM[i,j] = AMarray[j]/n
noise = np.random.normal(0,1,(n,n))
#AM = AM +noise/100
#TSVD method
svdclass=SVD(n)
g = iqewl
#f_tsvd = svdclass.f_tsvd(40,AM,g.T)
#f_tikhonov =svdclass.f_tikhonov(0.05,AM,g.T)
utb, utbs=svdclass.picardparameter(AM,g.T)
U, s, V =svdclass.svdmatrix(AM)

from pylab import *
x= np.arange(0,n)

ax1 = subplot(111)
ax1.set_yscale('log')
#ax.set_xscale('log')
ax1.scatter (x,abs(utbs),marker='o',label='TSVD Regularization',color='black')

ax1.set_xlabel('Normalize position in CIGS layer',fontsize=15)
ax1.set_ylabel('Charge collection probability',fontsize=15)
ax1.set_title('K=50',fontsize=20)
ax1.xaxis.set_tick_params(labelsize=15)
ax1.yaxis.set_tick_params(labelsize=15)
legend = ax1.legend(loc='upper right', shadow=True,prop={'size':15})
show()