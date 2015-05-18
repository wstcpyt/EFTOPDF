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
wavelengthfile = open('wavelength_normal_w.txt')
wl = np.array([])
for line in wavelengthfile:
    wl = np.append(wl,float(line))
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

#extended IQE from measurement
g= np.array([1.344,11.511,52.300,67.217,88.613,93.901,94.127,94.191,93.928,93.227,92.437,91.066,88.771,84.937,72.927,41.473,9.945])
g=g/100
xiqeextended = np.array([300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1050,1100])
fiqe = interp1d(xiqeextended, g,kind='cubic')
iqexx = np.linspace(300,1100, 81)
iqeyy = fiqe(iqexx)

n = 17
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
            if rawcount > 135*k and rawcount < 135*(k+1):
                sumAM = sumAM + rawobject
            rawcount = rawcount+1
        sumAM = sumAM/135.0
        AMarray = np.append(AMarray,sumAM)
    for j in range(0,n):
        #AMarray = AMarray[::-1]
        AM[i,j] = AMarray[j]/n
print(AM)
#noise = np.random.normal(0,1,(n,n))
#AM = AM +noise/1000
#TSVD method
svdclass=SVD(n)
g = iqewl
f_tsvd = svdclass.f_tsvd(3,AM,g.T)
f_tikhonov =svdclass.f_tikhonov(0.05,AM,g.T)
utb, utbs=svdclass.picardparameter(AM,g.T)
U, s, V =svdclass.svdmatrix(AM)

from pylab import *
x= np.arange(0,17)

ax1 = subplot(111)
#ax1.set_yscale('log')
#ax.set_xscale('log')
ax1.scatter (x/17.0,f_tsvd,marker='o',label='TSVD Regularization',color='black')
#ax1.scatter (x,s,marker='o',label=r"${\sigma _i}$",color='black')
#ax1.scatter (x,abs(utb),marker='o',label='$ {u_i^TP}$',color='red')
#ax1.scatter(x,abs(utbs),marker='o',label='${u_i^TP}/{\sigma _i}$',color='green')
ccx = np.arange(0,2300)
#ax1.plot(ccx/2300.0,cc,label='Exact profiles',color='red',linewidth=3.0)

#smoothcurce
f2 = interp1d(x, f_tsvd)
xx = np.linspace(0,16, 100)
yy = f2(xx)
# make a gaussian window
window = signal.gaussian(10, 20)
smoothed = signal.convolve(yy, window/window.sum(), mode='same')
#plt.plot(xx/80,smoothed,linewidth=3.0,label='Gaussian smooth')
#ax1.set_ylim(-0.1,1.7)
ax1.set_xlabel('Index for SVD',fontsize=15)
ax1.set_ylabel('Value (dimensionless)',fontsize=15)
#ax1.set_title('K=17',fontsize=20)
ax1.xaxis.set_tick_params(labelsize=15)
ax1.yaxis.set_tick_params(labelsize=15)
legend = ax1.legend(loc='upper right', shadow=True,prop={'size':15})
show()