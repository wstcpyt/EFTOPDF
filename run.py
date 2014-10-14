__author__ = 'yutongpang'
import numpy as np
from scipy import interpolate
from SVD import SVD
iqescr=0.9
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
        cc = np.append(cc,0.9)
    else:
        cc = np.append(cc,collectionfunction(z))

# calculate iqe vs waveklength
wavelengthfile = open('wavelength_normal.txt')
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
print(iqewl)
n = 81
s = (n,n)
#calculate the matrix
AM = np.zeros(s)
for i in range(0,n):
    search = open('pdf/%d_NORMAL_PDF.dat'%(wl[i]))
    AMarray = np.array([])
    AMfoward = np.array([])
    for line in search:
        AMfoward = np.append(AMfoward,float(line.split()[1]))
    AMreversed = AMfoward[::-1]
    modcount = 0
    for object in AMreversed:
        if modcount % 28 == 0:
            if modcount > 315 and modcount<2616:
                AMarray = np.append(AMarray,object)
        modcount = modcount + 1
    for j in range(0,n):
        #AMarray = AMarray[::-1]
        AM[i,j] = AMarray[j]/n

noise = np.random.normal(0,1,(n,n))
#AM = AM +noise/100
#TSVD method
svdclass=SVD(n)
g = iqewl
f_tsvd = svdclass.f_tsvd(4,AM,g.T)
f_tikhonov =svdclass.f_tikhonov(0.08,AM,g.T)
utb, utbs=svdclass.picardparameter(AM,g.T)


from pylab import *
x= np.arange(0,81)

ax1 = subplot(111)
#ax1.set_yscale('log')
#ax.set_xscale('log')
ax1.scatter(x,f_tsvd,marker='o',label='tkhonov',color='black')
#ax1.scatter(x,f_tsvd,marker='o',label='tkhonov',color='green')
ccx = np.arange(0,2300)
ax1.scatter(ccx/2300.0*81,cc,marker='o',label='tkhonov',color='red')

#smoothcurce
import scipy
from scipy.interpolate import interp1d
from scipy import signal
f2 = interp1d(x, f_tsvd)
xx = np.linspace(0,80, 100)
yy = f2(xx)
# make a gaussian window
window = signal.gaussian(7, 20)
smoothed = signal.convolve(yy, window/window.sum(), mode='same')
plt.plot(xx,smoothed)


show()