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
    f = (-(1150.0-z)**2+1150.0**2)/(1150.0**2)
    return f
#n=81 for pdf
cc =  np.array([])
for z in range(0,2300):
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

#extended IQE from measurement
g= np.array([1.344,11.511,52.300,67.217,88.613,93.901,94.127,94.191,93.928,93.227,92.437,91.066,88.771,84.937,72.927,41.473,9.945])
g=g/100
xiqeextended = np.array([300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1050,1100])
fiqe = interp1d(xiqeextended, g,kind='cubic')
iqexx = np.linspace(300,1100, 81)
iqeyy = fiqe(iqexx)

n = 81
s = (n,n)
#calculate the matrix
AM = np.zeros(s)
print(AM)
noise = np.random.normal(0,1,(n,n))
#AM = AM +noise/100
#TSVD method
svdclass=SVD(n)
g = iqewl
f_tsvd = svdclass.f_tsvd(35,AM,g.T)
f_tikhonov =svdclass.f_tikhonov(0.05,AM,g.T)
utb, utbs=svdclass.picardparameter(AM,g.T)


from pylab import *
x= np.arange(0,81)

ax1 = subplot(111)
#ax1.set_yscale('log')
#ax.set_xscale('log')
ax1.scatter (x/81.0,f_tsvd,marker='o',label='TSVD Regularization',color='black')
#ax1.scatter (x,iqeyy,marker='o',label='tkhonov',color='red')
#ax1.scatter(x,f_tsvd,marker='o',label='tkhonov',color='green')
ccx = np.arange(0,2300)
ax1.plot(ccx/2300.0,cc,label='Exact profiles',color='red',linewidth=3.0)

#smoothcurce
f2 = interp1d(x, f_tsvd)
xx = np.linspace(0,80, 100)
yy = f2(xx)
# make a gaussian window
window = signal.gaussian(10, 20)
smoothed = signal.convolve(yy, window/window.sum(), mode='same')
#plt.plot(xx/80,smoothed,linewidth=3.0,label='Gaussian smooth')
ax1.set_ylim(-0.1,1.7)
ax1.set_xlabel('Normalize position in CIGS layer',fontsize=15)
ax1.set_ylabel('Charge collection probability',fontsize=15)
ax1.xaxis.set_tick_params(labelsize=15)
ax1.yaxis.set_tick_params(labelsize=15)
legend = ax1.legend(loc='upper right', shadow=True,prop={'size':15})
show()