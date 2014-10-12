__author__ = 'yutongpang'
import numpy as np
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
    modcount = 0
    AMarray = np.array([])
    for line in search:
        if modcount % 27 == 0:
            if modcount > 315 and modcount<2516:
                AMarray = np.append(AMarray,float(line.split()[1]))
        modcount = modcount +1
    for j in range(0,n):
        AM[i,j] = AMarray[j]

#print('Reading.. CIGS-EF-Ga-p3-%dnm-Ex'%(wl))


from pylab import *
x= np.arange(0,81)
ax1 = subplot(111)
ax1.scatter(x,iqewl,marker='o',label='tkhonov',color='black')


show()