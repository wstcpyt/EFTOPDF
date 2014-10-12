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
#n=82 for pdf
cc =  np.array([])
for z in range(0,2300):
    if z < 351:
        cc = np.append(cc,0.9)
    else:
        cc = np.append(cc,collectionfunction(z))
file = open('pdf/900_NORMAL_PDF.dat')
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
sum = np.dot(gf,cc)

print(sum/2300)

from pylab import *
x= np.arange(0,2300)
ax1 = subplot(111)
ax1.scatter(x,gf,marker='o',label='tkhonov',color='black')


show()