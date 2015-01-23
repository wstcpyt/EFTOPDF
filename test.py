__author__ = 'yutongpang'
import numpy as np
a= np.arange(300,1110,10)
k= len(a)
print(k)

n=np.array([81,41,27,21,17,14,12,11,9,8,7])
k=np.array([18,10,5,4,4,4,4,4,4,4,4])


from pylab import *

ax1 = subplot(111)

#ax.set_xscale('log')
ax1.plot(n,k,marker='o',color='black')
ax1.set_xlabel('N(discretization size)',fontsize=15)
ax1.set_ylabel('K(Turning Points)',fontsize=15)
#ax1.set_title('K=50',fontsize=20)
ax1.xaxis.set_tick_params(labelsize=15)
ax1.yaxis.set_tick_params(labelsize=15)
show()
