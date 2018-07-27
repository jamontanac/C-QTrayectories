import numpy as np
import matplotlib.pylab as plt
datos=np.genfromtxt('thermo.txt')
tt=datos[:,0]
Energy=datos[:,1]
Heat=datos[:,2]
Work=datos[:,3]

fig=plt.figure(figsize=(6.5,4.5))
ax=fig.add_subplot(111)
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax.plot(tt,Energy,linestyle='-',color='darkgreen', label='Energy')
ax.plot(tt,Work,linestyle='-',color='navy', label='Work')
ax.plot(tt,Heat,linestyle='-',color='firebrick', label='Heat',alpha=0.9)
ax.legend(loc='upper left')
ax.set_xlabel(r'$t/\tau$',fontsize=14, fontweight='bold')
ax.set_ylabel('Energy',fontsize=14, fontweight='bold')
ax.set_xlim(tt.min(),tt.max())
plt.savefig('Energy.pdf')
plt.show()

