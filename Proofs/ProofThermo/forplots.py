import numpy as np
import matplotlib.pylab as plt
datos=np.genfromtxt('datos.txt')
tt=datos[:,0]
Xavg0=datos[:,1]
Yavg0=datos[:,2]
Zavg0=datos[:,3]

# fig=plt.figure(figsize=(10,4.5))
# ax=fig.add_subplot(111)
# ax.plot(t,X)
# plt.show()
f, fig=plt.subplots(3,figsize=(6,4.5),sharex=True)
fig[0].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
fig[1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
fig[2].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
fig[2].set_xlabel(r'$t/\tau$',fontsize=14, fontweight='bold')
fig[0].set_ylabel(r'$x$',rotation=0,fontsize=14, fontweight='bold')
fig[1].set_ylabel(r'$y$',rotation=0,fontsize=14, fontweight='bold')
fig[2].set_ylabel(r'$z$',rotation=0,fontsize=14, fontweight='bold')
fig[0].plot(tt,Xavg0,linestyle='-',color='coral')
fig[1].plot(tt,Yavg0,linestyle='-',color='magenta')
fig[2].plot(tt,Zavg0,linestyle='-',color='navy')
fig[0].set_xlim(tt.min(),tt.max())
fig[1].set_xlim(tt.min(),tt.max())
fig[2].set_xlim(tt.min(),tt.max())
plt.show()
