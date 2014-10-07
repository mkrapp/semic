import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys

mpl.rcParams['figure.figsize'] = 6,12

var_names = ['Surface Temperature (K)', 'Albedo', 'Snow Height (m)',
             'SMB (mm/day)', 'Melt (mm/day)', 'Accumulation (mm/day)']

data = np.loadtxt(sys.argv[1])
vali = np.loadtxt(sys.argv[2])

ndata = len(data[0,:])
nsteps = len(data[:,0])

f, ax = plt.subplots(6)
x = np.arange(nsteps)
lw = 2
if nsteps>365: lw=1

for y in range(ndata):
    fac = 1.
    if y >= 3: fac = 86.4e6
    label = None
    if y == 0: label = 'semic'
    ax[y].plot(x,fac*data[:,y],'k-',lw=lw,alpha=0.5,label=label)
    if y == 0: label = sys.argv[2]
    ax[y].plot(x,fac*vali[:,y],'r-',lw=lw,alpha=0.5,label=label)
    ax[y].set_xlim(x[0],x[-1])
    ax[y].grid()
    ax[y].set_title(var_names[y])
    ax[y].set_xlabel('days')
    if y == 0: ax[y].legend(loc=3,fancybox=True, framealpha=0.5,fontsize=12)
plt.tight_layout()
plt.show()

