import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats.stats import pearsonr
import sys

mpl.rcParams['figure.figsize'] = 8,12

var_names = ['Surface Temperature (K)', 'Albedo', 'Snow Height (m)',
             'SMB (mm/day)', 'Melt (mm/day)', 'Accumulation (mm/day)']

data = np.loadtxt(sys.argv[1])
vali = np.loadtxt(sys.argv[2])

ndata = len(data[0,:])
nsteps = len(data[:,0])

f, ax = plt.subplots(6,4)
#f = plt.figure(figsize=(6,12))

#ax2 = plt.subplot2grid((4,4), (0,1), colspan=3)
x = np.arange(nsteps)
lw = 2
if nsteps>365: lw=1

for y in range(ndata):
    ax[y,0] = plt.subplot2grid((6,4), (y,0), colspan=3)
    fac = 1.
    if y >= 3: fac = 86.4e6
    label = None
    if y == 0: label = 'semic'
    ax[y,0].plot(x,fac*data[:,y],'k-',lw=lw,alpha=0.5,label=label)
    if y == 0: label = sys.argv[2]
    ax[y,0].plot(x,fac*vali[:,y],'r-',lw=lw,alpha=0.5,label=label)
    ax[y,0].set_xlim(x[0],x[-1])
    ax[y,0].grid()
    ax[y,0].set_title(var_names[y])
    ax[y,0].set_xlabel('days')
    if y == 0: ax[y,0].legend(loc=3,fancybox=True, framealpha=0.5,fontsize=12)
    minmax_list = [(fac*data[:,y]).tolist(),(fac*vali[:,y]).tolist()]
    ymin = np.min(minmax_list)
    ymax = np.max(minmax_list)
    ax[y,-1].plot(np.linspace(ymin,ymax),np.linspace(ymin,ymax),'k-')
    corcoef = pearsonr(fac*data[:,y],fac*vali[:,y])[0]
    ax[y,-1].plot(fac*data[:,y],fac*vali[:,y],ls='',marker='o',ms=3,color='k',alpha=0.2,mew=0)
    ax[y,-1].text(0.05,0.80,'%4.2f'%corcoef,transform=ax[y,-1].transAxes,fontsize=10,
                  bbox=dict(fc='lightgray'),alpha=0.5)
    ax[y,-1].set_ylim(ymin,ymax)
    ax[y,-1].set_xlim(ymin,ymax)
    ax[y,-1].set_xlabel('semic')
    ax[y,-1].set_ylabel(sys.argv[2])
f.set_tight_layout(True)
plt.show()

