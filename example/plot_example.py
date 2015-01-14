import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats.stats import pearsonr
import sys

mpl.rcParams['figure.figsize'] = 8,12

var_names = ['Surface Temperature (K)', 'SW net (W/m2)', 'SMB (mm/day)',
             'Melt (mm/day)', 'Accumulation (mm/day)', 'SHF (W/m2)', 'LHF (W/m2)']

data_in = np.loadtxt(sys.argv[1])
vali_in = np.loadtxt(sys.argv[2])
bstart = int(sys.argv[3])


nbasins = len(data_in[0,:])/8
ndata = 8
nsteps = len(data_in[:,0])

if bstart >= nbasins:
    print 'basin number too large. Must be smaller than %i' % nbasins
    sys.exit()

data = data_in[:,bstart::nbasins]
vali = vali_in[:,bstart::nbasins]

# remove albedo (2nd entry)
data = np.delete(data,1,1)
vali = np.delete(vali,1,1)

f, ax = plt.subplots(7,4)
#f = plt.figure(figsize=(6,12))

x = np.arange(nsteps)
lw = 2
if nsteps>365: lw=1

for y in range(ndata-1):
    ax[y,0] = plt.subplot2grid((7,4), (y,0), colspan=3)
    fac = 1.
    if y >= 2 and y < 5: fac = 86.4e6
    label = None
    if y == 1: label = 'semic'
    ax[y,0].plot(x,fac*data[:,y],'k-',lw=lw,alpha=0.5,label=label)
    if y == 1: label = sys.argv[2]
    ax[y,0].plot(x,fac*vali[:,y],'r-',lw=lw,alpha=0.5,label=label)
    ax[y,0].set_xlim(x[0],x[-1])
    ax[y,0].grid()
    ax[y,0].set_title(var_names[y])
    ax[y,0].set_xlabel('days')
    if y == 1: ax[y,0].legend(loc=2,fancybox=True, framealpha=0.5,fontsize=12)
    minmax_list = [(fac*data[:,y]).tolist(),(fac*vali[:,y]).tolist()]
    ymin = np.min(minmax_list)
    ymax = np.max(minmax_list)
    ax[y,-1].plot(np.linspace(ymin,ymax),np.linspace(ymin,ymax),'k-')
    corcoef = pearsonr(fac*data[:,y],fac*vali[:,y])[0]
    x1 = fac*data[:,y]
    x1_avg = np.average(fac*data[:,y])
    x2 = fac*vali[:,y]
    x2_avg = np.average(fac*vali[:,y])
    crmse = np.sqrt(np.sum(((x1 - x1_avg)-(x2 - x2_avg))**2)/len(vali[:,y]))
    ax[y,-1].plot(fac*data[:,y],fac*vali[:,y],ls='',marker='o',ms=3,color='k',alpha=0.2,mew=0.1)
    ax[y,-1].text(0.05,0.80,'%4.2f'%crmse,transform=ax[y,-1].transAxes,fontsize=10,
                  bbox=dict(fc='lightgray'),alpha=0.5)
    ax[y,-1].set_ylim(ymin,ymax)
    ax[y,-1].set_xlim(ymin,ymax)
    ax[y,-1].set_xlabel('semic')
    ax[y,-1].set_ylabel('data')
f.set_tight_layout(True)
plt.show()

