'''
plot the best PSO positions.
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.mlab import griddata
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
from tools import get_params


def plot_topo(fnm,fnm2=None,show=True,savefig=False):
    
    params, par_names = get_params()

    data = np.loadtxt(fnm)
    print np.shape(data)
    if fnm2 is not None:
        data2 = np.loadtxt(fnm2)
        print np.shape(data2)
    
    mpl.rcParams['figure.figsize'] = 15,12
    mpl.rcParams['font.size'] = 10
    n = len(par_names)
    f, axarr = plt.subplots(n,n)
    fitness = data[:,0]
    print np.min(fitness), np.max(fitness)
    fitness = np.log(fitness[~np.isnan(fitness)])
    print np.shape(fitness)
    print np.min(fitness), np.max(fitness)
    j = 0
    for v in par_names:
        i = 0
        for p in par_names:
            if (j < i):
                px = data[:,i+1]
                py = data[:,j+1]
                px = px[~np.isnan(fitness)]
                py = py[~np.isnan(fitness)]
                xi = np.linspace(min(px), max(px),100)
                yi = np.linspace(min(py), max(py),100)
                zi = griddata(px, py, fitness, xi, yi)
                p = axarr[i,j].pcolormesh(xi,yi,zi,cmap='Paired')
                if fnm2 is not None:
                    px2 = data2[:,i+2]
                    py2 = data2[:,j+2]
                    axarr[i,j].plot(px2,py2,'ko-',alpha=0.5,mew=0.1)
                    axarr[i,j].plot(px2[-1],py2[-1],'rs',markersize=7)
                axarr[i,j].set_xlabel(par_names[i])
                axarr[i,j].set_ylabel(par_names[j])
                axarr[i,j].set_xlim(float(params[i][1]),float(params[i][2]))
                axarr[i,j].set_ylim(float(params[j][1]),float(params[j][2]))
                if i == 1 :
                    divider = make_axes_locatable(axarr[0,1])
                    cax = divider.append_axes("left", size="10%", pad=0.1)
                    cb = plt.colorbar(p,cax=cax)
                    cb.set_label("Fitness",fontsize=12)
            else:
                f.delaxes(axarr[i,j])
            plt.setp(axarr[i,j].get_xticklabels(), rotation=30)
            i += 1
        j += 1
    f.set_tight_layout(True)
    if show: plt.show()
    if savefig: f.savefig('parameter_matrix.pdf',type='PDF',bbox_inches='tight', pad_inches = 0.0)

if __name__ == "__main__":

    if len(sys.argv) == 2: plot_topo(sys.argv[1])
    if len(sys.argv) == 3: plot_topo(sys.argv[1],sys.argv[2])

