'''
plot the best PSO positions.
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.mlab import griddata
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys

def rdcsv(name):
    f = open(name,"r")
    data = []
    tmp = f.readline()
    while tmp !="":
        onerow = tmp.split(",")
        last = len(onerow)-1
        k = len(onerow[last])
        if k!= 1:
           onerow[last] = onerow[last][0:k-1]
        if len(onerow) !=  0:
           data.append(onerow)
        tmp = f.readline()
    f.close()
    return(data)

def plot_topo(fnm,show=True,savefig=False):
    params = rdcsv('topology.ranges')
    print params
    
    par_names=[]
    for i in range(0,len(params)):
        par_names.append(params[i][0])
    
    data = np.loadtxt(fnm)
    print np.shape(data)
    
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
                xi = np.linspace(min(px), max(px),25)
                yi = np.linspace(min(py), max(py),25)
                zi = griddata(px, py, fitness, xi, yi,interp='linear')
                p = axarr[i,j].pcolorfast(xi,yi,zi,cmap='jet')
                axarr[i,j].set_xlabel(par_names[i])
                axarr[i,j].set_ylabel(par_names[j])
                axarr[i,j].set_xlim(float(params[i][1]),float(params[i][2]))
                axarr[i,j].set_ylim(float(params[j][1]),float(params[j][2]))
                if i == 1 :
                    divider = make_axes_locatable(axarr[3,4])
                    cax = divider.append_axes("left", size="10%", pad=0.1)
                    cb = plt.colorbar(p,cax=cax)
                    cb.set_label("Fitness",fontsize=12)
            else:
                f.delaxes(axarr[i,j])
            plt.setp(axarr[i,j].get_xticklabels(), rotation=30)
            i += 1
        j += 1
    #f.subplots_adjust(right=0.8)
    #cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
    #f.colorbar(p, cax=cbar_ax)
    f.set_tight_layout(True)
    if show: plt.show()
    if savefig: f.savefig('parameter_matrix.pdf',type='PDF',bbox_inches='tight', pad_inches = 0.0)

if __name__ == "__main__":

    plot_topo(sys.argv[1])

