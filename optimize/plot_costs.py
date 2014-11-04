'''
plot the best PSO positions.
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
from tools import get_params

def plot_costs(fnm,show=True,savefig=False):

    params, par_names = get_params()
    print params
    
    data = np.loadtxt(fnm,unpack=True)
    
    mpl.rcParams['figure.figsize'] = 15,12
    mpl.rcParams['font.size'] = 10
    n = len(par_names)
    f, axarr = plt.subplots(n,n)
    j = 0
    for v in par_names:
        i = 0
        for p in par_names:
            if (j < i):
                px = data[i+2,:]
                py = data[j+2,:]
                axarr[i,j].plot(px,py,'ko-',alpha=0.1,mew=0)
                axarr[i,j].plot(px[-1],py[-1],'rs',alpha=0.5,markersize=7,mew=0)
                axarr[i,j].set_xlabel(par_names[i])
                axarr[i,j].set_ylabel(par_names[j])
                axarr[i,j].set_xlim(float(params[i][1]),float(params[i][2]))
                axarr[i,j].set_ylim(float(params[j][1]),float(params[j][2]))
            else:
                f.delaxes(axarr[i,j])
            plt.setp(axarr[i,j].get_xticklabels(), rotation=30)
            i += 1
        j += 1
    f.set_tight_layout(True)
    if show: plt.show()
    if savefig: f.savefig('parameter_matrix.pdf',type='PDF',bbox_inches='tight', pad_inches = 0.0)

if __name__ == "__main__":

    plot_costs(sys.argv[1])

