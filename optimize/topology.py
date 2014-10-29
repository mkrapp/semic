import numpy as np
from subprocess import call
import tempfile as tmp
import sys
import os.path
import itertools
import random
import math

import plot_topo as pt


def write_namelist(fname,particle):
    #print "write parameters to FORTRAN namelist " + fname
    f=open(fname,'w')
    f.write('&surface_physics\n')
    f.write('  boundary = "", "", "",\n')
    f.write('  tstic = 86400.,\n')
    f.write('  ceff = 2.0e6,\n')
    f.write('  csh = 1.5e-3,\n')
    f.write('  clh = 6.0e-4,\n')
    f.write('  albl = -999,\n')
#    f.write('  albr = -999,\n')
#    f.write('  alb_smin = -999,\n')
#    f.write('  alb_smax = -999,\n')
    f.write('  tmin = 263.15,\n')
    # loop through parameter values
    for j in range(0,k):
    	f.write('  ' + names[j] + ' = %.8f' % particle["pos"][j] + ',\n')
    f.write('  method = "ebm",\n')
    f.write('  alb_scheme = "slater",\n')
    f.write('/\n')
    f.write('&smb_output\n')
    f.write('  file_timser  = "",\n')
    f.write('  file_daily   = "",\n')
    f.write('  file_diag    = "",\n')
    f.write('  file_monthly = "",\n')
    f.write('/\n')
    f.close()


def create_search_space(search_space, inc):
    plen = len(search_space)
    slist = []
    for i, p in enumerate(search_space):
        plist = []
        steps = int((search_space[i][1]-search_space[i][0])/inc[i])+1
        for n in range(steps):
            plist.append(search_space[i][0] + n*inc[i])
        slist.append(plist)
    return list(itertools.product(*slist))

def random_vector(minmax):
    return [ minmax[i][0] + ((minmax[i][1] - minmax[i][0]) * random.random()) for i in range(len(minmax))]

def init_population_random(search_space,n):
    population = [{} for _ in range(n)]
    for i,p in enumerate(population):
        p["id"] = i
        p["pos"] = random_vector(search_space)
    return population

def init_population_lhs(search_space,n):
    k = len(search_space)
    pranges=[]
    for i in range(0,k):
        tmp=[]
        low=float(search_space[i][0])
        high=float(search_space[i][1])
        delta=(high-low)/float(n)
        for j in range(0,n):
            tmp.append(random.uniform(low+j*delta,low+(j+1)*delta))
        pranges.append(tmp)
         
    if (n<math.ceil(4.0/3.0*k)):
        print "n is too small"
    
    s=[]
    for i in range(0,k):
        s.append(range(0,n))
    
    result=[]
    for i in range(0,n):
        tmp=[]
        for j in range(0,k):
            a = random.sample(s[j],1)[0]
            tmp.append(a)
            s[j].remove(a)
        result.append(tmp)

    sample = []
    for l in range(0,len(result)):
        tmp=[]
        for j in range(0,len(result[l])):
            tmp.append(pranges[j][result[l][j]])
        sample.append(tmp)

    population = [{} for _ in range(n)]
    for i,p in enumerate(population):
        p["id"] = i
        p["pos"] = sample[i]
    return population


def run_particles(pop, nml_prefix):
    exe = 'run_particles.x'
    if not os.path.isfile(exe):
        print "\033[91mFile "+exe+" not found.\n Run:\033[0m make "+exe
        sys.exit()
    cost = {}
    print 'write namelist files'
    for p in pop:
        nml = nml_prefix+"%06d.nml" % p["id"]
        write_namelist(nml,p)
    i0 = 0
    i1 = len(pop)-1
    #call(["./"+exe,'"'+nml_prefix+'"',str(i0),str(i1)])
    call(["mpiexec","-n","3","./"+exe,'"'+nml_prefix+'"',str(i0),str(i1)])
    for p in pop:
        output = nml_prefix+"%06d.out" % p["id"]
        cost[p["id"]] = np.loadtxt(output)
    return cost

def search(search_space, pop_size, prefix):
    tmpdir = tmp.gettempdir()
    f = open(prefix+'cost.txt','w')
    f.write("#fitness ")
    for n in names:
        f.write(n+" ")
    f.write("\n")
    #pop = init_population_random(search_space, pop_size)
    pop = init_population_lhs(search_space, pop_size)
    try:
        cost = run_particles(pop,tmpdir+'/'+prefix)
        for i,candidate in enumerate(pop):
            candidate["cost"] = cost[candidate["id"]]
            f.write("%.6g " % candidate["cost"])
            [f.write("%.6g "%p) for p in candidate["pos"]]
            f.write("\n")
            #print candidate["pos"]
    except KeyboardInterrupt:
        print "\033[91mInterruption by user. Exiting...\033[0m"
    	pass
    f.close()


if __name__ == "__main__":
    
    #Get information on parameters and ranges from user file
    params = pt.rdcsv('namelist.ranges')
    #format is:
    #param-name,low,high
    
    #set dimension
    k=len(params)
    
    # get names of parameters
    names=[]
    for i in range(0,k):
        names.append(params[i][0])

    search_space = [[float(x) for x in p[1:]] for p in params]

    search(search_space, 10000, 'topo_')

    #pt.plot_topo('topo_cost.txt',show=True,savefig=True)
