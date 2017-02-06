import os.path
from subprocess import call
import tempfile as tmp
import numpy as np
import sys

def write_namelist(fname,names,particle):
    #print "write parameters to FORTRAN namelist " + fname
    f=open(fname,'w')
    f.write('&surface_physics\n')
    f.write('  boundary = "", "", "",\n')
    f.write('  tstic = 86400.,\n')
    f.write('  ceff = 2.0e6,\n')
    f.write('  csh = 2.0e-3,\n')
    f.write('  clh = 5.0e-4,\n')
#    f.write('  albi = 0.45,\n')
#    f.write('  albl = 0.15,\n')
#    f.write('  alb_smin = 0.79,\n')
#    f.write('  alb_smax = 0.91,\n')
#    f.write('  hcrit = 0.5,\n')
#    f.write('  amp = 2.4,\n')
#    f.write('  rcrit = 0.5,\n')
    f.write('  tmin = 263.15,\n')
    f.write('  tmax = 273.15,\n')
#    f.write('  tau_a = 0.008,\n')
#    f.write('  tau_f = 0.24,\n')
#    f.write('  w_crit = 5.0,\n')
#    f.write('  mcrit = 6.0e-8,\n')
    f.write('  afac  = -0.18,\n')
    f.write('  tmid  = 275.35,\n')
    f.write('  n_ksub = 3,\n')
    # loop through parameter values
    for j in range(0,len(names)):
    	f.write('  ' + names[j] + ' = %.8f' % particle["pos"][j] + ',\n')
    f.write('  alb_scheme = "none",\n')
    f.write('/\n')
    f.write('&smb_output\n')
    f.write('  file_timser  = "",\n')
    f.write('  file_daily   = "",\n')
    f.write('  file_diag    = "",\n')
    f.write('  file_monthly = "",\n')
    f.write('/\n')
    f.close()


def run_particles(pop, nml_prefix):
    exe = 'run_particles.x'
    if not os.path.isfile(exe):
        print "\033[91mFile "+exe+" not found.\n Run:\033[0m make "+exe
        sys.exit()
    cost = {}
    for p in pop:
        nml = nml_prefix+"%06d.nml" % p["id"]
        write_namelist(nml,names,p)
    i0 = 0
    i1 = len(pop)-1
    #call(["./"+exe,'"'+nml_prefix+'"',str(i0),str(i1)])
    call(["mpiexec","-n","2","./"+exe,'"'+nml_prefix+'"',str(i0),str(i1)])    
    for p in pop:
        output = nml_prefix+"%06d.out" % p["id"]
        cost[p["id"]] = np.loadtxt(output)
        if np.isnan(cost[p["id"]]): cost[p["id"]] = np.finfo('d').max
    return cost

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


def init_search(prefix,names):
    tmpdir = tmp.gettempdir()
    f = open(prefix+'cost.txt','w')
    f.write("# gen fitness ")
    for n in names:
        f.write(n+" ")
    f.write("\n")
    return f, tmpdir


def get_params():
    #Get information on parameters and ranges from user file
    params = rdcsv('namelist.ranges')


    #format is:
    #param-name,low,high
    
    #set dimension
    k=len(params)
    
    # get names of parameters
    names=[]
    for i in range(0,k):
        names.append(params[i][0])

    return params, names

params, names = get_params()
