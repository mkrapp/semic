import numpy as np
import itertools
import random
import math
from tools import run_particles, get_params, init_search


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


def search(search_space, pop_size, prefix):
    f, tmpdir = init_search(prefix,names)
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
    params, names = get_params()
    
    search_space = [[float(x) for x in p[1:]] for p in params]

    search(search_space, 1000, 'topo_')
