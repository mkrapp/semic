import numpy as np
from subprocess import call
import tempfile as tmp
import sys
import os.path
import argparse

import particle_swarm_optimization as pso
import evolution_strategies as es
import cultural_algorithm as ca
import random_search as rs
import plot_costs as pc


def write_namelist(fname,particle):
    #print "write parameters to FORTRAN namelist " + fname
    f=open(fname,'w')
    f.write('&surface_physics\n')
    f.write('  boundary = "", "", "",\n')
    f.write('  tstic = 86400.,\n')
    f.write('  ceff = 2.0e6,\n')
    f.write('  csh = 1.5e-3,\n')
    f.write('  clh = 3.0e-4,\n')
    f.write('  albl = -999,\n')
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

def generate_vel_space(search_space,frac=1./5.):
    vel_space = []
    for i,s in enumerate(search_space):
    	vel_space.append([-frac*abs(s[1]-s[0])/2.,frac*abs(s[1]-s[0])/2.])
    return vel_space

def run_particles(pop, nml_prefix):
    exe = 'run_particles.x'
    if not os.path.isfile(exe):
        print "\033[91mFile "+exe+" not found.\n Run:\033[0m make "+exe
        sys.exit()
    cost = {}
    for p in pop:
        nml = nml_prefix+"%03d.nml" % p["id"]
        write_namelist(nml,p)
    i0 = 0
    i1 = len(pop)-1
    call(["./"+exe,'"'+nml_prefix+'"',str(i0),str(i1)])
    for p in pop:
        output = nml_prefix+"%03d.out" % p["id"]
        cost[p["id"]] = np.loadtxt(output)
    return cost

def search_rs(search_space, max_iter, prefix, inc=None):
    tmpdir = tmp.gettempdir()
    f = open(prefix+'cost.txt','w')
    f.write("# iter fitness ")
    for n in names:
        f.write(n+" ")
    f.write("\n")
    best = None
    try:
        for i in range(max_iter):
            candidate = {}
            candidate["id"] = 0
            candidate["pos"] = rs.random_vector(search_space, inc)
            cost = run_particles([candidate],tmpdir+'/'+prefix)
            candidate["cost"] = cost[candidate["id"]]
            if best is None or candidate["cost"] < best["cost"]: best = candidate
            f.write("%i %.6g " % (i, best["cost"]))
            [f.write("%.6g "%p) for p in best["pos"]]
            f.write("\n")
            print " > iteration %i, best=%.4g" % (i, best["cost"])
            print names
            print best["pos"]
    except KeyboardInterrupt:
        print "\033[91mInterruption by user. Exiting...\033[0m"
    	pass
    f.close()
    return best


def search_pso(max_gens, search_space, vel_space, pop_size, c1, c2, prefix):
    tmpdir = tmp.gettempdir()
    f = open(prefix+'cost.txt','w')
    f.write("# gen fitness ")
    for n in names:
        f.write(n+" ")
    f.write("\n")
    pop = [pso.create_particle(search_space, vel_space, id=i) for i in range(pop_size)]
    gbest = pso.get_global_best(pop)
    try:
    	for gen in range(max_gens):
    	    cost = run_particles(pop,tmpdir+'/'+prefix)
    	    for particle in pop:
    	        pso.update_velocity(particle, gbest, vel_space)
    	        pso.update_position(particle, search_space)
    	        particle["cost"] = cost[particle["id"]]
    	        pso.update_best_position(particle)
    	    gbest = pso.get_global_best(pop, gbest)
            f.write("%i %.6g " % (gen, gbest["cost"]))
    	    [f.write("%.6g "%p) for p in gbest["pos"]]
    	    f.write("\n")
    	    print " > gen %i, fitness=%.6g" % (gen, gbest["cost"])
            print names
            print gbest["pos"]
    except KeyboardInterrupt:
        print "\033[91mInterruption by user. Exiting...\033[0m"
    	pass
    f.close()
    return gbest

def search_es(max_gens, search_space, pop_size, num_children, prefix):
    tmpdir = tmp.gettempdir()
    f = open(prefix+'cost.txt','w')
    f.write("# gen fitness ")
    for n in names:
        f.write(n+" ")
    f.write("\n")
    population = es.init_population(search_space, pop_size)
    best = sorted(population, key=lambda k: k["fitness"])[0]
    try:
    	for gen in range(max_gens):
    	    children = [es.mutate(population[i], search_space, id=i) for i in range(num_children)]
    	    cost = run_particles(children,tmpdir+'/'+prefix)
    	    for c in children:
    	        c["fitness"] = cost[c["id"]]
    	    union = children+population
    	    union = sorted(union, key=lambda k: k["fitness"])
    	    if union[0]["fitness"] < best["fitness"]: best = union[0]
            population = union[:pop_size]
            f.write("%i %.6g " % (gen, best["fitness"]))
            [f.write("%.6g "%p) for p in best["pos"]]
    	    f.write("\n")
    	    print " > gen %i, fitness=%.4g" % (gen, best["fitness"])
            print names
            print best["pos"]
    except KeyboardInterrupt:
        print "\033[91mInterruption by user. Exiting...\033[0m"
    	pass
    f.close()
    return best

def search_ca(max_gens, search_space, pop_size, num_accepted, prefix):
    tmpdir = tmp.gettempdir()
    f = open(prefix+'cost.txt','w')
    f.write("# gen fitness ")
    for n in names:
        f.write(n+" ")
    f.write("\n")
    # initialize
    pop = [{"pos" : ca.random_vector(search_space) , "id" : i } for i in range(pop_size)]
    belief_space = ca.initialize_beliefspace(search_space)  
    # evaluate
    cost = run_particles(pop,tmpdir+'/'+prefix)
    for p in pop:
        p["fitness"] = cost[p["id"]]
    best = sorted(pop, key=lambda k: k["fitness"])[0]
    # update situational knowledge
    ca.update_beliefspace_situational(belief_space, best)
    try:
        for gen in range(max_gens):
            # create next generation
            children = [ ca.mutate_with_inf(pop[i], belief_space, search_space, id=i) for i in range(pop_size)]
            # evaluate
            cost = run_particles(children,tmpdir+'/'+prefix)
            for c in children:
                c["fitness"] = cost[c["id"]]
            best = sorted(children, key=lambda k: k["fitness"])[0]
            # update situational knowledge
            ca.update_beliefspace_situational(belief_space, best)
            # select next generation    
            pop = [ca.binary_tournament(children + pop) for i in range(pop_size)]
            # update normative knowledge
            pop = sorted(pop, key=lambda k: k["fitness"])
            acccepted = pop[0:num_accepted]
            ca.update_beliefspace_normative(belief_space, acccepted)
            # user feedback
            f.write("%i %.6g " % (gen, belief_space["situational"]["fitness"]))
            [f.write("%.6g "%p) for p in belief_space["situational"]["pos"]]
            f.write("\n")
            print " > gen %i, f=%.4g" % (gen, belief_space["situational"]["fitness"])
            print names
            print belief_space["situational"]["pos"]
    except KeyboardInterrupt:
        print "\033[91mInterruption by user. Exiting...\033[0m"
    	pass
    f.close()
    return belief_space["situational"]

if __name__ == "__main__":
    
    #Get information on parameters and ranges from user file
    params = pc.rdcsv('namelist.ranges')
    #format is:
    #param-name,low,high
    
    #set dimension
    k=len(params)
    
    # get names of parameters
    names=[]
    for i in range(0,k):
        names.append(params[i][0])

    #
    # Parsing command line arguments
    #

    inc=None
    if len(params[0]) == 4:
        inc = [params[i][-1] for i in range(len(params))]
        # create parameter space with respect to ranges provided in the list 'namelist.ranges'
        search_space = [[float(x) for x in p[1:-1]] for p in params]
    else:
        search_space = [[float(x) for x in p[1:]] for p in params]

    # create velocity space
    vel_space = generate_vel_space(search_space,frac=1)
    

    def plot(args,fnm):
        if args.plot: pc.plot_costs(fnm) 
        if args.save and args.plot: pc.plot_costs(fnm,savefig=True) 
        if args.save and not args.plot: pc.plot_costs(fnm,show=False,savefig=True) 


    def run_es(args):
        max_gens = args.max_gens
        pop_size = args.pop_size
        num_children = args.num_children
        print 'Run Evolution Strategies for %i generations of a population with size %i having %i children.' \
              % (max_gens, pop_size, num_children)
        best = search_es(max_gens, search_space, pop_size, num_children, 'es_')
        print  "done! Solution: f = %.4g, s =" % best["fitness"], best["pos"]
        plot(args,'es_cost.txt')


    def run_pso(args):
        max_gens = args.max_gens
        pop_size = args.pop_size
        print 'Run Particle Swarm Optimization for %i generations of population with size %i ' \
              % (max_gens, pop_size)
        best = search_pso(max_gens, search_space, vel_space, pop_size, 0.0, 0.0, 'pso_')
        print  "done! Solution: f = %.4g, s =" % best["cost"], best["pos"]
        plot(args,'pso_cost.txt')


    def run_ca(args):
        max_gens = args.max_gens
        pop_size = args.pop_size
        num_accepted = args.num_accepted
        print 'Run Cultural Algorithm for %i generations of population with size %i and from this the top %i as accepted knowledge' \
              % (max_gens, pop_size, num_accepted)
        best = search_ca(max_gens, search_space, pop_size, num_accepted, 'ca_')
        print  "done! Solution: f = %.4g, s =" % best["fitness"], best["pos"]
        plot(args,'ca_cost.txt')


    def run_rs(args):
        max_iter = args.max_iter
        print 'Run Random Search for %i iterations' % max_iter
        best = search_rs(search_space, max_iter, 'rs_', inc)
        print  "Done! Best solution: f = %.4g, v =" % best["cost"], best["pos"]
        plot(args,'rs_cost.txt')

    # Create main parser
    parser = argparse.ArgumentParser(description='Parameter optimization for SEMIC.')
    # Create parsers for each of the available optimization techniques
    sub_parser  = parser.add_subparsers(title='Available optimization strategies',help='choose one of those')

    es_parser = sub_parser.add_parser('es',help='Evolution Strategies')
    es_parser.add_argument('--max-gens','-g', help='Maximum number of generations (default=%(default)s).',type=int, default=100)
    es_parser.add_argument('--pop-size','-p', help='Size of population (default=%(default)s).',type=int, default=30)
    es_parser.add_argument('--num-children','-c', help='Number of children from each population (default=%(default)s).',type=int, default=20)
    es_parser.add_argument('--plot', help='Plot the optimal parameter matrix.', action='store_true', default=False)
    es_parser.add_argument('--save', help='Save the optimal parameter matrix.', action='store_true', default=False)
    es_parser.set_defaults(func=run_es)
    
    ca_parser = sub_parser.add_parser('ca',help='Cultural Algorithm')
    ca_parser.add_argument('--max-gens','-g', help='Maximum number of generations (default=%(default)s).',type=int, default=200)
    ca_parser.add_argument('--pop-size','-p', help='Size of population (default=%(default)s).',type=int, default=100)
    ca_parser.add_argument('--num-accepted','-n', help='Number of top solutions accepted as cultural knowledge (default=%(default)s).',type=int, default=20)
    ca_parser.add_argument('--plot', help='Plot the optimal parameter matrix.', action='store_true', default=False)
    ca_parser.add_argument('--save', help='Save the optimal parameter matrix.', action='store_true', default=False)
    ca_parser.set_defaults(func=run_ca)

    pso_parser = sub_parser.add_parser('pso',help='Particle Swarm Optimization')
    pso_parser.add_argument('--max-gens','-g', help='Maximum number of generations (default=%(default)s).',type=int, default=100)
    pso_parser.add_argument('--pop-size','-p', help='Size of population (default=%(default)s).',type=int, default=30)
    pso_parser.add_argument('--plot', help='Plot the optimal parameter matrix.', action='store_true', default=False)
    pso_parser.add_argument('--save', help='Save the optimal parameter matrix.', action='store_true', default=False)
    pso_parser.set_defaults(func=run_pso)

    rs_parser = sub_parser.add_parser('rs',help='Random Search')
    rs_parser.add_argument('--max-iter','-m', help='Number of iterations (default=%(default)s)',type=int, default=100)
    rs_parser.add_argument('--plot', help='Plot the optimal parameter matrix plot.', action='store_true', default=False)
    rs_parser.add_argument('--save', help='Save the optimal parameter matrix plot.', action='store_true', default=False)
    rs_parser.set_defaults(func=run_rs)

    # Parse and evaluate
    args = parser.parse_args()
    args.func(args)
