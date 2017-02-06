import numpy as np
import argparse
from optimization_algorithms import particle_swarm_optimization as pso
from optimization_algorithms import evolution_strategies as es
from optimization_algorithms import cultural_algorithm as ca
from optimization_algorithms import random_search as rs
from optimization_algorithms import latin_hypercube_search as lhs
from optimization_algorithms import harmony_search as hs
import plot_costs as pc
from tools import run_particles, get_params, init_search


def search_rs(search_space, max_iter, prefix):
    f, tmpdir = init_search(prefix,names)
    pop = []
    for i in range(max_iter):
        candidate = {}
        candidate["id"] = i
        candidate["pos"] = rs.random_vector(search_space)
        pop.append(candidate)
    best = None
    cost = run_particles(pop,tmpdir+'/'+prefix)
    for i,candidate in enumerate(pop):
        candidate["cost"] = cost[candidate["id"]]
        if best is None or candidate["cost"] < best["cost"]: best = candidate
        f.write("%i %.6g " % (i, best["cost"]))
        [f.write("%.6g "%p) for p in best["pos"]]
        f.write("\n")
        print " > iteration %i, best=%.4g" % (i, best["cost"])
        print [names[i] + ": %.6g " % p for i,p in enumerate(best["pos"])]
    f.close()
    return best

def search_lhs(search_space, max_iter, prefix):
    f, tmpdir = init_search(prefix,names)
    pop = []
    lh = lhs.latin_hypercube(search_space,max_iter)
    for i in range(max_iter):
        candidate = {}
        candidate["id"] = i
        candidate["pos"] = lh[i]
        pop.append(candidate)
    best = None
    cost = run_particles(pop,tmpdir+'/'+prefix)
    for i,candidate in enumerate(pop):
        candidate["cost"] = cost[candidate["id"]]
        if best is None or candidate["cost"] < best["cost"]: best = candidate
        f.write("%i %.4g " % (i, best["cost"]))
        [f.write("%.4g "%p) for p in best["pos"]]
        f.write("\n")
        print " > iteration %i, best=%.4g" % (i, best["cost"])
        print [names[i] + ": %.4g " % p for i,p in enumerate(best["pos"])]
    f.close()
    return best

def search_pso(max_gens, search_space, pop_size, c1, c2, prefix):
    f, tmpdir = init_search(prefix,names)
    vel_space = pso.generate_vel_space(search_space,frac=1)
    pop = [pso.create_particle(search_space, vel_space, id=i) for i in range(pop_size)]
    cost = run_particles(pop,tmpdir+'/'+prefix)
    for p in pop:
        p["cost"] = cost[p["id"]]
        p["b_cost"] = p["cost"]
    gbest = pso.get_global_best(pop)
    part_best = {}
    part_best["pos"] = gbest["pos"]
    part_best["cost"] = gbest["cost"]
    try:
    	for gen in range(max_gens):
    	    for particle in pop:
    	        pso.update_velocity(particle, gbest, vel_space)
    	        pso.update_position(particle, search_space)
    	    cost = run_particles(pop,tmpdir+'/'+prefix)
    	    for particle in pop:
    	        particle["cost"] = cost[particle["id"]]
    	        pso.update_best_position(particle)
    	    gbest = pso.get_global_best(pop, gbest)
            if gbest["cost"] < part_best["cost"]:
                part_best["pos"] = list(gbest["pos"])
                part_best["cost"] = gbest["cost"]
            f.write("%i %.6g " % (gen, part_best["cost"]))
    	    [f.write("%.6g "%p) for p in part_best["pos"]]
    	    f.write("\n")
    	    print " > gen %i, fitness=%.6g" % (gen, part_best["cost"])
            print [names[i] + ": %.6g " % p for i,p in enumerate(part_best["pos"])]
    except KeyboardInterrupt:
        print "\033[91mInterruption by user. Exiting...\033[0m"
    	pass
    f.close()
    return part_best

def search_es(max_gens, search_space, pop_size, num_children, prefix):
    f, tmpdir = init_search(prefix,names)
    population = es.init_population(search_space, pop_size, id=True)
    cost = run_particles(population,tmpdir+'/'+prefix)
    for p in population:
        p["fitness"] = cost[p["id"]]
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
            print [names[i] + ": %.6g " % p for i,p in enumerate(best["pos"])]
    except KeyboardInterrupt:
        print "\033[91mInterruption by user. Exiting...\033[0m"
    	pass
    f.close()
    return best

def search_ca(max_gens, search_space, pop_size, num_accepted, prefix):
    f, tmpdir = init_search(prefix,names)
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
            print [names[i] + ": %.6g " % p for i,p in enumerate(belief_space["situational"]["pos"])]
    except KeyboardInterrupt:
        print "\033[91mInterruption by user. Exiting...\033[0m"
    	pass
    f.close()
    return belief_space["situational"]

def search_hs(search_space, max_iter, mem_size, consid_rate, adjust_rate, nrange, prefix):
    f, tmpdir = init_search(prefix,names)
    memory = [hs.create_random_harmony(search_space,id=i) for i in range(mem_size*3)]
    cost = run_particles(memory,tmpdir+'/'+prefix)
    for m in memory:
        m["fitness"] = cost[m["id"]]
    memory = sorted(memory, key=lambda k: k["fitness"])[:mem_size]
    best = memory[0]
    try:
        for iter in range(max_iter):
            harm = hs.create_harmony(search_space, memory, consid_rate, adjust_rate, nrange)
            harm["id"] = 0
            cost = run_particles([harm],tmpdir+'/'+prefix)
            harm["fitness"] = cost[harm["id"]]
            if harm["fitness"] < best["fitness"]: best = harm 
            memory.append(harm)
            memory = sorted(memory, key=lambda k: k["fitness"])
            del memory[-1]
            f.write("%i %.6g " % (iter, best["fitness"]))
            [f.write("%.6g "%p) for p in best["pos"]]
    	    f.write("\n")
            print " > iteration %i, fitness=%.4g" % (iter, best["fitness"])
            print [names[i] + ": %.6g " % p for i,p in enumerate(best["pos"])]
    except KeyboardInterrupt:
        print "\033[91mInterruption by user. Exiting...\033[0m"
    	pass
    f.close()
    return best

if __name__ == "__main__":
    

    params, names = get_params()

    #
    # Parsing command line arguments
    #

    search_space = [[float(x) for x in p[1:]] for p in params]

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
        best = search_pso(max_gens, search_space, pop_size, 0.0, 0.0, 'pso_')
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
        best = search_rs(search_space, max_iter, 'rs_')
        print  "Done! Best solution: f = %.4g, v =" % best["cost"], best["pos"]
        plot(args,'rs_cost.txt')


    def run_lhs(args):
        max_iter = args.max_iter
        print 'Run Latin Hypercube Search for %i iterations' % max_iter
        best = search_lhs(search_space, max_iter, 'lhs_')
        print  "Done! Best solution: f = %.4g, v =" % best["cost"], best["pos"]
        plot(args,'rs_cost.txt')


    def run_hs(args):
        mem_size = args.mem_size
        consid_rate = args.consid_rate
        adjust_rate = args.adjust_rate
        nrange = args.range
        max_iter = args.max_iter
        best = search_hs(search_space, max_iter, mem_size, consid_rate, adjust_rate, nrange, 'hs_')
        print  "done! Solution: f = %.4g, s =" % best["fitness"], best["pos"]
        plot(args,'hs_cost.txt')
        

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

    lhs_parser = sub_parser.add_parser('lhs',help='Latin Hypercube Search')
    lhs_parser.add_argument('--max-iter','-m', help='Number of iterations (default=%(default)s)',type=int, default=100)
    lhs_parser.add_argument('--plot', help='Plot the optimal parameter matrix plot.', action='store_true', default=False)
    lhs_parser.add_argument('--save', help='Save the optimal parameter matrix plot.', action='store_true', default=False)
    lhs_parser.set_defaults(func=run_lhs)

    hs_parser = sub_parser.add_parser('hs',help='Harmony Search')
    hs_parser.add_argument('--max-iter','-m', help='Number of iterations (default=%(default)s)',type=int, default=500)
    hs_parser.add_argument('--consid_rate','-c', help='Harmony Memory consideration rate, typically between 0.7 and 0.95 (default=%(default)s).',type=float, default=0.95)
    hs_parser.add_argument('--adjust_rate','-a', help='Pitch adjustment rate, typically between 0.1 and 0.5 (default=%(default)s).',type=float, default=0.5)
    hs_parser.add_argument('--range','-r', help='Random adjustemt to Pitch adjustemnt (default=%(default)s).',type=float, default=0.05)
    hs_parser.add_argument('--mem-size','-s', help='Size of memory (default=%(default)s).',type=int, default=20)
    hs_parser.add_argument('--plot', help='Plot the optimal parameter matrix plot.', action='store_true', default=False)
    hs_parser.add_argument('--save', help='Save the optimal parameter matrix plot.', action='store_true', default=False)
    hs_parser.set_defaults(func=run_hs)

    # Parse and evaluate
    args = parser.parse_args()
    args.func(args)
