'''

This is the Python version of the evolution strategies algorithm 
presented in

    "Clever Algorithms: Nature-Inspired Programming Recipes".


Find the book and the Ruby source codes on GitHub: 
  
   https://github.com/jbrownlee/CleverAlgorithms
'''
import random
import numpy as np

def objective_function(vector):
   return sum([x ** 2.0 for x in vector])

def random_vector(minmax):
    return [ minmax[i][0] + ((minmax[i][1] - minmax[i][0]) * random.random()) for i in range(len(minmax))]

def random_gaussian(mean=0.0, stdev=1.0):
    return random.gauss(mean,stdev)

def mutate_problem(vector, stdevs, search_space):
    child = []
    for i, v in enumerate(vector):
        child.append(v + stdevs[i] * random_gaussian())
	if child[i] < search_space[i][0]: child[i] = search_space[i][0]
	if child[i] > search_space[i][1]: child[i] = search_space[i][1]
    return child

def mutate_strategy(stdevs):
    tau = np.sqrt(2.0*len(stdevs))**-1.0
    tau_p = np.sqrt(2.0*np.sqrt(len(stdevs)))**-1.0
    child = [stdevs[i] * np.exp(tau_p*random_gaussian() + tau*random_gaussian()) 
             for i in range(len(stdevs))]
    return child

def mutate(par, minmax, id=None):
    child = {}
    if id is not None: child["id"] = id
    child["pos"] = mutate_problem(par["pos"], par["strategy"], minmax)
    child["strategy"] = mutate_strategy(par["strategy"])
    return child

def init_population(minmax, pop_size, id=False):
    strategy = [[0, (minmax[i][1]-minmax[i][0]) * 0.05] for i in range(len(minmax))]
    pop = [{} for i in range(pop_size)]
    for i, p in enumerate(pop):
        p["pos"] = random_vector(minmax)
        p["strategy"] = random_vector(strategy)
        if id: p["id"] = i
    return pop

def search(max_gens, search_space, pop_size, num_children):
    population = init_population(search_space, pop_size)
    for p in population:
        p["fitness"] = objective_function(p["pos"])
    population = sorted(population, key=lambda k: k["fitness"])
    best = population[0]
    for gen in range(max_gens):
        children = [mutate(population[i], search_space) for i in range(num_children)]
        for c in children:
	        c["fitness"] = objective_function(c["pos"])
        union = children+population
        union = sorted(union, key=lambda k: k["fitness"])
        if union[0]["fitness"] < best["fitness"]: best = union[0]
        population = union[:pop_size]
        print " > gen %i, fitness=%.4g" % (gen, best["fitness"])
    return best

if __name__ == "__main__":

    search_space = [[-5,5],[-5,5]]

    max_gens = 100
    pop_size = 30
    num_children = 20

    best = search(max_gens, search_space, pop_size, num_children)
    print  "done! Solution: f = %.4g, s =" % best["fitness"], best["pos"]
