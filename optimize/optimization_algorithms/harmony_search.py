'''

This is the Python version of the Harmony Search algorithm 
presented in

    "Clever Algorithms: Nature-Inspired Programming Recipes".


Find the book and the Ruby source codes on GitHub: 
  
   https://github.com/jbrownlee/CleverAlgorithm:s
'''
import random

def objective_function(vector):
   return sum([x ** 2.0 for x in vector])

def rand_in_bounds(min, max):
    return min + ((max-min) * random.random()) 

def random_vector(minmax):
    return [rand_in_bounds(minmax[i][0], minmax[i][1]) for i in range(len(minmax))]

def create_random_harmony(search_space, id=None):
    harmony = {}
    harmony["pos"] = random_vector(search_space)
    if id is not None: harmony["id"] = id
    return harmony

def create_harmony(search_space, memory, consid_rate, adjust_rate, nrange):
    vector = [0] * len(search_space)
    for i in range(len(search_space)):
        if random.random() < consid_rate:
            value = memory[random.randrange(len(memory))]["pos"][i]
            if random.random() < adjust_rate: value = value + nrange*rand_in_bounds(-1.0, 1.0) 
            if value < search_space[i][0]: value = search_space[i][0] 
            if value > search_space[i][1]: value = search_space[i][1] 
            vector[i] = value
    else:
        vector[i] = rand_in_bounds(search_space[i][0], search_space[i][1])
    return {"pos" : vector}

def search(bounds, max_iter, mem_size, consid_rate, adjust_rate, nrange):
    memory = [create_random_harmony(search_space) for i in range(mem_size*3)]
    for m in memory:
        m["fitness"] = objective_function(m["pos"])
    memory = sorted(memory, key=lambda k: k["fitness"])[:mem_size]
    best = memory[0]
    for iter in range(max_iter):
        harm = create_harmony(bounds, memory, consid_rate, adjust_rate, nrange)
        harm["fitness"] = objective_function(harm["pos"])
        if harm["fitness"] < best["fitness"]: best = harm 
        memory.append(harm)
        memory = sorted(memory, key=lambda k: k["fitness"])
        del memory[-1]
        print " > iteration %i, fitness=%.4g" % (iter, best["fitness"])
    return best

if __name__ == "__main__":
    # problem configuration
    problem_size = 3
    search_space = [[-5,5],[-5,5]]
    # algorithm configuration
    mem_size = 20
    consid_rate = 0.95
    adjust_rate = 0.7
    nrange = 0.05
    max_iter = 500
    # execute the algorithm
    best = search(search_space, max_iter, mem_size, consid_rate, adjust_rate, nrange)
    print  "done! Solution: f = %.4g, s =" % best["fitness"], best["pos"]
