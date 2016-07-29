'''

This is the Python version of the cultural algorithm 
presented in

    "Clever Algorithms: Nature-Inspired Programming Recipes".


Find the book and the Ruby source codes on GitHub: 
  
   https://github.com/jbrownlee/CleverAlgorithms

'''
import random


def objective_function(vector):
    return sum([x ** 2.0 for x in vector])


def rand_in_bounds(min, max):
    return min + ((max-min) * random.random()) 

def random_vector(minmax):
    return [rand_in_bounds(minmax[i][0], minmax[i][1]) for i in range(len(minmax))]

def mutate_with_inf(candidate, beliefs, minmax, id=None):
    v = []
    for i,c in enumerate(candidate["pos"]):
        v.append(rand_in_bounds(beliefs["normative"][i][0],beliefs["normative"][i][1]))
        if v[i] < minmax[i][0]: v[i] = minmax[i][0]
        if v[i] > minmax[i][1]: v[i] = minmax[i][1]
    c = {"pos" : v}
    if id is not None: c["id"] = id
    return c
  
def binary_tournament(pop):
    i, j = random.randint(0,len(pop)-1), random.randint(0,len(pop)-1)
    while j==i: j = random.randint(0,len(pop)-1)
    return pop[i] if (pop[i]["fitness"] < pop[j]["fitness"]) else pop[j]

def initialize_beliefspace(search_space):
    belief_space = {}
    belief_space["situational"] = None
    belief_space["normative"] = [search_space[i] for i in range(len(search_space))]
    return belief_space

def update_beliefspace_situational(belief_space, best):
    curr_best = belief_space["situational"]
    if curr_best is None  or best["fitness"] < curr_best["fitness"]:
        belief_space["situational"] = best

def update_beliefspace_normative(belief_space, acc):
    for i, bounds in enumerate(belief_space["normative"]):
        seq = [x["pos"][i] for x in acc]
        bounds[0] = min(seq)
        bounds[1] = max(seq)

def search(max_gens, search_space, pop_size, num_accepted):
    # initialize
    pop = [{"pos" : random_vector(search_space) } for i in range(pop_size)]
    belief_space = initialize_beliefspace(search_space)  
    # evaluate
    for p in pop:
        p["fitness"] = objective_function(p["pos"])
    best = sorted(pop, key=lambda k: k["fitness"])[0]
    # update situational knowledge
    update_beliefspace_situational(belief_space, best)
    for gen in range(max_gens):
        # create next generation
        children = [ mutate_with_inf(pop[i], belief_space, search_space) for i in range(pop_size)]
        # evaluate
        for c in children:
            c["fitness"] = objective_function(c["pos"])
        best = sorted(children, key=lambda k: k["fitness"])[0]
        # update situational knowledge
        update_beliefspace_situational(belief_space, best)
        # select next generation    
        pop = [binary_tournament(children + pop) for i in range(pop_size)]
        # update normative knowledge
        pop = sorted(pop, key=lambda k: k["fitness"])
        acccepted = pop[0:num_accepted]
        update_beliefspace_normative(belief_space, acccepted)
        # user feedback
        print " > gen %i, f=%.4g" % (gen, belief_space["situational"]["fitness"])
    return belief_space["situational"]

if __name__ == "__main__":
    # problem configuration
    problem_size = 2
    search_space = [[-5,5],[-5,5]]
    # algorithm configuration
    max_gens = 200
    pop_size = 100
    num_accepted = int(round(pop_size*0.20))
    # execute the algorithm
    best = search(max_gens, search_space, pop_size, num_accepted)
    print  "done! Solution: f = %.4g, s =" % best["fitness"], best["pos"]
