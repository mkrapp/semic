'''

This is the Python version of the particle swarm optmization 
presented in

    "Clever Algorithms: Nature-Inspired Programming Recipes".


Find the book and the Ruby source codes on GitHub: 
  
   https://github.com/jbrownlee/CleverAlgorithms

Inertia/Constriction coefficient for w, c1, and c2 as in 
   Particle swarm optimization
   R Poli, J Kennedy, T Blackwell; Swarm intelligence; 2007; Springer
'''
import random

def objective_function(vector):
   return sum([x ** 2.0 for x in vector])

def random_vector(minmax):
    return [ minmax[i][0] + ((minmax[i][1] - minmax[i][0]) * random.random()) for i in range(len(minmax))]

def generate_vel_space(search_space,frac=1./5.):
    vel_space = []
    for i,s in enumerate(search_space):
    	vel_space.append([-frac*abs(s[1]-s[0])/2.,frac*abs(s[1]-s[0])/2.])
    return vel_space

def create_particle(search_space, vel_space, id=None):
    particle = {}
    if id is not None: particle["id"] = id
    particle["pos"]    = random_vector(search_space)
    particle["b_pos"]  = particle["pos"]
    particle["vel"]    = random_vector(vel_space)
    return particle

def get_global_best(population, current_best=None):
    best = sorted(population, key=lambda k: k["cost"])[0]
    if current_best is None or best["cost"] < current_best["cost"]:
        current_best = {}
        current_best["pos"]  = best["pos"]
        current_best["cost"] = best["cost"]
    return current_best

def update_velocity(part, gbest, vel_space, c1=1.49618, c2=1.49618, w=0.7298):
    for i,v in enumerate(part["vel"]):
        v1 = c1 * random.random() * (part["b_pos"][i] - part["pos"][i])
        v2 = c2 * random.random() * (gbest["pos"][i] - part["pos"][i])
        part["vel"][i] = w*v + v1 + v2
        #if part["vel"][i] > max_v: part["vel"][i] = max_v 
        #if part["vel"][i] < -max_v: part["vel"][i] = -max_v 
        'Instead of prescribing a global max velocity, we choose the bounds'
        'of the velocity space as upper/lower velocity limit'
        if part["vel"][i] > vel_space[i][1]: part["vel"][i] = vel_space[i][1] 
        if part["vel"][i] < vel_space[i][0]: part["vel"][i] = vel_space[i][0] 

def update_position(part, bounds):
    for i,p in enumerate(part["pos"]):
        part["pos"][i] = p + part["vel"][i]
        if part["pos"][i] > bounds[i][1]:
            part["pos"][i] = bounds[i][1] - abs(part["pos"][i] - bounds[i][1])
            part["vel"][i] *= -1.0
        if part["pos"][i] < bounds[i][0]:
            part["pos"][i] = bounds[i][0] + abs(part["pos"][i] - bounds[i][0])
            part["vel"][i] *= -1.0

def update_best_position(part):
    if part["cost"] < part["b_cost"]:
        part["b_cost"] = part["cost"]
        part["b_pos"] = list(part["pos"])

def search(max_gens, search_space, pop_size, c1, c2):
    # create velocity space
    vel_space = generate_vel_space(search_space,frac=1)
    pop = [create_particle(search_space, vel_space) for i in range(pop_size)]
    for p in pop:
        p["cost"] = objective_function(p["pos"])
        p["b_cost"] = p["cost"]
    gbest = get_global_best(pop)
    part_best = gbest
    for gen in range(max_gens):
        for particle in pop:
            update_velocity(particle, gbest, vel_space, c1, c2)
            update_position(particle, search_space)
            particle["cost"] = objective_function(particle["pos"])
            update_best_position(particle)
        gbest = get_global_best(pop, gbest)
        if gbest["cost"] < part_best["cost"]:
            part_best = dict(gbest)
        print " > gen %i, fitness=%.4g" % (gen, part_best["cost"])
    return part_best

if __name__ == "__main__":

    search_space = [[-5,5],[-5,5]]
    
    pop_size = 50
    max_gens = 100
    max_vel = 100.0
    c1, c2 = 2.0, 2.0
    best = search(max_gens, search_space, pop_size, c1,c2)
    print  "done! Solution: f = %.4g, s =" % best["cost"], best["pos"]
