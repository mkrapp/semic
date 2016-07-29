'''

This is the Python version of the Random Search algorithm 
presented in

    "Clever Algorithms: Nature-Inspired Programming Recipes".


Find the book and the Ruby source codes on GitHub: 
  
   https://github.com/jbrownlee/CleverAlgorithm:s
'''
import random

def objective_function(vector):
   return sum([x ** 2.0 for x in vector])

def random_vector(minmax):
    return [minmax[i][0] + ((minmax[i][1] - minmax[i][0]) * random.random()) for i in range(len(minmax))]

def search(search_space, max_iter):
    best = None
    for i in range(max_iter):
        candidate = {}
        candidate["pos"] = random_vector(search_space)
        candidate["cost"] = objective_function(candidate["pos"])
        if best is None or candidate["cost"] < best["cost"]: best = candidate
        print " > iteration %i, best=%.4g" % (i, best["cost"])
    return best

if __name__ == "__main__":

    search_space = [[-5,5],[-5,5]]
    problem_size = 2
    max_iter = 100
    best = search(search_space, max_iter)
    print  "Done! Best solution: f = %.4g, v =" % best["cost"], best["pos"]
