'''

This is my Python version of a Latin-Hypercube Search algorithm 
not presented in

    "Clever Algorithms: Nature-Inspired Programming Recipes".

but based on the Random Search algorithm.

Find the book and the Ruby source codes on GitHub: 
  
   https://github.com/jbrownlee/CleverAlgorithm:s
'''
import random

def objective_function(vector):
   return sum([x ** 2.0 for x in vector])

def latin_hypercube(minmax,max_iter):
    segmentSize = float(len(minmax)) / float(max_iter)
    pos = []
    for j in range(max_iter):
        segmentMin = j * segmentSize
        point = segmentMin + (random.random() * segmentSize)
        pointValue = [(point * (minmax[i][1] - minmax[i][0])) + minmax[i][0] for i in range(len(minmax))]
        pos.append(pointValue)
    return pos
#    return [minmax[i][0] + ((minmax[i][1] - minmax[i][0]) * random.random()) for i in range(len(minmax))]

def search(search_space, max_iter):
    best = None
    lhs = latin_hypercube(search_space,max_iter)
    for i in range(max_iter):
        candidate = {}
        candidate["pos"] = lhs[i]
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
