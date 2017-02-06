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
    k = len(minmax)
    pranges=[]
    for i in range(k):
        tmp=[]
        low=float(minmax[i][0])
        high=float(minmax[i][1])
        delta=(high-low)/float(max_iter)
        for j in range(max_iter):
            tmp.append(random.uniform(low+j*delta,low+(j+1)*delta))
        pranges.append(tmp)

    s=[]
    for i in range(k):
        s.append(range(max_iter))

    result=[]
    for i in range(max_iter):
        tmp=[]
        for j in range(k):
            a = random.sample(s[j],1)[0]
            tmp.append(a)
            s[j].remove(a)
        result.append(tmp)

    sample = []
    for l in range(len(result)):
        tmp=[]
        for j in range(len(result[l])):
            tmp.append(pranges[j][result[l][j]])
        sample.append(tmp)
    return sample

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
