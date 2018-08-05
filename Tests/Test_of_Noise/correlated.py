#!/usr/bin/env python
import numpy as np
import random
import math
import matplotlib.pyplot as plt
def get_gaussian_random():
    m = 0
    while m == 0:
        m = round(np.random.random() * 100)
        # print m
    numbers = np.random.random(int(m))
    summation = float(np.sum(numbers))
    gaussian = (summation - m/2) / math.sqrt(m/12.0)
    
    return gaussian
def generate_known_gaussian(dimensions):
    count = 1000
    ret = []
    for i in range(count):
        current_vector = []
        for j in range(dimensions):
            g= random.gauss(0.0,1.0)
            current_vector.append(g)
        ret.append(tuple(current_vector))
    return ret
def learn_mean_cov(pts):
    learned_mean = np.matrix([[0.0], [0.0]])
    learned_cov  = np.zeros( (2, 2) )
    count = len(pts)
    for pt in pts:
        learned_mean += pt
        learned_cov  += pt*pt.transpose()

    learned_mean /= count
    learned_cov /= count
    learned_cov -= learned_mean * learned_mean.transpose()

    print(learned_mean)
    print(learned_cov)
def main():
    known = generate_known_gaussian(2)
    target_mean = np.matrix([ [0.0], [0.0]])
    target_cov  = np.matrix([[1.0, 0.7],[0.7, 1.0]])
    [eigenvalues, eigenvectors] = np.linalg.eig(target_cov)
    l = np.matrix(np.diag(np.sqrt(eigenvalues)))
    Q = np.matrix(eigenvectors) * l
    print Q
    x1_tweaked = []
    x2_tweaked = []
    tweaked_all = []
    for i, j in known:
        original = np.matrix( [[i], [j]]).copy()
        tweaked = (Q*original) + target_mean
        # print tweaked
        x1_tweaked.append(float(tweaked[0]))
        x2_tweaked.append(float(tweaked[1]))
        tweaked_all.append( tweaked )
    # learn_mean_cov(tweaked_all)
    # plt.scatter(x1_tweaked, x2_tweaked)
    # plt.show()
if __name__== "__main__":
    main()
