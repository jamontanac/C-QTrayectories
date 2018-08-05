#!/usr/bin/env python
import numpy as np
import random
import math
import matplotlib.pyplot as plt

Datos=np.genfromtxt("a.txt")
X=Datos[:,0]
Y=Datos[:,1]
# plt.xlim([-5,5])
# plt.ylim([-5,5])
# for i in range(len(X)):
#     print X[i], Y[i]
plt.scatter(X,Y)
plt.show()
