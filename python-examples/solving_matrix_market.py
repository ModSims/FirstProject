from scipy.io import mmread
import random
import os
import numpy as np
import matplotlib.pyplot as plt
import ModSims

def generate_vector(matrix, randomize = False):
    if randomize:
        return [random.uniform(0,1) for i in range(matrix.shape[0])]
    else:
        return [random.uniform(2,9) for i in range(matrix.shape[0])]

A = mmread('bcsstm03.mtx').toarray()
x = generate_vector(A, randomize = True)
b = generate_vector(A, randomize = False)

print(A)

solver = ModSims.Jacobi()
errors = solver.solve(A, x, b)

path = os.path.dirname(os.path.abspath(__file__))
plt.plot(errors, 'r-')
plt.yscale("log")
if not os.path.exists(path + "/../out/plots"):
    os.makedirs(path + "/../out/plots")
plt.savefig(path + "/../out/plots/matrix_market.png")