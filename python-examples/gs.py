import os
import numpy as np
import matplotlib.pyplot as plt
import PyModSims

path = os.path.dirname(os.path.abspath(__file__))

A = np.array([[0.7, -0.4], [-0.2, 0.5]])
b = np.array([0.3, 0.3])
x = np.array([21.0, -19.0])
solver = PyModSims.GaussSeidel()
errors = solver.solve(A, x, b)

plt.plot(errors, 'r-')
plt.yscale("log")
if not os.path.exists(path + "/../out/plots"):
    os.makedirs(path + "/../out/plots")
plt.savefig(path + "/../out/plots/gs.png")
