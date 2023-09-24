import numpy as np
import matplotlib.pyplot as plt
import ModSims

A = np.array([[0.7, -0.4], [-0.2, 0.5]])
b = np.array([0.3, 0.3])
x = np.array([21.0, -19.0])
solver = ModSims.Jacobi()
errors = solver.solve(A, b, x)

print(f"Number of iterations: {len(errors)}")

plt.plot(errors, 'r-')
plt.yscale("log")
plt.savefig("jac.png")