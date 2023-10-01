import os
import numpy as np
import matplotlib.pyplot as plt
from ModSims import Jacobi

path = os.path.dirname(os.path.abspath(__file__))

def test_gauss_seidel_solver():
    A = np.array([[0.7, -0.4], [-0.2, 0.5]])
    b = np.array([0.3, 0.3])
    x = np.array([21.0, -19.0])

    solver = Jacobi()
    errors = solver.solve(A, b, x)

    assert len(errors) >= 45
    assert len(errors) <= 55
