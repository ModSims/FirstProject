import os
import numpy as np
import matplotlib.pyplot as plt
from PyModSims import GaussSeidel

path = os.path.dirname(os.path.abspath(__file__))

def test_gauss_seidel_solver():
    A = np.array([[0.7, -0.4], [-0.2, 0.5]])
    x = np.array([21.0, -19.0])
    b = np.array([0.3, 0.3])

    dt = 0.001
    max_iterations = 100
    omega = 1.0
    solver = GaussSeidel.simulate(A, x, b, max_iterations, omega, dt)
    errors = solver.getResiduals()

    assert len(errors) >= 23
    assert len(errors) <= 27
