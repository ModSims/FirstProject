import os
import numpy as np
import matplotlib.pyplot as plt
from PyModSims import Richardson

path = os.path.dirname(os.path.abspath(__file__))

def test_richardson_solver():
    A = np.array([[0.7, -0.4], [-0.2, 0.5]])
    x = np.array([21.0, -19.0])
    b = np.array([0.3, 0.3])

    dt = 0.001
    max_iterations = 100
    omega = 1.666 # lambda_max / lambda_min = 5/3
    solver = Richardson.simulate(A, x, b, max_iterations, omega, dt)
    errors = solver.getResiduals()

    assert len(errors) >= 50
    assert len(errors) <= 55
