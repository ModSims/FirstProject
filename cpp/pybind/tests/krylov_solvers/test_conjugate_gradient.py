import os
import numpy as np
from PyModSims import ConjugateGradient

path = os.path.dirname(os.path.abspath(__file__))

def test_conjugate_gradient_solver():
    A = np.array([[0.7, -0.4], [-0.4, 0.7]])
    x = np.array([21.0, -19.0])
    b = np.array([0.3, 0.3])

    dt = 0.001
    solver = ConjugateGradient.simulate(A, x, b, dt)
    errors = solver.getAlphas()

    assert len(errors) == 2
