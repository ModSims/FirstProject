import numpy as np
from PyModSims import Spectral

class TestSpectralMatrix1:
    def setup_method(self):
        self.A = np.array([[0.7, -0.4], [-0.2, 0.5]])
        
    def test_get_spectral_radius_via_eigen(self):
        expected_radius = 0.9
        spectral_radius = Spectral.getSpectralRadiusViaEigen(self.A)
        assert spectral_radius == expected_radius

    def test_get_gelfands_spectral_approximation(self):
        max_iterations = 100
        spectral_radius = Spectral.getGelfandsSpectralApproximation(self.A, max_iterations)
        assert 0.88 < spectral_radius < 0.92