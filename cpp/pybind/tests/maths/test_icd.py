import numpy as np
from PyModSims.Matrix import IncompleteCholeskyDecomposition as ICD

class TestMatrixMatrix1:
    def setup_method(self):
        self.A = np.array([[0.7, -0.4], [-0.2, 0.5]])
        self.L = np.array([[0.83666,  0.0], [-0.239046, 0.665475]])
        self.F = np.array([[1.19523, 0], [-4.1833, 1.50269]])
        self.icd = ICD(self.A)
        
    def test_L(self):
        assert np.isclose(self.icd.getL(), self.L).all()

    def test_F(self):
        print(self.icd.getF())
        assert np.isclose(self.icd.getF(), self.F).all()
