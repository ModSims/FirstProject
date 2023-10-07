import numpy as np
from PyModSims import Matrix

class TestMatrixMatrix1:
    def setup_method(self):
        self.A = np.array([[0.7, -0.4], [-0.2, 0.5]])
        
    def test_get_condition_number(self):
        expected_cond = 3.666
        cond = Matrix.getConditionNumber(self.A)
        assert cond <= expected_cond

    def test_is_symmetric(self):
        assert not Matrix.isSymmetric(self.A)

    def test_is_diagonally_dominant(self):
        assert Matrix.isDiagonallyDominant(self.A)

    def test_is_SPD(self):
        assert not Matrix.isSPD(self.A)

class TestMatrixMatrix2:
    def setup_method(self):
        self.A = np.array([[0.7, -0.2], [-0.2, 0.5]])
        
    def test_get_condition_number(self):
        expected_cond = 3.666
        cond = Matrix.getConditionNumber(self.A)
        assert cond <= expected_cond

    def test_is_symmetric(self):
        assert Matrix.isSymmetric(self.A)

    def test_is_diagonally_dominant(self):
        assert Matrix.isDiagonallyDominant(self.A)

    def test_is_SPD(self):
        assert Matrix.isSPD(self.A)