import pytest as pt
import numpy as np

from ..utils import check_distance_matrix


class TestCheckDistanceMatrix:

    not_square_matrix_test_data = [
            pt.param(np.zeros((1, 2)), marks=pt.mark.xfail),
            pt.param(np.ones((2, 3)), marks=pt.mark.xfail),
            pt.param(np.zeros((5, 4)), marks=pt.mark.xfail),
            pt.param(np.ones((3, 9)), marks=pt.mark.xfail),
            pt.param(np.zeros((10, 2)), marks=pt.mark.xfail),
            pt.param(np.ones((5, 6)), marks=pt.mark.xfail),
            ]

    @pt.mark.parametrize("matrix", not_square_matrix_test_data)
    def test_raise_value_error_for_not_square_matrices(self, matrix):
        check_distance_matrix(matrix)
