"""Tests for the distances module"""

import numpy as np
from ..distances import (general_distance_matrix, l2_distance_matrix,
                         l1_distance_matrix)

class TestGeneralDistanceFucntion:
    def test_l2_distance(self):
        pass

    def test_l1_distance(self):
        pass

    def test_cos_distance(self):
        pass

class TestL2DistanceFunction:
    def test_2d_canonical_base_vectors(self):
        X = np.array([
            [0, 1],
            [1, 0]
            ])
        expected_distances = np.array([
            [0, np.sqrt(2)],
            [np.sqrt(2), 0]
            ])
        distances = l2_distance_matrix(X, X)
        assert np.array_equal(expected_distances, distances)

    def test_3d_canonical_base_vectors(self):
        X = np.array([
            [0, 1],
            [1, 0]
            [0, 0]
            ])
        expected_distances = np.array([
            [0, np.sqrt(2)],
            [np.sqrt(2), 0]
            ])
        distances = l2_distance_matrix(X, X)
        assert np.array_equal(expected_distances, distances)

    def test_inline_vectors(self):
        X = np.array([
            [0, 0, 0, 0],
            [1, 2, 3, 4]
            ])
        expected_distances = np.array([
            [0, 1, 2, 3]
            [1, 0, 1, 2]
            [2, 1, 0, 1]
            [3, 2, 1, 0]
            ])
        distances = l2_distance_matrix(X, X)
        assert np.array_equal(expected_distances, distances)

    def test_pitagoric_triangle(self):
        X = np.array([
            [0, 4],
            [3, 0]
            ])
        expected_distances = np.array([
            [0, 5],
            [5, 0]
            ])
        distances = l2_distance_matrix(X, X)
        assert np.array_equal(expected_distances, distances)

class TestL1DistanceFunction:
    def test_canonical_base_vectors(self):
        pass

    def test_square_rooted_vectors(self):
        pass
