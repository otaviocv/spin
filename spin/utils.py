import numpy as np


def check_distance_matrix(distances):
    """ Perform all test to check if the distance matrix provided respects
    all constraints a distance matrix must have.

    Parameters
    ----------
    distances: array, shape (n_points, n_points)
        The distances symmetric square matrix.

    """

    # square matrix check
    if distances.shape[0] != distances.shape[1]:
        raise ValueError("Distance matrix provided is not square. "
                         f"The shape provided is {distances.shape}")

    # symmetry check
    if not is_symmetric(distances):
        raise ValueError("Distance matrix provided is not symmetric.")


def is_symmetric(matrix, rtol=1e-5, atol=1e-8):
    """ Check if a matrix is symmetric.

    Parameters
    ----------
    matrix: array, shape (n, n)
        matrix to be checked.
    rtol: float, optional (default=1e-5)
        Relative tolerance.
    atol: float, optional (default=1e-8)
        Absolute tolerance.
    """
    return np.allclose(a, a.T, rtol=rtol, atol=atol)
