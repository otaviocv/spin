import numpy as np

def neighborhood(distances, weight_matrix, max_iter=100):
    permutation = np.identity(distances.shape[0])
    W = weight_matrix
    M = distances.dot(W)
    for i in range(max_iter):
        new_M = distances.dot(W)
        new_permutation = argmin_(M)
        if permutation.dot(M).trace() != new_permutation.dot(new_M).trace():
            W = permutation.T.dot(W)
    return permutation

def argmin_(M):
    return M

def side_to_side(distances, strictly_increasing_vector, max_iter=100):
    """ Side To Side SPIN algorithm

    Parameters
    ----------
    distances: np.array, shape [n, n]
        Distance symmetric square matrix.
    strictly_increasing_vector: np.array, shape [n]
        A vector with strictly increasing elements with the same dimension as
        the distance matrix.
    max_iter: int, default=100
        Maximum number of iterations.

    Returns
    -------
    permutation: np.array, shape [n, n]
        Permutation matrix with the same dimensions of the distance matrix.
    """
    X = strictly_increasing_vector
    permutation = np.identity(distances.shape[0])
    for i in range(max_iter):
        S = distances.dot(X)
        reverse_index_sort = (S).argsort()[::-1]
        new_permutation = np.identity(distances.shape[0])[reverse_index_sort]
        if np.all(new_permutation.dot(S) == permutation.dot(S)):
            break
        permutation = new_permutation
        X = permutation.dot(X)
    return permutation


def cost_function(distances, permutation, weight):
    return (permutation.dot(distances)
                       .dot(permutation.T)
                       .dot(weight)).trace()

def weight_matrix(strictly_increasing_vector):
    return strictily_increasing_vector.dot(strictly_increasing_vector.T)


def general_distance_matrix(X, dist_function):
    n = X.shape[1]
    dist_matrix = np.zeros((n, n), dtype=float)
    for i in range(0, n):
        for j in range(i, n):
            dist = dist_function(X[:, i], X[:, j])
            dist_matrix[i,j] = dist
            dist_matrix[j,i] = dist
        

def l2_distance_matrix(X, Y):
    dists = -2 * X.T.dot(Y) + \
            np.sum(X**2, axis=0) + \
            np.sum(Y**2, axis=0).reshape(1, -1).T
    dists[dists < 0] = 0
    return np.sqrt(dists)

def l1_distance_matrix(X, Y):
    pass