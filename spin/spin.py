import numpy as np

def side_to_side(distances, strictly_increasing_vector, max_iter=100):
    X = strictly_increasing_vector
    permutation = np.identity(distances.shape[0])
    for i in range(max_iter):
        S = distances.dot(X)
        reverse_index_sort = (S).argsort()[::-1]
        new_permutation = np.identity(distances.shape[0])[reverse_index_sort]
        if new_permutation.dot(S) == permutation.dot(S):
            break
        permutation = new_permutation
        X = permutation.dot(X)
    return permutation.dot(distances).dot(permutation.T)


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
    return np.sqrt(dists)

def l1_distance_matrix(X, Y):
    pass
