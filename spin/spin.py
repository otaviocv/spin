import numpy as np


def general_distance_matrix(X, dist_function):
    n = X.shape[1]
    dist_matrix = np.zeros((n, n), dtype=float)
    for i in range(0, n):
        for j in range(i, n):
            dist = dist_function(X[:, i], X[:, j])
            dist_matrix[i,j] = dist
            dist_matrix[j,i] = dist
        

def l2_distance_matrix(X):
    dists = (np.dot(X.T, X))
    return dists

def l1_distance_matrix(X):
    pass
