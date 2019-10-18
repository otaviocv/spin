"""Main SPIN module with main classes and SPIN algorithms."""
import numpy as np
from .utils import check_distance_matrix

class SideToSideSPIN():
    """Side to side SPIN clustering method.

    Parameters
    ----------
    random_starts : int, optional (default=5)
        The number of different initial random permutations that will
        generated.
    max_iter : int, optional (default=100)
        The maximum number of iterations of each round of sorting.
    verbose : boolean, optional (default=False)
        Flag indicating to show logs and information during the SPIN process.
    
    Attributes
    ----------
    distances_ : array, shape (n, n)
        The original distances matrix provided.
    permutation_ : array, shape (n, n)
        Permutation matrix that can be applied to the original distances matrix
        to get to the ordered distances matrix.
    ordered_distances_ : array, shape (n, n)
        Distances matrix reordered by the permutation matrix. Before run this
        is the original distance matrix.

    References
    ----------
    D. Tsafrir, I. Tsafrir, L. Ein-Dor, O. Zuk, D.A. Notterman, E. Domany,
        Sortiug points into neighborhoods (SPIN): data analysis and
        visualization by ordering distance matrices, Bioinformatics, Volume 21,
        Issue 10, , Pages 2301–2308,
        https://doi.org/10.1093/bioinformatics/bti329

    """
    def __init__(self, random_starts=5, max_iter=100, verbose=False):
        self.random_starts = random_starts
        self.max_iter = max_iter
        self.verbose = verbose

    def run(self, X):
        if X.shape[0] != X.shape[1]:
            raise ValueError("The SPIN method only works with square matrices."
                             f"You provided a matrix of shape {X.shape}.")
        print("Setup")
        self.size_ = X.shape[0]
        self.distances_ = X
        self.permutation_ = np.identity(self.size_)
        self.ordered_distances_ = self.permutation_.dot(X) \
                                                   .dot(self.permutation_.T)
        assert np.array_equal(self.distances_, self.ordered_distances_)
        self.increasing_vector_ = np.array([i-(self.size_+1)/2
                                           for i in range(self.size_)]).reshape(-1, 1)
        self.weight_matrix_ = self.increasing_vector_ \
                                  .dot(self.increasing_vector_.T)
        print(self.weight_matrix_)
        self.energy_ = spin_energy(self.ordered_distances_,
                                   self.weight_matrix_)

        print(f"Initial energy: {self.energy_}")
        print("Actual spin")
        for i in range(self.random_starts):
            initial_permutation = random_permutation(self.size_)
            print(initial_permutation[:5, :5])
            permutation = side_to_side(self.distances_,
                                       self.increasing_vector_,
                                       initial_permutation,
                                       self.max_iter,
                                       self.verbose)
            if np.array_equal(permutation, initial_permutation):
                print("They are equal.")
            ordered_distances = permutation.dot(self.distances_) \
                                           .dot(permutation.T)
            energy = spin_energy(ordered_distances, self.weight_matrix_)
            print(f"{i}: {energy}")
            if energy < self.energy_:
                self.permutation_ = permutation
                self.ordered_distances_ = ordered_distances
                self.energy_ = energy

def random_permutation(size):
    """Random permutation matrix.
    
    Parameters
    ----------
    size : int
        The dimension of the random permutation matrix.

    Returns
    -------
    random_permutation : array, shape (size, size)
        An identity matrix with its rows random shuffled.

    """
    identity = np.identity(size)
    index = np.arange(0, size)
    np.random.shuffle(index)
    random_permutation = identity[index]
    return random_permutation


def side_to_side(distances, strictly_increasing_vector, initial_permutation,
                 max_iter=100, verbose=False):
    """Side To Side SPIN algorithm.

    Parameters
    ----------
    distances : np.array, shape [n, n]
        Distance symmetric square matrix.
    strictly_increasing_vector : np.array, shape [n]
        A vector with strictly increasing elements with the same dimension as
        the distance matrix.
    initial_permutation : array, shape [n ,n]
        The initial permutation matrix.
    max_iter : int, default=100
        Maximum number of iterations.
    verbose : bool
        Verbosity flag, if it is true print useful information about the
        process.

    Returns
    -------
    permutation : np.array, shape [n, n]
        Permutation matrix with the same dimensions of the distance matrix.

    """
    X = strictly_increasing_vector
    permutation = initial_permutation.copy()
    for i in range(max_iter):
        print(".", end="")
        S = distances.dot(X).flatten()
        reverse_index_sort = (S).argsort()[::-1]
        new_permutation = np.identity(distances.shape[0])[reverse_index_sort]
        if np.all(new_permutation.dot(S) == permutation.dot(S)):
            break
        permutation = new_permutation
        X = permutation.dot(X)
    return permutation

def spin_energy(ordered_distances, weight_matrix):
    """SPIN matrix energy.

    Metrices with better sorting have lower energies.

    Parameters
    ----------
    ordered_distances : array, shape (n, n)
        The ordered distances matrix. ordered_distances = PDP.T
    weight_matrix : array, shape (n, n)
        The weight matrix to weight the ordered distances matrix.

    Returns
    -------
    energy : float
        The energy of the associated ordered distance matrix and the
        weight matrix.
    """
    energy = np.trace(ordered_distances.dot(weight_matrix))
    return energy

class NeighborhoodSPIN():
    """Neighborhood SPIN clustering method.

    Parameters
    ----------
    initial_sigma : float, optional (default=2e10)
        Initial sigma value. This parameter controls the weight matrix
        dispersion.
    update_factor : float, optional (default=0.5)
        The number that will update the sigma value at each iteration. Sigma will
        be updated by sigma = sigma * update_factor.
    max_iter : int, optional (default=100)
        The maximum number of iterations of each round of sorting.
    verbose : boolean, optional (default=False)
        Flag indicating to show logs and information during the SPIN process.

    Attributes
    ----------
    distances_ : array, shape (n, n)
        The original distances matrix provided.
    permutation_ : array, shape (n, n)
        Permutation matrix that can be applied to the original distances matrix
        to get to the ordered distances matrix.
    ordered_distances_ : array, shape (n, n)
        Distances matrix reordered by the permutation matrix. Before run this
        is the original distance matrix.

    References
    ----------
    D. Tsafrir, I. Tsafrir, L. Ein-Dor, O. Zuk, D.A. Notterman, E. Domany,
        Sortiug points into neighborhoods (SPIN): data analysis and
        visualization by ordering distance matrices, Bioinformatics, Volume 21,
        Issue 10, , Pages 2301–2308,
        https://doi.org/10.1093/bioinformatics/bti329

    """
    def __init__(intial_sigma=2e10, update_factor=0.5, max_iter=100,
                 verbose=False):
        self.intial_sigma = intial_sigma
        self.update_factor = update_factor
        self.max_iter = max_iter
        self.verbose = verbose
        

class SPIN():
    """SPIN Clustering method.

    Parameters
    ----------
    method : str, {"sts", "neighborhood"}, optional, (default="sts")
        String determining which method to calculate permutations, either "sts"
        for side to side method or "neighborhood" for neighborhood method.
    verbose : boolean, optional (default=False)
        Flag indicating to show logs and information during the SPIN process.
    max_iter : int, optional (default=100)
        Maximum number of iterations to perform. It is needed for both methods.

    """

    def __init__(self, method="sts", verbose=False, max_iter=100):
        if method in {"sts", "neighborhood"}:
            self.method = method
        else:
            raise ValueError("The allowed methods are 'sts' and "
                             f" 'neighborhood', you provided {method}")
        self.verbose = verbose
        self.max_iter = max_iter

    def run(self, distances):
        """Calculate the permutation matrix.

        Calculate permutation matrix and apply to the original data with the
        specified method.

        Parameters
        ----------
        distances : array, shape (n_points, n_points)
            The distances symmetric square matrix.

        """
        check_distance_matrix(distances)

        self.distances_ = distances
        n_points = len(self.distances_)
        if self.method == "sts":
            self.permutation_ = side_to_side(self.distances_,
                                             self.increasing_vector_,
                                             self.max_iter,
                                             self.verbose)
            self.ordered_distances_ = self.permutation_ \
                                          .dot(self.distances_) \
                                          .dot(self.permutation_.T)
        elif self.method == "neighborhood":
            self.weight_matrix_ = initial_weight_matrix(n_points)
            self.permutation_ = neighborhood(self.distances_,
                                             self.weight_matrix_,
                                             self.max_iter,
                                             self.verbose)
            self.ordered_distances_ = self.permutation_ \
                                          .dot(self.distances_) \
                                          .dot(self.permutation_.T)
        return self

def initial_weight_matrix(size, sigma=1e2):
    """Initial weight matrix for neighborhood method.

    This initial matrix is initialized with exponential coefficients and then
    turned into a doubly stochastic matrix.
    
    Parameters
    ----------
    size : int
        The size of the initial weight matrix.
    sigma : float
        Coefficient to control dispersion of the weigth metrix coefficients.

    Returns
    -------
    weight_matrix : array, shape (size, size)
        The initial weight matrix. It is a square matrix.

    """
    rows_index_matrix, columns_index_matrix = np.indices((size, size))
    diff_index_matrix = rows_index_matrix - columns_index_matrix
    exp_arg_index_matrix = -(diff_index_matrix**2)/(size*sigma)
    non_normalized_weight_matrix = np.exp(exp_arg_index_matrix)
    print(non_normalized_weight_matrix)
    weight_matrix = sinkhorn_knopp_normalization_alogrithm(
            non_normalized_weight_matrix
            )
    return weight_matrix

def sinkhorn_knopp_normalization_alogrithm(matrix, tolerance=1e-5,
                                           max_iter=1000):
    """ The Sinkhorn Knopp algorithm to turn matrices into doubly stochastic.
    
    Parameters
    ----------
    matrix : array
        The matrix that will be normalized.
    tolerance : float
        The tolerance in the matrix approximation.
    max_iter : int
        If the tolerance is not reached this argument will set the maximun
        number of iterations.

    Returns
    -------
    norm_matrix : array
        The normalized version from the original matrix.

    References
    ----------
    Sinkhorn, Richard. A Relationship Between Arbitrary Positive Matrices and
        Doubly Stochastic Matrices. Ann. Math. Statist. 35 (1964), no. 2,
        876--879.  doi:10.1214/aoms/1177703591.
        https://projecteuclid.org/euclid.aoms/1177703591

    Sinkhorn, Richard, and Paul Knopp. "Concerning nonnegative matrices and
        doubly stochastic matrices." Pacific Journal of Mathematics 21.2
        (1967): 343-348.
        http://www.yaroslavvb.com/papers/sinkhorn-concerning.pdf
    """
    norm_matrix = matrix.copy()
    for i in range(max_iter):
        print(".", end="")
        col_sum = norm_matrix.sum(axis=0)
        norm_matrix = norm_matrix/col_sum
        row_sum = norm_matrix.sum(axis=1).reshape(-1, 1)
        norm_matrix = norm_matrix/row_sum

        if (np.all(np.abs(norm_matrix.sum(axis=1) - 1) < tolerance) and
            np.all(np.abs(norm_matrix.sum(axis=0) - 1) < tolerance)):
            break
    return norm_matrix


def neighborhood(distances, weight_matrix, max_iter=100, verbose=False):
    """Neighborhood SPIN algorithm.

    Parameters
    ----------
    distances : np.array, shape [n, n]
        Distance symmetric square matrix.
    weight_matrix : np.array, shape [n, n]
        A initial weight matrix to update permutaions matrix.
    max_iter : int, default=100
        Maximum number of iterations.
    verbose : bool
        Verbosity flag, if it is true print useful information about the
        process.

    Returns
    -------
    permutation : np.array, shape [n, n]
        Permutation matrix with the same dimensions of the distance matrix.

    """
    permutation = np.identity(distances.shape[0])
    W = weight_matrix
    M = distances.dot(W)
    for i in range(max_iter):

        if verbose:
            if i%100 == 0:
                print(f"\niter={i} ", end="")
            else:
                print(".", end="")

        new_M = distances.dot(W)
        new_permutation = np.identity(distances.shape[0])[np.argmin(new_M,
                                                                    axis=0)]
        if permutation.dot(M).trace() != new_permutation.dot(new_M).trace():
            W = permutation.T.dot(W)
        else:
            break
    return permutation


