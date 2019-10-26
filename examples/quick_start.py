from spin import NeighborhoodSPIN
from spin.distances import l2_distance_matrix
from sklearn.datasets import load_iris
import matplotlib.pyplot as plt

# Loading data
print("Loading data...")
data = load_iris()
X = data.data
y = data.target

print("Calculating distances...")
distances = l2_distance_matrix(X.T, X.T) # The l2 distance matrix computes the
                                         # l2 norm for two sets of column
                                         # vectors, so lets compute all
                                         # distances from all data points.
print(distances)

print("Running SPIN...")
spin = NeighborhoodSPIN()
spin.run(distances)

plt.matshow(spin.ordered_distances_)
plt.show()
