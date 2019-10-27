# SPIN
*Sorting Points Into Neighborhoods Python Package*

[![Travis](https://api.travis-ci.org/otaviocv/spin.svg?branch=master)](https://travis-ci.org/otaviocv/spin)

## Introduction

**S**orting **P**oints **I**nto **N**eighborhoods is a clustering algorithm
different from the common options available. The method is taken from 

> D. Tsafrir, I. Tsafrir, L. Ein-Dor, O. Zuk, D.A. Notterman, E. Domany, Sorting
> points into neighborhoods (SPIN): data analysis and visualization by ordering
> distance matrices, Bioinformatics, Volume 21, Issue 10, , Pages 2301–2308,
> https://doi.org/10.1093/bioinformatics/bti329

This python relies on another implementation from
[jonataseduardo](https://github.com/jonataseduardo/SPIN).

## Installation

To intall the spin package just run

```
$ pip install spin-clustering
```

## Quick start

To run the any of the SPIN methods is super easy. Just try the following code.
*You will need some data, if you don't have any install scikit-learn and use one
of their examples, in this case we are using the iris dataset.*

```
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
```
