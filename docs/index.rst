.. SPIN documentation master file, created by
   sphinx-quickstart on Sat Nov  9 09:57:11 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Welcome to SPIN's documentation!
================================

This is the Sorting Points Into Neighborhoods clustering method documentation.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   :ref:`introduction`
   api_reference

.. _introduction: 

Introduction
============

Sorting Points Into Neighborhoods, aka SPIN, is clustering
technique that only relies on the data and does not map any function on the
original space of the data points. :ref:`Introduction`

This implementation is based on the paper

    D. Tsafrir, I. Tsafrir, L. Ein-Dor, O. Zuk, D.A. Notterman, E. Domany,
    Sorting points into neighborhoods (SPIN): data analysis and visualization
    by ordering distance matrices, Bioinformatics, Volume 21, Issue 10, , Pages
    2301–2308, https://doi.org/10.1093/bioinformatics/bti329 

and in a previous implementation from jonataseduardo_.

The SPIN clustering method relies on defining an energy of a given sort of a
distance matrix. A distance matrix is a matrix which each entry i, j is the
distance between the sample i and the sample j. At each iteration SPIN tries
to reduce the matrix energy by changing the samples order, by the end natural
structures, if there are any, present on that shows up on the distance matrix
plot. Check the examples_ to understand the SPIN method better.

Installation
============

To install spin just use pip.

.. code-block:: shell

   $ pip install spin-clustering

If you can run the following import you have successfully installed the SPIN
library.

.. code-block:: python

   import spin


Quick Start
===========
To start playing around with SPIN there is a quick sample using the iris
dataset.

.. code-block:: python

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
    distances = l2_distance_matrix(X.T, X.T) # The l2 distance matrix computes
                                             # the l2 norm for two sets of
                                             # column vectors, so lets compute
                                             # all distances from all data
                                             # points.
    print(distances)

    print("Running SPIN...")
    spin = NeighborhoodSPIN()
    spin.run(distances)

    plt.matshow(spin.ordered_distances_)
    plt.show()


.. _jonataseduardo: https://github.com/jonataseduardo/SPIN
.. _examples: https://github.com/otaviocv/spin/tree/master/examples

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
