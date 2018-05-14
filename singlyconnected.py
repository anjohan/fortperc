from walk import walk
import numpy as np
from scipy.ndimage import measurements
import multiprocessing as mp


def M_SC_one_sample(p, L):
    binary_matrix = np.random.uniform(size=(L, L)) < p
    label_matrix, num_clusters = measurements.label(binary_matrix)

    intersect_labels = np.intersect1d(label_matrix[0, :], label_matrix[-1, :])
    percolating_cluster = intersect_labels.max()

    if percolating_cluster == 0:
        label_matrix = label_matrix.T
        intersect_labels = np.intersect1d(label_matrix[0, :],
                                          label_matrix[-1, :])
        percolating_cluster = intersect_labels.max()

    if percolating_cluster == 0:
        return 0

    percolating_matrix = label_matrix == percolating_cluster

    left, right = walk(percolating_matrix)

    return np.count_nonzero(np.logical_and(left, right))


def M_SC(p, L, num_samples):
    pLs = [(p, L)] * num_samples
    with mp.Pool(4) as pool:
        result = pool.starmap(M_SC_one_sample, pLs)

    return np.sum(result) / num_samples
