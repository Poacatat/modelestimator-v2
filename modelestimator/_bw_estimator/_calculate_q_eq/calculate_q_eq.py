import os
import sys
import numpy as np
from .comp_posterior_JC import comp_posterior_JC
from .matrix_weight import matrix_weight
from .estimate_q import estimate_q, find_zero_eigenvalue_eigenvector, find_eigens
from .simple_estimation import simple_estimation


def calculate_q_eq(count_matrix_list, threshold):
    #   VL = Inverse of right eigenvectorts
    #   VR = Right eigenvectors
    #   EQ = Right eigenvector corresponding to zero eigenvalue but normalized so it sums to 1
    vl, vr, eq = find_eigens(count_matrix_list)

    #   Get a first simple estimate of Q using a Jukes-Cantor model
    dist_samples = np.arange(1, 400, 5)
    posterior = comp_posterior_JC(count_matrix_list, dist_samples)   # posterior.shape = (10, 80). Rows are identical to Octave but in different order
    pw = matrix_weight(count_matrix_list, posterior, dist_samples)    
    w = posterior.sum(axis=0)
    q = estimate_q(pw, w, vl, vr, eq, dist_samples)

    #   Set loop variables
    difference = 1+threshold
    iterations = 0
    max_iterations = 10
    
    #   Calculate Q
    while (iterations < max_iterations and difference > threshold):
        iterations += 1
        q_new = simple_estimation(count_matrix_list, q, vl, vr, eq, dist_samples)
        difference = np.linalg.norm(q_new - q)
        q = q_new
    return q, eq
