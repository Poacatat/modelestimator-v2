import numpy as np
import scipy
from scipy.special import logsumexp

### Private functions

def _log_lik(COUNT_MATRIX, PT):
    return np.sum(COUNT_MATRIX * PT, dtype=np.float32)

# Compute the posterior probability of observing a set of replacements
#
# The integration code demands that the samples are uniformly distributed.
# Numerical integration using simple linear interpolation. 
#
# prePt is a list of pre-computed matrices Pt=expm(Q*t).
def _my_posterior_pre(COUNT_MATRIX, PRE_PT, DIST_SAMPLES):
    num_samples = len(DIST_SAMPLES)
    L = np.empty(num_samples, dtype=np.float32)

    for i, pre_pt_element in enumerate(PRE_PT):
        L[i] = _log_lik(COUNT_MATRIX, pre_pt_element)

    weights = np.ones(num_samples, dtype=np.float32)
    weights[0] = np.float32(0.5)
    weights[-1] = np.float32(0.5)

    log_weights = L + np.log(weights, dtype=np.float32)
    log_norm = logsumexp(log_weights).astype(np.float32) + np.log(
        DIST_SAMPLES[1] - DIST_SAMPLES[0], dtype=np.float32
    )

    posterior_vec = np.exp(log_weights - log_norm, dtype=np.float32)
    return posterior_vec


### Interface

# Given an estimate of Q, compute posterior probabilities for all
# distances for all seq pairs. 
#
# Similar to previous comp_posterior, but not re-computing matrix
# exponentials all the time. 
#
def comp_posterior(COUNT_MATRIX_LIST, Q, EQ, DIST_SAMPLES):   
    NUMBER_OF_DIST_SAMPLES = len(DIST_SAMPLES)
    q32 = np.asarray(Q, dtype=np.float32)
    eq32 = np.asarray(EQ, dtype=np.float32)
    PRE_PT = np.array([
        np.log(
            np.clip(
                np.diag(eq32) @ scipy.linalg.expm(q32 * np.float32(DIST_SAMPLE)).astype(np.float32),
                np.float32(1e-12),
                None,
            ),
            dtype=np.float32,
        )
        for DIST_SAMPLE in DIST_SAMPLES], dtype=np.float32)


    MATRIX_LIST_LENGTH = len(COUNT_MATRIX_LIST)
    PD = np.empty((MATRIX_LIST_LENGTH, NUMBER_OF_DIST_SAMPLES), dtype=np.float32)
    PD = np.array([_my_posterior_pre(COUNT_MATRIX, PRE_PT, DIST_SAMPLES) for COUNT_MATRIX in COUNT_MATRIX_LIST], dtype=np.float32)

    return PD
        
