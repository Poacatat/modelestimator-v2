import numpy as np
import scipy
from scipy.special import logsumexp

### Private functions
def _precompute_pre_pt(q_matrix, eq_vector, dist_samples):
    """
    Pre-compute log(diag(eq) @ expm(Q * t)) for each distance sample.
    Uses an eigen-based exponential (cheap for 20x20) and falls back to expm
    if the spectrum looks unstable.
    """
    diag_eq = np.diag(eq_vector)

    def _as_real_if_close(matrix, tol=1e-6):
        real_matrix = np.real_if_close(matrix, tol=tol)
        if np.iscomplexobj(real_matrix):
            raise ValueError("Matrix exponential produced complex values")
        return real_matrix

    try:
        eigenvalues, vr = scipy.linalg.eig(q_matrix, left=False, right=True)
        vl = scipy.linalg.inv(vr)

        max_imag = max(
            np.max(np.abs(np.imag(eigenvalues))),
            np.max(np.abs(np.imag(vr))),
            np.max(np.abs(np.imag(vl))),
        )
        if max_imag > 1e-6:
            raise ValueError("Non-negligible imaginary part in eigendecomposition")

        eigenvalues = np.real(eigenvalues)
        vr = np.real(vr)
        vl = np.real(vl)

        pre_pt = []
        for dist_sample in dist_samples:
            exp_diag = np.exp(eigenvalues * dist_sample, dtype=np.float64)
            pt = diag_eq @ (vr @ (exp_diag[:, None] * vl))
            pt = _as_real_if_close(pt)
            pre_pt.append(np.log(np.clip(pt, np.float64(1e-12), None)))
        return np.asarray(pre_pt, dtype=np.float32)
    except Exception:
        pre_pt = []
        for dist_sample in dist_samples:
            pt = diag_eq @ scipy.linalg.expm(q_matrix * dist_sample)
            pt = _as_real_if_close(pt, tol=1e-3)
            pre_pt.append(np.log(np.clip(pt, np.float64(1e-12), None)))
        return np.asarray(pre_pt, dtype=np.float32)

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
    q64 = np.asarray(Q, dtype=np.float64)
    eq64 = np.asarray(EQ, dtype=np.float64)
    pre_pt = _precompute_pre_pt(q64, eq64, DIST_SAMPLES)

    MATRIX_LIST_LENGTH = len(COUNT_MATRIX_LIST)
    PD = np.empty((MATRIX_LIST_LENGTH, NUMBER_OF_DIST_SAMPLES), dtype=np.float32)
    PD = np.array([_my_posterior_pre(COUNT_MATRIX, pre_pt, DIST_SAMPLES) for COUNT_MATRIX in COUNT_MATRIX_LIST], dtype=np.float32)

    return PD
        
