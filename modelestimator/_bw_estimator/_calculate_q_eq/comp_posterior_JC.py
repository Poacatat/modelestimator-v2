import numpy as np
from scipy.stats import binom
import sys
np.set_printoptions(threshold=sys.maxsize)

### Private functions
def _jc_posterior_ng(COUNT_MATRIX, DIST_SAMPLES):
    
    MATRIX_SUM = COUNT_MATRIX.sum()
    # here we for some reason get the total amount of transitions?
    MATRIX_DIAGONAL_SUM = COUNT_MATRIX.diagonal().sum()
    # maybe using .trace() is faster.
    # trace right?
    # in this case i guess it is the amount of times the amino acid does not evolve
    #DTYPECHANGE
    P = np.exp(- DIST_SAMPLES / 100, dtype=np.float32)
    # is the correct jukes cantor?
    # P_0(t) = 1/4 + 3/4 * exp(-4/3 * t)
   
    
    likelihood = binom.pmf(MATRIX_DIAGONAL_SUM, MATRIX_SUM, P)
    # fattar inte riktigt denna rad

#   This code is not commented in Octave
#    if (any(isnan(likelihood)))
#      # In case the binomial is tricky to compute, approx with normal distr!
#      likelihood = normpdf(k, tot .* p, tot .* p .* (1 .- p)); 
#    endif
    
    likelihood[0] /= 2
    likelihood[-1] /= 2
    # trapezoidal rule? numerisk integration?

    # DIST_SAMPLES[1] - DIST_SAMPLES[0] is the step size, 5

    # So we dived by delta x * sum, trap rule, so we basically normalise the likelihood
    POSTERIOR_VEC = likelihood / ( likelihood.sum() * (DIST_SAMPLES[1] - DIST_SAMPLES[0]) )

    return POSTERIOR_VEC    
    
### Interface
def comp_posterior_JC(COUNT_MATRIX_LIST, DIST_SAMPLES):
    NUMBER_OF_COUNT_MATRICES = len(COUNT_MATRIX_LIST)
    NUMBER_OF_DIST_SAMPLES = len(DIST_SAMPLES)
    #DTYPECHANGE
    PD = np.empty((NUMBER_OF_COUNT_MATRICES, NUMBER_OF_DIST_SAMPLES))
    #DTYPECHANGE
    PD = np.array([_jc_posterior_ng(COUNT_MATRIX, DIST_SAMPLES) for COUNT_MATRIX in COUNT_MATRIX_LIST], dtype=np.float32)
    return PD
        