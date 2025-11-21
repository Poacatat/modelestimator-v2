import numpy as np
from scipy.linalg import eig

def find_zero_eigenvalue_eigenvector(matrix):
    '''
    Find the index of the eigenvector corresponding to Q's zero eigenvalue.
    This is recognized as the row (because we will be looking at the 'right'
    eigenvectors, not the usual left) with all positive or all negative elements.
    '''
    zero_eigenvector_candidates = []
    for index, eigen_vector in enumerate(matrix):
        if all(eigen_vector > 0) or all(eigen_vector < 0):
         zero_eigenvector_candidates.append((eigen_vector, index))
    if (len(zero_eigenvector_candidates) > 1):
        raise ValueError("More than one candidate for null-vector")
    if (len(zero_eigenvector_candidates) == 0):
        raise ValueError("No candidate for null-vector!")
        
    eigen_vector, index = zero_eigenvector_candidates[0]
    
    return np.asarray(eigen_vector, dtype=np.float32), index


def find_eigens(count_matrix_list):
    p_sum = sum((0.5 * (matrix + matrix.T)).astype(np.float32) for matrix in count_matrix_list) 
   
    # this is a bit wierd
    # anywho, it makes the matrix symmetric, which is maybe needed for the eigenvalue calculation

    # Make every row sum to 1
    row_sums = np.linalg.norm(p_sum, axis=1, ord=1, keepdims=1).astype(np.float32)
    
    p_sum = np.divide(p_sum, row_sums, out=np.zeros_like(p_sum, dtype=np.float32), where=row_sums != 0).astype(np.float32)   # Only divide where the row sum is non-zero

    
    try:
        eigen_values, vr = eig(p_sum.astype(np.float32), left=False, right=True)   #   Calculate eigenvalues and the right eigenvectors of p_sum
        
    except ValueError:
        raise ValueError("Unable to calculate eigenvalues")

    if not np.all(np.isreal(eigen_values)):
        raise ValueError("An eigenvalue is complex")
    eigen_values = np.real(eigen_values).astype(np.float32)
    

    # So so far we take the matrix make it symmetric, then normalise it, 
    # then calculate both the left and right eigenvector

    vr = np.real(vr).astype(np.float32)
    vl = np.linalg.inv(vr).astype(np.float32)

    eq,_ = find_zero_eigenvalue_eigenvector(vl)
    eq /= np.linalg.norm(eq, ord=1).astype(np.float32)
    eq = np.abs(eq, dtype=np.float32)
    


    

    
    #normalises the vector, with only positiv entrire
    # eq is therefore the zero, or closest to zero eigenvalues right eigenvector, normalised to sum to 1

    return vl, vr, eq


# Scale Q so that the average mutation rate is 0.01
def scale_q(Q, EQ):
        SCALE_FACTOR = np.dot(EQ, (-np.diag(Q))).astype(np.float32)

        if(np.isclose(SCALE_FACTOR, 0)):
            raise ZeroDivisionError('No Q diagonal cause a problem in estimate_q.py:scale_q')

        Q = Q.astype(np.float32, copy=False)
        Q /= SCALE_FACTOR
        return Q.astype(np.float32, copy=False)


#   Sometimes, when data is sparse, Q estimates come out with
#   off-diagonal entries being negative. Not good.
def fix_negatives(Q):
    Q = Q.astype(np.float32, copy=False)
    # Replace negative elements with smallest absolute value
    MINIMUM_ELEMENT = np.min(np.abs(Q))
    Q[Q<0] = MINIMUM_ELEMENT
    #   Recalculate diagonal to be negative rowsums
    np.fill_diagonal(Q, 0)
    ROW_SUMS = Q.sum(axis=1)
    np.fill_diagonal(Q, -ROW_SUMS)

    return Q


def recover_q(L, VR, VL):
    Q = np.float32(0.01) * (VR.astype(np.float32) @ np.diag(L.astype(np.float32)) @ VL.astype(np.float32))
    return Q.astype(np.float32)

# Alternative: Estimate eigenvalues of Q using weighted points
#
# The estimated eigenvalues are returned in L
def _weighted_estimate_eigenvals(PW, W, VL, VR, DIST_SAMPLES):

    # The X and Y matrises will contain 20 least-squares problems, one for each eigenvalue
    NUMBER_OF_DIST_SAMPLES = len(DIST_SAMPLES)
    X = np.empty(NUMBER_OF_DIST_SAMPLES, dtype="float32").reshape((80,1))
    Y = np.empty((20, NUMBER_OF_DIST_SAMPLES), dtype="float32")


    # Find the eigenvector corresponding to eigenvalue = 1
    _, NULL_VECTOR_INDEX = find_zero_eigenvalue_eigenvector(VL)
    # Gather some datapoints
    for i, DIST_SAMPLE in enumerate(DIST_SAMPLES):
        ELAMBDA = np.diag(VL @ PW[i] @ VR)
        X[i] = (DIST_SAMPLE / 100) * W[i]

        for li in range(20):
            if (li == NULL_VECTOR_INDEX):
                continue

            if (ELAMBDA[li] > 0):    # Skip complex value data points!
                Y[li, i] = np.real(np.log(ELAMBDA[li])) * W[i]
            else:
                X[i] = 0   # No disturbance if set to 0!
                Y[li, i] = 0
    L = np.zeros(20, dtype="float32")
    for i,_ in enumerate(L):
        if(i == NULL_VECTOR_INDEX):
            L[i] = 0
        else:
            tempY = Y[i,:].reshape(80,1)
            res = np.linalg.lstsq(X, tempY, rcond = None)[0][0][0]
            L[i] = np.float32(res)
    return L


### Interface
def estimate_q(PW, W, VL, VR, EQ, DIST_SAMPLES):
    L = _weighted_estimate_eigenvals(PW, W, VL, VR, DIST_SAMPLES)
    
    Q = recover_q(L, VR, VL)
    Q = fix_negatives(Q)
    Q = scale_q(Q, EQ)

    return Q
