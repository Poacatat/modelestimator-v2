from ._calculate_q_eq.match_closest_pair import choose_close_pairs
from ._calculate_q_eq.create_count_matrices import create_count_matrices
from ._calculate_q_eq.calculate_q_eq import calculate_q_eq
import numpy
import sys
#numpy.set_printoptions(threshold=sys.maxsize)
#import time


def bw_estimator(threshold, msa_list, compare_indels_flag = False):
    '''
    threshold  -  Decides when to stop iterations
    msa_list   -  A list of alignments 
    compare_indels_flag  -  decides if indels should be included when comparing likeness of sequences

    Returns: a rate matrix and a corresponding steady state distribution on amino acids,
             as estimated with the BW method.
    '''
    aggregated_count_matrix_list = []
   
    # MSA list then being just a list of fa file, the > lines being the breaks between sequences
    for msa in msa_list:
        #start = time.perf_counter()
        close_pairs = choose_close_pairs(msa, compare_indels_flag)
        #one1 = time.perf_counter()
        count_matrix_list = create_count_matrices(close_pairs)
        #two2 = time.perf_counter()
        aggregated_count_matrix_list.extend(count_matrix_list)
        #three3 = time.perf_counter()
    #so there we go, one stp closer to the truth. 
    Q, eq = calculate_q_eq(aggregated_count_matrix_list, threshold)
    #end = time.perf_counter()
    # 
    '''
    print(f"Choosing close pairs took {one1 - start:0.4f} seconds")
    print(f"Creating count matrices took {two2 - one1:0.4f} seconds")
    print(f"Aggregating count matrices took {three3 - two2:0.4f} seconds")
    print(f"Calculating Q and eq took {end - three3:0.4f} seconds")
    print(f"Total time: {end - start:0.4f} seconds")
    '''
    return Q, eq

