from modelestimator._bw_estimator.bw_estimator import bw_estimator
from modelestimator._handle_input.handle_input_file import handle_input_file
from modelestimator._bootstraper.bootstraper import bootstraper

def controller(FORMAT, BOOTSTRAP, RESAMPLINGS, threshold, FILE_NAME_LIST):
    MULTIALIGNMENT_LIST = []
    
    for FILE in FILE_NAME_LIST:
        MULTIALIGNMENT = handle_input_file(FILE, FORMAT)
        MULTIALIGNMENT_LIST.append(MULTIALIGNMENT) 
        
    if threshold == None:
        threshold = 0.001
    
    if BOOTSTRAP:
        MULTIALIGNMENT = MULTIALIGNMENT_LIST[0]
        BOOTSTRAP_NORM,_ = bootstraper(RESAMPLINGS, threshold, MULTIALIGNMENT)
        OUTPUT_STRING = "Bootstrap norm = " + str(BOOTSTRAP_NORM)
    else:
        Q, EQ = bw_estimator(threshold, MULTIALIGNMENT_LIST)

        # Used to output the Q matrix:
        # OUTPUT_STRING = "# Q =\n" + str(Q) + "\n# EQ =\n" + str(EQ)   

        # ...but if outputting the R matrix (which is symmetrix and q_ij = r_ij*p_j)
        # we get compatibility with other tools, like PhyML.
        OUTPUT_STRING = ''
        for row in range(1,20):
            for col in range(0, row):
                r = Q[row, col] / EQ[col]
                OUTPUT_STRING += f'{r:10.5} '
            OUTPUT_STRING += '\n'

        EQ_STR = ''
        for elem in range(0,20):
            EQ_STR += f'{EQ[elem]:10.5} '

        OUTPUT_STRING += '\n' + EQ_STR

    return OUTPUT_STRING
