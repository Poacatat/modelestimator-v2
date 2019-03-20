import argparse
import sys

from modelestimator.version import __version__

from modelestimator._bw_estimator.bw_estimator import bw_estimator
from modelestimator._handle_input.handle_input_file import handle_input_file
from modelestimator._bootstraper.bootstraper import bootstraper


USAGE_STRING ="""Usage: python -m modelestimator <format> <options> infiles

<format> should be either FASTA, STOCKHOLM or PHYLIP
Output is a rate matrix and residue distribution vector.
        
Options:
    -threshold or -t <f> Stop when consecutive iterations do not change by more
                     than <f>. Default is 0.001.
    -bootstrap or -b <r> Perform bootstrapping on multialignment with <r> resamplings.
                         Only one infile should be given in this mode. Returns
                         bootstrap norm.
"""

description='''
Infer a amino-acid replacement-rate matrix from one or more protein multisequence
alignments. The output is a rate matrix and an associated steady-state amino-acid
distribution vector.
'''
example_usage='''
Example usage:
    modelestimator fasta -t 0.001 file1.fa file2.fa file3.fa
    modelestimator fasta -b 200 file.fa
'''

def setup_argument_parsing():
    '''
    Create an argument parser.
    '''
    ap = argparse.ArgumentParser(description=description, epilog=example_usage,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('infiles', nargs='+', help="One or more infiles, containing protein multialignments in FASTA format or chosen according to the -f/--format option.")
    ap.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    ap.add_argument('-f', '--format', choices=['fasta', 'clustal', 'nexus', 'phylip', 'stockholm'], default='fasta',
                    help="Specify sequence format of input files. Default: %(default)s")
    # ap.add_argument('-f', '--format', choices=['guess', 'fasta', 'clustal', 'nexus', 'phylip', 'stockholm'], default='guess',
    #                 help="Specify what sequence type to assume. Be specific if the file is not recognized automatically. When reading from stdin, the format is always guessed to be FASTA. Default: %(default)s")
    ap.add_argument('-t', '--threshold', type=float, default=0.001,
                    help="Stop when consecutive iterations do not change by more than THRESHOLD. Default: %(default)f.")
    ap.add_argument('-b', '--bootstrap', type=int, default=None,
                    help="Perform bootstrapping on multialignment with BOOTSTRAP resamplings. Only one infile should be given in this mode. Returns bootstrap norm.")
    return ap


def controller(msa_format, bootstraps, threshold, filenames):
    '''
    Control the actions of the program.
    In: 
        msa_format - One of the allowed input formats, recognised by BioPython.
        bootstraps - None or an integer for the number of replacements in the bootstrap procedure.
        threshold  - Determines when to stop iterations, the minimum improvement before termination.
        
    '''
    MULTIALIGNMENT_LIST = []
    
    for FILE in filenames:
        MULTIALIGNMENT = handle_input_file(FILE, msa_format)
        MULTIALIGNMENT_LIST.append(MULTIALIGNMENT) 
        
    if threshold == None:
        threshold = 0.001
    
    if bootstraps:
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

def main():
    ap = setup_argument_parsing()
    args = ap.parse_args()
            
    try:
        OUTPUT_STRING = controller(args.format, args.bootstrap, args.threshold, args.infiles)
    except Exception as e:
        print("Error: ", e)
        exit()

    print(OUTPUT_STRING)
