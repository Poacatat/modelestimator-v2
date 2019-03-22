import argparse
import sys

from modelestimator.version import __version__

from modelestimator._bw_estimator.bw_estimator import bw_estimator
from modelestimator.io import handle_input_file, format_model_output
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
    ap.add_argument('-a', '--application', choices=['iqtree', 'matlab', 'mrbayes', 'octave', 'paml', 'phyml'], default='paml',
                    help="Choose output format to suit the application you want to use for inference. The 'iqtree', 'paml' and 'phyml' options are identical. The 'matlab' and 'octave' optins are for import into MatLab-compatible programs and are presenting the actual Q matrix rather than the R matrix used by PAML/PhyML, etc. Default: %(default)s")
    ap.add_argument('-f', '--format', choices=['fasta', 'clustal', 'nexus', 'phylip', 'stockholm'], default='fasta',
                    help="Specify sequence format of input files. Default: %(default)s")
    # ap.add_argument('-f', '--format', choices=['guess', 'fasta', 'clustal', 'nexus', 'phylip', 'stockholm'], default='guess',
    #                 help="Specify what sequence type to assume. Be specific if the file is not recognized automatically. When reading from stdin, the format is always guessed to be FASTA. Default: %(default)s")
    ap.add_argument('-t', '--threshold', type=float, metavar='T', default=0.001,
                    help="Stop when consecutive iterations do not change by more than T. Default: %(default)g")
    ap.add_argument('-b', '--bootstrap', type=int, metavar='N', default=None,
                    help="Perform bootstrapping on multialignment with N resamplings. Only one infile should be given in this mode. Returns bootstrap norm.")
    return ap


def controller(msa_format, output_format, bootstraps, threshold, filenames):
    '''
    Control the actions of the program.
    In: 
        msa_format - One of the allowed input formats, recognised by BioPython.
        bootstraps - None or an integer for the number of replacements in the bootstrap procedure.
        threshold  - Determines when to stop iterations, the minimum improvement before termination.
        
    '''
    multialignment_list = []
    
    for file in filenames:
        multialignment = handle_input_file(file, msa_format)
        multialignment_list.append(multialignment) 
        
    if threshold == None:
        threshold = 0.001
    
    if bootstraps:
        multialignment = multialignment_list[0]
        bootstrap_norm,_ = bootstraper(bootstraps, threshold, multialignment)
        output_string = "Bootstrap norm = " + str(bootstrap_norm)
    else:
        Q, EQ = bw_estimator(threshold, multialignment_list)
        output_string = format_model_output(Q, EQ, output_format)

    return output_string

def main():
    ap = setup_argument_parsing()
    args = ap.parse_args()
            
    try:
        output_string = controller(args.format, args.application, args.bootstrap, args.threshold, args.infiles)
    except Exception as e:
        print("Error: ", e)
        exit()

    print(output_string)
