import sys
import os
import io
from contextlib import redirect_stdout
import numpy as np
from modelestimator.main import main
from .integration_scripts import load_reference_output


def test_integration_1():
    #   Load REFERENCE_OUTPUT_STRING
    reference_output_string_path = "tests/test_integration_1/reference-output-string.txt"
    reference_output_string_path = os.path.join(sys.path[0], reference_output_string_path)

    with open(reference_output_string_path) as text_file:
            REFERENCE_OUTPUT_STRING = text_file.read()

    #   Calculate OUTPUT_STRING
    reference_input_string_path = "/tests/test_integration_1/test_integration_1_20seqs_1000long_50pam.fa"
    reference_input_string_path = sys.path[0] + reference_input_string_path
    sys.argv = ["-", "-f", "fasta", reference_input_string_path]

    f = io.StringIO()
    with redirect_stdout(f):
        main()

    OUTPUT_STRING = f.getvalue()
    REFERENCE_OUTPUT_LIST = load_reference_output(REFERENCE_OUTPUT_STRING)
    OUTPUT_LIST = load_reference_output(OUTPUT_STRING)
    THRESHOLD = 0.001
    
    assert len(REFERENCE_OUTPUT_LIST) == len(OUTPUT_LIST)
    for i in range(len(REFERENCE_OUTPUT_LIST)):
        assert(np.allclose(np.array(OUTPUT_LIST[i]), np.array(REFERENCE_OUTPUT_LIST[i]), atol=THRESHOLD))
