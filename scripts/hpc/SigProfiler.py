#!/usr/bin/env python

from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerExtractor import sigpro as sig
import argparse
import os

parser = argparse.ArgumentParser(description = "This script runs the SigProfiler package.")
parser.add_argument("-i", "--input", required = True, help = "The input mut_mat.")
parser.add_argument("-o", "--output_dir", required = True, help = "The output directory.")
parser.add_argument("-c", "--context_type", required = True, help = "The context type.")
parser.add_argument("-p", "--penalty", required = True, help = "The nnls_remove_penalty")
parser.add_argument("--min_nr_signatures", required = False, default = 5, help = "The minimum amount of signatures that is extracted")
parser.add_argument("--max_nr_signatures", required = False, default = 20, help = "The maximum amount of signatures that is extracted")
args = parser.parse_args()

if __name__ ==  '__main__':

#    # Determine location input
#    cwd = os.getcwd()
#    input_mut_mat = cwd +"/simulated_data/mut_mat_" + args.name + ".txt"
#
#    # Determine output location
#    output_dir = cwd + "/sigprofiler_output_" + args.name
#    if not os.path.isdir(output_dir):
#        os.mkdir(output_dir)



    sig.sigProfilerExtractor("matrix",
    args.output_dir,
    args.input,
    context_type = args.context_type,
    nnls_remove_penalty = args.penalty,
    resample = True,
    nmf_replicates = 200,
    minimum_signatures = int(args.min_nr_signatures),
    maximum_signatures = int(args.max_nr_signatures),
    make_decomposition_plots = True,
    refit_denovo_signatures = True,
    get_all_signature_matrices = False,
    export_probabilities = True,
    cpu = 4)
