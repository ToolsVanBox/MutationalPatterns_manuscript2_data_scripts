#!/usr/bin/env python

from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerExtractor import sigpro as sig
import argparse
import os

parser = argparse.ArgumentParser(description = "This script runs the SigProfiler package.")
parser.add_argument("-n", "--name", required = True, help = "The experiment name.")
args = parser.parse_args()

# Create matrices
#matrices = matGen.SigProfilerMatrixGeneratorFunc(
#"test1",
#"GRCh37",
#"/Users/freekmanders/surfdrive/Shared/Boxtel_General/Scripts/Git_submission/Freek_MutationalPatterns/test_samples_signatureanalyzer",
#exome=False,
#bed_file=None,
#chrom_based=False,
#plot=True,
#tsb_stat=False,
#seqInfo=False)

# Perform signature extraction.
#data = sig.importdata("matrix")
if __name__ ==  '__main__':
    
    # Determine location input
    cwd = os.getcwd()
    input_mut_mat = cwd +"/simulated_data/mut_mat_" + args.name + ".txt"
    
    # Determine output location
    output_dir = cwd + "/sigprofiler_output_" + args.name
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    
    
    sig.sigProfilerExtractor("matrix",
    output_dir,
    input_mut_mat,
    context_type = "96",
    resample = True,
    nmf_replicates = 2,
    minimum_signatures = 5,
    maximum_signatures = 20,
    make_decomposition_plots = True,
    refit_denovo_signatures = True,
    get_all_signature_matrices = False,
    export_probabilities = True,
    cpu = -1)
