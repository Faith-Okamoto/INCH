#!/usr/bin/env python

"""
Command-line script to categorize samples one of a set of given haplotypes
"""

import argparse
from . import myutils as myutils
from inch import __version__
import os
from cyvcf2 import VCF
import sys

def main():
	parser = argparse.ArgumentParser(
		prog="inch",
		description="Command-line script to categorize samples as a haplotype"
	)

	# Input
	parser.add_argument("founders", help="VCF with founder genotypes", type=str)

	parser.add_argument("-c", "--chr", "chr", help="chromosome name", type=str)

	parser.add_argument("-p", "progeny", help="VCF with progeny genotypes", type=str)


	# Output
	parser.add_argument("-o", "--out", help="Write output to file. " \
		"Default: stdout", metavar="FILE", type=str, required=False)
	
	parser.add_argument("--version", help="Print the version and quit", \
		action="version", version = '{version}'.format(version=__version__))

	# Parse args
	args = parser.parse_args()

	# Set up output file
	if args.out is None:
		outf = sys.stdout
	else: outf = open(args.out, "w")
	# https://scikit-learn.org/stable/modules/generated/sklearn.metrics.pairwise_distances.html
	# Peform computation
	outf.close()
	sys.exit(0)

if __name__ == "__main__":
    main()