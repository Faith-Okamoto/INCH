#!/usr/bin/env python

"""
Command-line script to categorize samples one of a set of given haplotypes
"""

import argparse
from inch import myutils
from inch import __version__
import os
import sys

def main():
	parser = argparse.ArgumentParser(
		prog="inch",
		description="Command-line script to categorize samples as a haplotype"
	)

	# Input
	parser.add_argument("founders", help="VCF with founder genotypes", metavar="FOUNDER VCF")

	parser.add_argument("-c", "--chr", help="Chromosome name", metavar="CHR")
	parser.add_argument("-g", "--groups", help="Founder groups", metavar="GROUP", nargs="+")

	parser.add_argument("-p", "--pca", help="Run PCA on founders", nargs='?', const=1, type=int, metavar="NUM")
	parser.add_argument("-m", "--matrix", help="Calculate distance matrix for founders", action="store_true")
	parser.add_argument("-d", "--descendents", metavar="DESCENDENT VCF", help="VCF with descendent genotypes")

	# Output
	parser.add_argument("-o", "--out", help="Write output to file. " \
		"Default: stdout", metavar="FILE")
	
	parser.add_argument("--version", help="Print the version and quit", \
		action="version", version = '{version}'.format(version=__version__))

	# Parse args
	args = parser.parse_args()

	if sum(map(bool, [args.pca is not None, args.matrix, args.descendents])) != 1:
		myutils.ERROR("Please specify exactly one of --pca, --matrix, or --descendents.")
	
	if args.pca is not None and args.out is None:
		myutils.ERROR("The --pca option must be used with --out.")

	if not os.path.exists(args.founders):
		myutils.ERROR("{founders} does not exist".format(founders=args.founders))
	
	if args.descendents is not None:
		if not os.path.exists(args.descendents):
			myutils.ERROR("{descendents} does not exist".format(descendents=args.descendents))
		desc = args.descendents
	else:
		desc = None
	
	

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