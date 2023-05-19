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
	parser.add_argument("founders", help="VCF with founder genotypes", 
		     metavar="FOUNDER VCF")

	parser.add_argument("-c", "--chr", help="Chromosome name", metavar="CHR")
	parser.add_argument("-g", "--groups", help="Founder groups", 
		     metavar="GROUP", nargs="+")

	parser.add_argument("-p", "--pca", help="Run PCA on founders", 
		     metavar="NUM", nargs='?', const=1, type=int)
	parser.add_argument("-m", "--matrix", 
		     help="Calculate distance matrix for founders", action="store_true")
	parser.add_argument("-d", "--descendents", metavar="DESCENDENT VCF", 
		     help="VCF with descendent genotypes")

	# Output
	parser.add_argument("-o", "--out", help="Write output to file. " \
		"Default: stdout", metavar="FILE")
	
	parser.add_argument("--version", help="Print the version and quit", \
		action="version", version = '{version}'.format(version=__version__))

	# Parse args
	args = parser.parse_args()

	if sum(map(bool, [args.pca is not None, args.matrix, args.descendents])) != 1:
		myutils.ERROR("Please specify exactly one of --pca, --matrix, or --descendents.")

	if not os.path.exists(args.founders):
		myutils.ERROR("{founders} does not exist".format(founders=args.founders))
	
	if args.descendents is not None and not os.path.exists(args.descendents):
		myutils.ERROR("{desc} does not exist".format(desc=args.descendents))
	
	if args.out is not None and not os.path.exists(os.path.dirname(args.out)):
		myutils.ERROR("Directory for {out} does not exist".format(out=args.out))

	if args.groups is not None and args.pca is not None:
		myutils.ERROR("Groups cannot be used in conjuction with PCA analysis")
	
	if args.groups is not None:
		myutils.ERROR("Not yet implemented")
	
	# Set up output file
	if args.out is None:
		outf = sys.stdout
	else: outf = open(args.out, "w")
	
	if args.matrix:
		founders = myutils.extract_genotypes(args.founders, args.chr)[0]
		outf.write(str(myutils.dist_matrix(founders, None)))
	if args.descendents is not None:
		outf.write(myutils.identify_founders(
			args.founders, args.descendents, args.chr).to_string())
	if args.pca is not None:
		evectors, evalues = myutils.pca(args.founders, args.chr, args.pca)
		outf.write(' '.join(str(round(eval, ndigits = 4)) for eval in evalues))
		outf.write("\n" + str(evectors))
	outf.close()
	sys.exit(0)

if __name__ == "__main__":
    main()