#!/usr/bin/env python

"""
Command-line script to categorize samples one of a set of given haplotypes
"""

# command-line argument handling
import argparse
# helpers from elsewhere in the module
from inch import myutils, __version__
# basic utilities
from os import path
import sys

def main():
	parser = argparse.ArgumentParser(
		prog = 'inch',
		description = 'Command-line script to categorize samples as a haplotype'
	)

	# input
	parser.add_argument('founders', help = 'VCF with founder genotypes', 
		     metavar = 'FOUNDER VCF')

	# extra analysis options
	parser.add_argument('-c', '--chr', help = 'Chromosome to filter from VCF', 
		     metavar = 'CHR')
	parser.add_argument('-g', '--groups', help = 'Founder groups', 
		     metavar = 'GROUP', nargs = '+')

	# what analysis type to run
	parser.add_argument('-p', '--pca', help = 'Run PCA on founders. Default: 2', 
		     metavar = 'NUM', nargs = '?', const = 2, type = int)
	parser.add_argument('-m', '--matrix',
		     help = 'Calculate distance matrix for founders',
			 action = 'store_true')
	parser.add_argument('-d', '--descendents', metavar = 'DESCENDENT VCF', 
		     help = 'VCF with descendent genotypes')

	# output
	parser.add_argument('-o', '--out',
		     help = 'Write output to file. Default: stdout', metavar = 'FILE')
	
	parser.add_argument('--version', help = 'Print the version and quit',
		action='version', version = '{version}'.format(version=__version__))

	# parse args
	args = parser.parse_args()

	# check legality of arguments
	if sum(map(bool, [args.pca is not None, args.matrix, args.descendents])) != 1:
		myutils.ERROR('Please specify exactly one of -p, -m, or -d.')
	if not path.exists(args.founders):
		myutils.ERROR('{founders} does not exist'.format(args.founders))
	if args.descendents is not None and not path.exists(args.descendents):
		myutils.ERROR('{desc} does not exist'.format(args.descendents))
	if args.out is not None and not path.exists(path.dirname(args.out)):
		myutils.ERROR('Directory for {out} does not exist'.format(args.out))
	if args.groups is not None and args.pca is not None:
		myutils.ERROR('Groups cannot be used in conjuction with PCA analysis')
	
	
	# perform appropriate analysis and save result as a string
	if args.matrix:
		out = str(myutils.dist_matrix(args.founders, args.chr, args.groups))
	if args.descendents is not None:
		out = myutils.identify_founders(
			args.founders, args.descendents, args.chr, args.groups).to_string()
	if args.pca is not None:
		e_vecs, e_vals = myutils.pca(args.founders, args.chr, args.pca)
		out = ' '.join([str(round(e_val, ndigits = 4)) for e_val in e_vals]) \
			+ '\n' + str(e_vecs)
	
	# write string to output file
	outf = sys.stdout if args.out is None else open(args.out, 'w')
	outf.write(out)
	outf.close()

	sys.exit(0)

if __name__ == '__main__':
    main()