"""
Utilities for inch.py

Has functions corresponding to each of inch.py's central analysis flags, an
error function, and various helpers to make the analysis functions work.
"""

# basic utilities
import gzip
from typing import Tuple
from itertools import chain
# used to access stderr and force-kill the program
import sys
# used for handling large amounts of data
import pandas as pd
from pandas.api.types import is_numeric_dtype
# a function for computing distance matrices
from sklearn.metrics import pairwise_distances
# a class to do PCA
from sklearn.decomposition import PCA

# these are the first 9 columns in a well-formed VCF file which has samples
STANDARD_VCF_COLS = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 
                     'QUAL', 'FILTER', 'INFO', 'FORMAT']
# a nonintuitive name for the chromosome column; artifact of how file is read
CHR_COL = STANDARD_VCF_COLS[0]
# a simple encoding of DNA bases to numbers
BASES = {'A' : 1, 'C': 2, 'G': 3, 'T': 4}

def ERROR(msg):
    """
    Print an error message and die (copied from Dr. Gymrek's example project)

    Parameters
    ----------
    msg : str
        Error message to print
    """

    sys.stderr.write('[ERROR]: {msg}\n'.format(msg = msg))
    sys.exit(1)

def dist_matrix(founders: str, chr: str, groups: list[str]) -> pd.DataFrame:
    """
    Calculate pairwise distances between all founders, perhaps grouped

    Parameters
    ----------
    founders : str
        VCF file with founder genotypes
    chr : str
        Chromosome to use from the VCF file (None means to use all)
    groups : list[str]
        Founders to group together during distance computation.
        Each list item is a group; IDs with a group are comma-separated.
    
    Returns
	-------
	matrix : pd.DataFrame
        n x n matrix (labeled) of Hamming distances between founders/groups
    """

    geno = _get_geno(founders, chr)[0]
    # process groups early to avoid unnecessary computation if they error
    if groups is not None:
        groups = _make_groups(set(geno.columns), groups)
    # founder v founder distances
    matrix = _geno_dists(geno, geno)
    return matrix if groups is None else _merge_matrix_groups(matrix, groups)

def pca(founders: str, chr: str, n_pc: int) -> Tuple[pd.DataFrame, list[float]]:
    """
    Run PCA on samples

    Parameters
    ----------
    founders : str
        VCF file with founders genotypes
    chr : str
        Chromosome to use from the VCF file (None means to use all)
    n_pc : int
        How many principle components to calculate
    
    Returns
	-------
	eigenvectors : pd.DataFrame
        k samples x n PCs table of each sample's weight along each PC
    eigenvalues : list[float]
        eigenvalues for each PC in decreasing order (PC1, PC2, ...)
    """

    geno = _get_geno(founders, chr)[0]
    # PCA requires the rows to be samples and the columns to be features
    geno = geno.transpose()
    # the PCA model must be pre-initialized with how many components to extract
    pca = PCA(n_components = n_pc)
    pca.fit(geno)
    # extract weights for each sample along the eigenvectors
    eigenvec = pd.DataFrame(pca.transform(geno), index = geno.index,
                            columns = ['PC' + str(i + 1) for i in range(n_pc)])
    return eigenvec, [evalue for evalue in pca.explained_variance_]

def identify_founders(founders: str, descendents: str, 
                      chr: str, groups: list[str]) -> pd.Series:
    """
    Identify which founder a descendent matches best

    Parameters
    ----------
    founders : str
        VCF file with founder genotypes
    descendents : str
        VCF file with descendent genotypes
    chr : str
        Chromosome to use from the VCF files (None means to use all)
    groups : list[str]
        Founders to group together during final assignment.
        Each list item is a group; IDs with a group are comma-separated.
    
    Returns
	-------
	matches : pd.Series
        The founder ID best matching each descendent (descendent IDs in index)
    """

    founders, founder_pos, founder_chr = _get_geno(founders, chr)
    # process groups early to avoid unnecessary computation if they error
    if groups is not None:
        groups = _make_groups(set(founders.columns), groups)
        # dictionary with each founder ID pointing to its group's string
        groups = {id : ','.join(group) for group in groups for id in group}
    desc, desc_pos, desc_chr = _get_geno(descendents, chr)

    if desc_chr != founder_chr:
        ERROR('Founder and descendents have different chromosomes')
    if not set(founder_pos).intersection(set(desc_pos)):
        ERROR('Founders and descendents share no positions')
    
    # filter down to only shared positions
    founders = founders[founder_pos.isin(desc_pos)]
    desc = desc[desc_pos.isin(founder_pos)]
    
    # select closest founder to each desc using all founder v desc distances
    matches = _geno_dists(desc, founders).idxmin(axis = 1)
    return matches if groups is None else matches.replace(groups)

def _make_groups(founder_ids: set[str], groups: list[str]) -> list[list[str]]:
    """
    Process input groups to create final final groups

    Parameters
    ----------
    founder_ids : set[str]
        All founder IDs
    groups : list[str]
        Founders to group together. Must not be None.
        Each list item is a group; IDs with a group are comma-separated.
    
    Returns
	-------
	groups : list[list[str]]
        Sublists are groups and sublist items are founder IDs. 
        All IDs are in exactly one group.
    """

    # turn input groups into a ragged array
    groups = [group.split(',') for group in groups]
    flat_groups = list(chain.from_iterable(groups))

    if len(flat_groups) != len(set(flat_groups)):
        ERROR('Some groups contain duplicate founder IDs')
    if not set(flat_groups).issubset(founder_ids):
        ERROR('Some groups contain nonexistant founder IDs')
    
    # make a new group for each ID which appears in no group
    for id in founder_ids:
        if id not in flat_groups:
            groups.append([id])
    return groups

def _geno_dists(row_geno: pd.DataFrame, col_geno: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate pairwise distances between sample genotypes

    Parameters
    ----------
    row_geno:
        k positions x n_row table of genotypes for samples to be output rows
    col_geno:
        k positions x n_col table of genotypes for samples to be output columns
    
    Returns
	-------
	matrix : pd.DataFrame
        n_row x n_col matrix (labeled) of Hamming distances between samples
    """

    # sklearn requires samples to be rows and features to be columns
    matrix = pairwise_distances(row_geno.transpose(), col_geno.transpose(), 
                                metric = 'hamming')
    # label samples before returning the matrix
    return pd.DataFrame(matrix, index = row_geno.columns, 
                        columns = col_geno.columns)

def _merge_matrix_groups(matrix: pd.DataFrame, 
                         groups: list[list[str]]) -> pd.DataFrame:
    """
    Merge groups in a distance matrix by averaging distance between members

    Parameters
    ----------
    matrix : pd.DataFrame
        n x n matrix (labeled) of Hamming distances between founders
    groups : list[list[str]]
        Sublists are groups and sublist items are founder IDs. Must not be None.
        All IDs are in exactly one group.
    
    Returns
	-------
	merged : pd.DataFrame
        n groups x n groups matrix (labeled) of Hamming distances between groups
    """

    # dictionary with each group's string pointing to its list of IDs
    groups = {','.join(group) : group for group in groups}
    merged = pd.DataFrame(0, index = groups.keys(), columns = groups.keys())

    def dist(row, col):
        # calculate distance at row,col in the merged distance matrix
        # main diagonal is manually set to 0 and not average in-group distance
        return 0 if row == col else \
            matrix.loc[groups[row]][groups[col]].mean(axis = None)
    
    # average between each group's members in the merged distance matrix
    return merged.apply(lambda x: [dist(row, x.name) for row in x.index])

def _is_gz_file(file: str) -> bool:
    """
    Check if a file is gzip-compressed using the magic number

    Parameters
    ----------
    file : str
        filename to check
    
    Returns
	-------
	is_gz : bool
        if this file is gzip compressed
    """

    # check if the first two bytes are the magic number
    with open(file, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

def _get_vcf_names(file: str) -> list[str]:
    """
    Extract column names from a VCF file

    Method from https://www.biostars.org/p/416324/#9480044

    Parameters
    ----------
    file : str
        VCF filename
    
    Returns
	-------
	names : list[str]
        column names from the VCF file
    """

    # open file with normal or gzip reader, as appropriate
    ifile = gzip.open(file, 'rt') if _is_gz_file(file) else open(file, 'rt')
    # placeholder value to allow a later check if any names were found
    names = None

    for line in ifile:
        # the line starting with #CHROM is the one with column names
        if line.startswith('#CHROM'):
            names = [x for x in line.split()]
            break
    ifile.close()

    if names is None:
        ERROR('VCF file {name} has no header line'.format(file = file))
    if names[:9] != STANDARD_VCF_COLS:
        ERROR('VCF column names are malformed; the first 9 are not as expected')
    return names

def _load_vcf(file: str, chr: str) -> pd.DataFrame:
    """
    Load a VCF file into a Pandas DataFrame

    Parameters
    ----------
    file : str
        VCF filename
    chr : str
        Chromosome to use from the VCF file (None means to use all)
    
    Returns
	-------
	vcf : pd.DataFrame
        k positions x 9 + n samples table with all of the VCF's data 
    """

    # read VCF with pandas, splitting on whitespace and ignoring header lines
    vcf = pd.read_csv(file, comment='#', delim_whitespace = True, 
                      header = None, names = _get_vcf_names(file))
    
    if chr is None and vcf[CHR_COL].nunique() > 1:
        ERROR('More than one chromosome detected in VCF file. ' \
              'Must specify one chromosome to use.')
    # subset to a single chromosome if specified
    if chr is not None:
        # index must be reset for later looping over rows
        vcf = vcf[vcf[CHR_COL] == chr].reset_index(drop = True)
        # check if there are now no rows
        if vcf.shape[0] == 0:
            ERROR('No variants on chromosome {chr}'.format(chr = chr))
            
    return vcf

def _to_code(geno: str) -> int:
    """
    Convert a genotype to a number

    Parameters
    ----------
    geno: str
        An allele's genotype (from a VCF file, REF or ALT columns)
    
    Returns
	-------
	code : int
        A numeric code corresponding to this genotype
    """
        
    # missingness due to an upstream deletion
    if geno == '*':
        return -1
    # build up numeric code for allele as if ACGT is 1234 in base-5
    code = 0
    for base in geno:
        if base not in BASES:
            ERROR('Invalid base {base} in REF or ALT'.format(base = base))
        code = code * 5 + BASES[base]
    return code

def _get_geno(file: str, chr: str) -> Tuple[pd.DataFrame, pd.Series, str]:
    """
    Extract unambiguous numeric genotypes from a VCF file

    Parameters
    ----------
    file : str
        VCF filename
    chr : str
        Chromosome to use from the VCF file (None means to use all)
    
    Returns
	-------
	geno : pd.DataFrame
        k positions x n samples table of numeric genotypes
    positions : pd.Series
        list of k positions in order
    chr : str
        Chromosome used from VCF file
    """

    vcf = _load_vcf(file, chr)
    if not (vcf['FORMAT'].str.split(':').str[0] == 'GT').all():
        ERROR('GT must be the first FORMAT field for all positions')

    # load and encode all alleles for each position into a ragged array
    alleles = (vcf['REF'] + ',' + vcf['ALT']).str.split(',')
    alleles = [[_to_code(geno) for geno in row] for row in alleles]

    def get_geno_code(row, col, geno):
        # convert a GT code (relative to REF and ALT) to an unambiguous number
        if geno != '.' and geno != '*' and not geno.isdigit():
            ERROR('Invalid genotype code: ' + geno)       
        # unknown (.) have sample-unique codes; no matching on missingness
        return -1 - col if geno == '.' else alleles[row][int(geno)]
    
    # flag for if any genotypes were detected to be nonhaploid
    non_haploid = False

    # loop over all sample indices in the vcf's columns
    for col in range(9, vcf.shape[1]):
        geno = vcf.iloc[:, col]
        # if a column has only 0/1/etc., it will be numeric, but this needs str
        if is_numeric_dtype(geno):
            geno = geno.astype("string")
        # GT is the first item, before any :
        geno = geno.str.split(':').str[0]
        # | and / are used to separate alleles on different chromosomes
        if geno.str.contains('|').any() or geno.str.contains('/').any():
            non_haploid = True
            # graceful handling of double genotypes for these known-haploid chrs
            geno = geno.str.split('|').str[0].str.split('/').str[0]
        # update this sample's genotypes to have unambiguous numeric codes
        vcf.iloc[:, col] = [get_geno_code(row, col, geno[row]) 
                            for row in range(vcf.shape[0])]
    
    if non_haploid:
        sys.stderr.write('\nNon haploid genotypes detected in {file}. ' \
                         'First allele used.\n'.format(file = file))
    return vcf.iloc[:, 9:], vcf['POS'], vcf[CHR_COL][0]