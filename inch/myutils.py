import sys
import pandas as pd
from sklearn.metrics import pairwise_distances
from sklearn.decomposition import PCA
from itertools import chain
import gzip
# error handling code from https://github.com/gymreklab/cse185-demo-project/blob/main/mypileup/myutils.py

STANDARD_VCF_COLS = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 
                     'QUAL', 'FILTER', 'INFO', 'FORMAT']
BASES = {"A" : 0, "C": 1, "G": 2, "T": 3}

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def ERROR(msg):
    """
    Print an error message and die

    Parameters
    ----------
    msg : str
        Error message to print
    """
    sys.stderr.write(bcolors.FAIL + "[ERROR]: " + bcolors.ENDC + "{msg}\n".format(msg=msg) )
    sys.exit(1)

def make_groups(founder_ids, groups):
    if groups is None:
        return None
    groups = [group.split(",") for group in groups]
    if not all([set(group).issubset(founder_ids) for group in groups]):
        ERROR("Some groups contain nonexistant founder IDs")
    flat_groups = list(chain.from_iterable(groups))
    if len(flat_groups) != len(set(flat_groups)):
        ERROR("Some groups contain duplicate founder IDs")
    for id in founder_ids:
        if id not in flat_groups:
            groups.append([id])
    return groups

def identify_founders(founders, desc, chr, groups):
    founders, founder_pos, founder_chr = extract_genotypes(founders, chr)
    if groups is not None:
        groups = make_groups(set(founders.columns), groups)
        groups = {id : ",".join(group) for group in groups for id in group}
    desc, desc_pos, desc_chr = extract_genotypes(desc, chr)
    if desc_chr != founder_chr:
        ERROR("Founder and descendents have different chromosomes")
    if not set(founder_pos).intersection(set(desc_pos)):
        ERROR("Founders and descendents share no positions")
    founders = founders[founder_pos.isin(desc_pos)]
    desc = desc[desc_pos.isin(founder_pos)]
    matrix = dist_matrix(founders, desc, None).transpose()
    matches = matrix.idxmin(axis=1)
    return matches if groups is None else matches.replace(groups)

def dist_matrix(founders, desc, groups):
    groups = make_groups(set(founders.columns), groups)
    founder_ids = founders.columns.tolist()
    founders = founders.transpose()
    if desc is not None:
        desc_ids = desc.columns.tolist()
        desc = desc.transpose()
    matrix = pairwise_distances(founders, desc, metric="hamming")
    matrix = pd.DataFrame(matrix, index = founder_ids,
                          columns = founder_ids if desc is None else desc_ids)
    return matrix if groups is None else merge_groups(matrix, groups)

def merge_groups(matrix, groups):
    groups = {",".join(group) : group for group in groups}
    merged = pd.DataFrame(0, index = groups.keys(), columns = groups.keys())
    def dist(row, col):
        return 0 if row == col else matrix.loc[groups[row]][groups[col]].mean(axis=None)
    return merged.apply(lambda x: [dist(row, x.name) for row in x.index])

def is_gz_file(file):
    with open(file, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

def get_vcf_names(file):
    ifile = gzip.open(file, "rt") if is_gz_file(file) else open(file, "rt")
    names = None
    for line in ifile:
        if line.startswith("#CHROM"):
            names = [x for x in line.split()]
            break
    ifile.close()
    if names is None:
        ERROR("VCF file {name} has no header line".format(name=file))
    return names

def load_vcf(file, chr):
    vcf = pd.read_csv(file, comment='#', delim_whitespace = True, 
                      header = None, names = get_vcf_names(file))
    if vcf.columns[:9].tolist() != STANDARD_VCF_COLS:
        ERROR("VCF column names are malformed")
    if chr is None and vcf['#CHROM'].nunique() > 1:
        ERROR("More than one chromosome detected in VCF file")
    if chr is not None:
        vcf = vcf[vcf['#CHROM'] == chr]
        if vcf.shape[0] == 0:
            ERROR("No variants on chromosome {chr}".format(chr=chr))
    return vcf

def to_code(geno):
    if geno == "*":
        return -1
    code = 0
    for base in geno:
        if base not in BASES:
            ERROR("Invalid base {base} in REF or ALT alleles".format(base))
        code *= len(BASES)
        code += BASES[base]
    return code

def extract_genotypes(file, chr):
    vcf = load_vcf(file, chr)
    if not (vcf["FORMAT"].str.split(':').str[0] == "GT").all():
        ERROR("GT is not the first FORMAT field for all positions")
    samples = vcf.columns[9:]
    alleles = (vcf["REF"] + "," + vcf["ALT"]).str.split(",")
    alleles = [[to_code(geno) for geno in row] for row in alleles]

    def get_geno_code(row, col, geno):
        if geno != "." and geno != "*" and not geno.isdigit():
            ERROR("Invalid genotype code: " + geno)       
        return -1 - col if geno == "." else alleles[row][int(geno)]
    
    non_haploid = False
    for col in range(9, vcf.shape[1]):
        geno = vcf.iloc[:, col].str.split(':').str[0]
        if geno.str.contains('|').any() or geno.str.contains('/').any():
            non_haploid = True
            geno = geno.str.split('|').str[0].str.split('/').str[0]
        vcf.iloc[:, col] = [get_geno_code(row, col, geno[row]) for row in range(vcf.shape[0])]
    if non_haploid:
        sys.stderr.write("Non haploid genotypes detected in {file}. First allele used.".format(file = file))

    return vcf.loc[:, samples.tolist()], vcf["POS"], vcf["#CHROM"][0]

def pca(file, chr, ncomp):
    geno = extract_genotypes(file, chr)[0]
    samples = geno.columns
    geno = geno.transpose()
    pca = PCA(n_components=ncomp)
    pca.fit(geno)
    eigenvec = pd.DataFrame(pca.transform(geno))
    eigenvec.columns = ["PC" + str(i + 1) for i in range(ncomp)]
    eigenvec.index = samples
    return eigenvec, pca.explained_variance_