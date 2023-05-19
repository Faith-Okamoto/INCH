import sys
import pandas as pd
from sklearn.metrics import pairwise_distances
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

def dist_matrix(founders, desc, outf, chr):
    founders, pos = extract_genotypes_vcf(founders, chr)
    samples = founders.columns.tolist()
    founders = founders.transpose()
    matrix = pairwise_distances(founders, desc, metric="hamming")
    matrix = pd.DataFrame(matrix, columns = samples)
    matrix.index = samples
    outf.write(str(matrix))

def is_gz_file(filename):
    with open(filename, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

def open_vcf(filename):
    if is_gz_file(filename):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "rt")

def get_vcf_names(filename):
    ifile = open_vcf(filename)
    names = None
    for line in ifile:
        if line.startswith("#CHROM"):
            names = [x for x in line.split()]
            break
    ifile.close()
    if names is None:
        ERROR("VCF file {file} has no header line".format(file=filename))
    return names

def load_vcf(filename, chr):
    vcf = pd.read_csv(filename, comment='#', delim_whitespace = True, 
                      header = None, names = get_vcf_names(filename))
    if vcf.columns[:9].tolist() != STANDARD_VCF_COLS:
        ERROR("VCF column names are malformed")
    if chr is None and vcf['#CHROM'].nunique() > 1:
        ERROR("More than one chromosome detected in VCF file")
    if chr is not None:
        vcf = vcf[vcf['#CHROM'] == chr]
        if vcf.shape[0] == 0:
            ERROR("No variants on chromosome {chr}".format(chr = chr))
    return vcf

def extract_genotypes_vcf(filename, chr):
    vcf = load_vcf(filename, chr)
    if not (vcf["FORMAT"].str.split(':').str[0] == "GT").all():
        ERROR("GT is not the first FORMAT field for all positions")
    samples = vcf.columns[9:]
    alleles = (vcf["REF"] + "," + vcf["ALT"]).str.split(",")
    if not all([set(row).issubset(BASES.keys()) for row in alleles]):
        ERROR("Invalid base in REF or ALT alleles")

    def get_allele(row, geno):
        if geno != "." and not geno.isdigit():
            ERROR("Invalid genotype code: " + geno)
        return -1 if geno == "." else BASES[alleles[row][int(geno)]]
    
    non_haploid = False
    for sample in samples:
        geno = vcf[sample].str.split(':').str[0]
        if geno.str.contains('|').any() or geno.str.contains('/').any():
            non_haploid = True
            geno = geno.str.split('|').str[0].str.split('/').str[0]
        vcf[sample] = [get_allele(i, geno[i]) for i in range(vcf.shape[0])]
    if non_haploid:
        print("Non haploid genotypes detected. First allele used.")

    return vcf.loc[:, samples.tolist()], vcf["POS"]
    
if __name__ == "__main__":
    dist_matrix("example-files/descendents.vcf", None, sys.stdout, None)
