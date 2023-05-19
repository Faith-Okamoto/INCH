import sys
import pandas as pd
import sklearn as sk
import gzip
# error handling code from https://github.com/gymreklab/cse185-demo-project/blob/main/mypileup/myutils.py

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

def dist_matrix(founders, desc, outf):
    founders = parse_vcf(founders)
    print(founders)
    #matrix = sk.metrics.pairwise_distances(founders, desc)
    print("hi")
    #print(matrix)
    #outf.write(matrix)

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
    for line in ifile:
        if line.startswith("#CHROM"):
            names = [x for x in line.split('\t')]
            break
    ifile.close()
    return names

def parse_vcf(filename):
    return pd.read_csv(filename, comment='#', delim_whitespace = True, 
                       header = None, names = get_vcf_names(filename))
    
