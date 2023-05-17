# INCH
Identification of Nonrecombinant Chromosome Haplotypes (CSE 185 Project)

This is a class project which aims to match samples to their most likely founder
haplotype. It assumes no recombination, thus, only the Y and MT chromosomes
should be used. Founders must have known genotypes.

This README is pattened after the example one, which may be found at
https://github.com/gymreklab/cse185-demo-project/blob/main/README.md

## Install instructions

Installation requires the (?) libraries to be installed. You can install these 
with `pip`:

```
pip install pyfaidx pysam
```

Once required libraries are installed, you can install INCH with the following
command:

```
python3 setup.py install
```

Note: if you do not have root access, you can run the commands above with
additional options to install locally:

```
pip install --user <libraries>
python3 setup.py install --user
```

If the install was successful, typing `inch --help` should show a useful 
message.

## Basic usage

The basic usage of INCH is:

```
inch [other options] [-c CHR] [--pca/--matrix/-d descendents.vcf] founders.vcf
```

To run INCH on a small test example (using files in this repo):

```
inch -c Y -d example-files/descendents.vcf example-files/founders.vcf
```

This should produce the output below:
```
D1  F1
D2  F2
D3  F1
D4  F2
D5  F1
```

## INCH options

The only required input to INCH is a VCF file. Users must specify exactly one of
the three following flags:
- `-p NUM`, `--pca NUM`: compute NUM number of PCs for the given founders. The 
  `-o` option must also be used; eigenvalue, eigenvector, and image files will
  be written using that prefix. Default NUM is 5.
- `-m`, `--matrix`: compute a distance matrix for the given founders.
- `-d FILE`, `--descendents FILE`: match descendents to founders.

Users may additionally specify the options below:
- `-c CHR`, `--chr CHR`: select only data for a specific chromosome. By default,
assumes all data is from the same chromosome.
- `-o FILE`, `--output FILE`: Write output to file. By default, output is
  written to stdout.
- `-g GROUPS`, `--groups GROUPS`: group multiple founders together. The format
  is `F1,F2 F3,F4,F5` - note the commas separating founders within a group and
  the whitespace separating groups. Any founders not mentioned will be put in a
  group by themselves.
  - If used with the `-d` flag, a descendent will be assigned to the group
    containing the founder it matches best.
  - If used with the `-m` flag, average distances between the members of each 
    group will be computed.

## File format

The format for the `-d` option is a two-column tab-separated file without a
header line. The first column is descendent IDs and the second is the ID of the
founder it matches best.

## Contributors

This project was entirely coded by Faith Okamoto, with some idea generation by
other members of the Palmer Lab. (Notably Thiago Sanches)

## Testing

To run tests (no tests right now)