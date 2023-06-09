# INCH
Identification of Nonrecombinant Chromosome Haplotypes (CSE 185 Project)

This is a class project which aims to match samples to their most likely founder
haplotype. It assumes no recombination, thus, only the Y and MT chromosomes
should be used. Founders must have known genotypes.

This README is pattened after the example one, which may be found at
https://github.com/gymreklab/cse185-demo-project/blob/main/README.md

## Install instructions

Installation requires the `numpy`, `pandas`, `scikit-lean`, and `scikit-allel`
libraries to be installed. You can install these  with `pip`:

```
pip install numpy pandas scikit-learn scikit-allel
```

Next, create a local copy of this repository:
```
git clone https://github.com/Faith-Okamoto/INCH
cd INCH
```

Once required libraries are installed, you can install INCH with the following
command:

```
python setup.py build
python setup.py install
```

Note: if you do not have root access, you can run the commands above with
additional options to install locally:

```
pip install --user numpy pandas scikit-learn scikit-allel
git clone https://github.com/Faith-Okamoto/INCH.git
cd INCH
python setup.py build
python setup.py install --user
```

If the install was successful, typing `inch --help` should show a useful 
message. If that does nothing, check where INCH was installed to, and make sure
that it is on your `PATH` environment variable (run `echo $PATH`).

Code has been tested with these package versions:
- `numpy` 1.24.3
- `pandas` 1.5.3
- `scikit-learn` 1.2.2
- `scikit-allel` 1.3.6

## Basic usage

The basic usage of INCH is:

```
inch {--pca 2|--matrix|--descendents descendents.vcf} founders.vcf [other options]
```

To run INCH on a small test example (using files in this repo):

```
inch -d example-files/descendents.vcf example-files/founders.vcf
```

This should produce the output below:

```
D1      F1
D2      F2
D3      F1
D4      F2
D5      F1
```

Some other commands which may be instructive to run:

```
# distance between founders: they do not match
inch -m example-files/founders.vcf 
# PCA on founders: the first dimension can capture everything
inch -p 2 example-files/founders.vcf 
# distance between descendents: see the groups
inch -m example-files/descendents.vcf
# PCA on descendents: see the groups again
inch -p 3 example-files/descendents.vcf
```

## INCH options

The only required input to INCH is a VCF file. This file may be optionally 
compressed using gzip. Users must specify exactly one of these three flags:
- `-p NUM`, `--pca NUM`: compute NUM number of PCs for the given founders.
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
  - This flag may not be used with the `-p` flag.
- `--dump-matrix FILE`: Write the descendents-vs-founders distance matrix
  produced as an intermediate step to a file. Only used with `--descendents`.

## File format

File formats differ depending on which analysis was requested:
- The format for the `-d` option is a two-column TSV file without a header line.
  The first column is descendent IDs and the second is the ID of the founder or
  group that the descendent matches best.
- The format for the `-m` option is an `n x n` TSV matrix with row and column
  labels of founder IDs, where `n` is the number of founders or founder groups.
  Each cell is the Hamming distance between the founders in that row and column.
- The format for the `-p` option is a file with two parts. The first line has
  tab-separated eigenvalues for the PCA eigenvectors, starting with PC1. After
  that is an `n x p` TSV table, where `n` is the number of founders and `p` is
  the number of PCs used. Column and row labels are included. Each cell is the
  weight of that row's founder under that column's principle component.

## Usage on real data

To run on a larger test example, I have provided some real data from the Y
chromosome of [HS rats](https://ratgenes.org/cores/core-b/). This is courtesy of
[Palmer Lab](https://palmerlab.org/). Get these files with:

```
curl https://palmerlab.s3.sdsc.edu/Faith_CS185_data/data.tar.gz > data.tar.gz
tar -xzvf data.tar.gz
```

A new directory `project-data` will be created. Some suggested commands are:

```
# see relationships between founders
inch -m project-data/founders.vcf.gz
# separate into the two major Y chromosome groups
inch -d project-data/descendents.vcf.gz project-data/founders.vcf.gz -g ACI,BN,MR BUF,F344,M520,WKY,WN
```

## Contributors

This project was entirely coded by Faith Okamoto, with some idea generation by
other members of the Palmer Lab. (Notably Thiago Sanches.) The founder 
identification method was run by Daniel Munro.

## Testing

There is a test suite in `inch/test_utils.py`. It uses `pytest`.