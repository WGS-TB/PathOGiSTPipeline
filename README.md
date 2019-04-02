# PathOGiST

## Build Status
[![Build Status](https://travis-ci.org/WGS-TB/PathOGiST.svg?branch=master)](https://travis-ci.org/WGS-TB/PathOGiST)

## Bugfixing Protocol
1. Raise Github issue
2. Fix the issue
3. Create new unit tests 
4. Run unit tests
5. Close issue

## Tasks
- [ ] Lawyering up
  - [x] Add a license (e.g. MIT)
  - [ ] In each package file, add a header with copyright information
- [ ] Documentation


## Installation
*Note*: PathOGiST is currently not compatible with OSX

We recommend you create a conda environment for PathOGiST, and install PathOGiST through conda.
First set up Bioconda as per the instructions [here](https://bioconda.github.io/).
PathOGiST requires Python 3.5 or newer:
```bash
conda create --name pathogist 
```
And then activate the environment and install PathOGiST:
```bash
source activate pathogist
conda install pathogist
```
When inside the `pathogist` conda environment, you can then simply run `PATHOGIST -h`, for example.
Note that you will need to install CPLEX separately, as CPLEX is proprietary software.

## Subcommands

### Entire Pipeline (`PATHOGIST run`)
This subcommand runs the PathOGiST pipeline from start to finish 
(i.e. distance matrix creation -> correlation clustering -> consensus clustering).

The main input file is a YAML configuration file, which you can create with the command
```bash
PATHOGIST run [path to where you want your config] --new_config
```
The configuration file will look like so:
```bash
---
# PathOGiST configuration file.
# This configuration file is in YAML file format, release 1.2 (Third Edition).
# Google yaml for the specification.
# Do not remove the sections 'genotyping', 'distances', 'thresholds', 'all_constraints', or 'output'.
# Remove all key-value pairs in the sections 'genotyping' or 'distances' if you want them blank.
# Remove all list items in the section 'fine_clusterings' if you want that blank, too.
# Keys from 'genotyping' and 'distances' sections should not overlap.
genotyping:
    # Raw genotyping data from which to create distance matrices.
    # Accepted values are paths to text files containing paths to call files.
    # Currently only compatible with SNPs from snippy, MLSTs from MentaLiST, and CNVs from PRINCE.
    # Example:
    # MLST: /path/to/mlst_calls.txt
    SNP:
    MLST:
    CNV:
distances:
    # Paths to pre-constructed distance matrices in tsv format.
    # You can also specify SNP, MLST and CNV distance matrices here if you pre-constructed them,
    # but then they shouldn't appear in the section 'genotyping'.
    kmer:
    spoligotyping:
fine_clusterings:
    # The genotyping datatypes that are considered to be the "finest".
    - MLST
    - SNP
    - kmer
thresholds:
    # Threshold values for performing correlation clustering on genotyping data types given above.
    # Every key appearing in the sections 'genotyping' and 'distances' should appear here with a value. 
    SNP:
    MLST:
    CNV:
    kmer:
    spoligotyping:
# Put the path for the final output consensus clustering (which will be in tsv format) here.
output:
# Use all constraints when performing correlation and consensus clustering
all_constraints: False
# Solver to use for the clustering algorithm; choices are `pulp` or `cplex`
solver: pulp
...
```
Modify the configuration by adding the paths to your input files.
You can add your own keys to the YAML configuration file, and delete the default keys which aren't relevant to your experiment.

The inputs to the `genotyping` entries should be a file which contains absolute paths to your call files.
For example, `mlst_calls.txt` should look something like:
```bash
/absolute/path/to/SRR00001.calls
/absolute/path/to/SRR00002.calls
/absolute/path/to/SRR00003.calls
```
The output of PathOGiST is a TSV file containing the file consensus cluster assignment for each sample.

### Correlation Clustering (`PATHOGIST correlation`)
This subcommand is for clustering bacterial samples based on a distance matrix.

The inputs to correlation clustering are:
* A distance matrix in the form of a TSV file
* A threshold cutoff value for the construction of the similarity matrix
The output is a TSV file containing the cluster assignments of the samples described by the distance matrix.

You can run correlation clustering with the following command:
```bash
PATHOGIST correlation [distance matrix] [threshold] [output path]
```

### Distance Matrix Creation (`PATHOGIST distance`)
This subcommand is used for creating distance matrices from genotyping calls, e.g. SNPs, MLSTs, CNVs, etc.
Currently, this subcommand is only compatible with SNP calls from Snippy, MLST calls from MentaLiST, and CNV calls from Prince.
The input is:
* A text file containing paths to genotyping call files.

The output is a distance matrix represented as a TSV file.

You can run this subcommand like so:
```bash
PATHOGIST distance [path/to/calls_file.tsv] [one of SNP/MLST/CNV] [output path]
```

### Consensus Clustering (`PATHOGIST consensus`)
The input for consensus clustering is three files:
* A text file containing paths to distance matrices in `.tsv` format.
* A text file containing paths to clustering assignments in `.tsv` format.
* A text file containing the names of the clusterings which are 'finest'.

The output is a TSV file containing the cluster assignments of the samples which are common to all the input distance matrices.

You can run consensus clustering with the following command:
```bash
PATHOGIST consensus [distances] [clusterings] [fine_clusterings] [output path]
```

Each line of the input files should correspond to a specific data type, e.g. SNPs, MLSTs, or CNVs.
Absolute paths to distance matrices and cluster assignments should be prepended with the name of the clustering and an equal sign, i.e. `[name]=[absolute path to file]`.
An example:

_Distances file_
```bash
SNP=/path/to/snp_dist
MLST=/path/to/mlst_dist
CNV=/path/to/cnv_dist
```
_Clusterings file_
```bash
SNP=/path/to/snp_clust
MLST=/path/to/mlst_clust
CNV=/path/to/cnv_clust
```
_Fine clusterings file_
```bash
SNP
```
