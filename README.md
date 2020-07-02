# sequence_qc

Package for doing various ad-hoc quality control steps from MSK-ACCESS generated FASTQ or BAM files

[![Build Status](https://travis-ci.com/msk-access/sequence_qc.svg?branch=master)](https://travis-ci.com/msk-access/sequence_qc)
[![PyPi](https://img.shields.io/pypi/v/sequence_qc.svg?)](https://pypi.python.org/pypi/sequence_qc)
[![Anaconda](https://anaconda.org/ionox0/sequence-qc/badges/version.svg)](https://anaconda.org/ionox0/sequence-qc/)

* Free software: Apache Software License 2.0
* Documentation: https://msk-access.gitbook.io/sequence-qc/

### Installation

From pypi:

``pip install sequence_qc``

From conda:

``conda install -c ionox0 -c conda-forge -c bioconda sequence-qc``

### Usage
```
Usage: calculate_noise [OPTIONS]

  Calculate noise level of given bam file, across the given positions in
  `bed_file`.

Options:
  --ref_fasta TEXT           Path to reference fasta, containing all regions
                             in bed_file  [required]
  --bam_file TEXT            Path to BAM file for calculating noise
                             [required]
  --bed_file TEXT            Path to BED file containing regions over which to
                             calculate noise  [required]
  --threshold FLOAT          Alt allele frequency past which to ignore
                             positions from the calculation
  --truncate INTEGER         Whether to exclude trailing bases from reads that
                             only partially overlap the bed file (0 or 1)
  --min_mapq INTEGER         Exclude reads with a lower mapping quality
  --min_basq INTEGER         Exclude bases with a lower base quality
  --help                     Show this message and exit.
```

### Versioning

To increase the version number use the following command:

``bumpversion (major|minor|patch) --tag``
