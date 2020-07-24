---
description: Generate a pileup and noise QC metrics from a Bam file
---

# Noise Calculation

## Description

Noise is calculated by looking at positions from the `bed_file` for which there is no alt allele that exceeds `threshold`, and dividing the total number of non-reference bases by the total number of bases at such positions. 

## Usage:

```text
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

