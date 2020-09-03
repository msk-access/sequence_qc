---
description: Generate a pileup and noise QC metrics from a Bam file
---

# Noise Calculation

## Description

Noise is calculated by looking at positions from the `bed_file` , and setting the genotype for each position to the base with the highest count. Then, for positions where there is no alternate allele that exceeds `threshold`, we divide the total number of non-genotype bases by the total number of bases.

## Usage

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
                             
  --sample_id TEXT           Prefix to include in all output file names 

  --threshold FLOAT          Alt allele frequency past which to ignore
                             positions from the calculation

  --truncate INTEGER         Whether to exclude trailing bases from reads that
                             only partially overlap the bed file (0 or 1)

  --min_mapq INTEGER         Exclude reads with a lower mapping quality
  --min_basq INTEGER         Exclude bases with a lower base quality
  --help                     Show this message and exit.
```

## Parameters

| Parameter | Description | Default |
| :--- | :--- | :--- |
| **ref\_fasta** \(string\) | Path to reference fasta which was used for mapping Bam |  |
| **bam\_file** \(string\) | Path to Bam file for which to do calculation |  |
| **output\_prefix** \(string\) | Prefix used for output files \(normally a sample ID\) |  |
| **bed\_file** \(string\) | Path to bed file which contains regions for which to calculate noise |  |
| **threshold** \(float\) | This value will be used as a definition of "noisy" positions. For the default of `0.02`this means that only positions with alt alleles at less than 2% allele frequency will contribute to the major\_allele\_count and minor\_allele\_count. | 0.02 |
| **truncate** \(bool\) | If set to 1, bases from reads that only partially overlap the regions in `bed_file` will be included in the calculation. | True |
| **min\_mapq** \(int\) | Exclude reads with a lower mapping quality | 1 |
| **min\_basq** \(int\) | Exclude reads with a lower base quality | 1 |

## Outputs Description

* `pileup.tsv` Pileup file of all positions listed in the bed file
* `noise_positions.tsv` Pileup file limited to positions with at least one alt allele below the noise threshold
* `noise_acgt.tsv` Noise file with the following columns \(calculated from single base changes, excluding N and deletions\):

| Column | Description |
| :--- | :--- |
| sample\_id | Taken from `sample_id` parameter |
| minor\_allele\_count | Total number of bases below threshold that do not match the sample's genotype |
| major\_allele\_count | Total number of bases from positions that meet threshold criteria that support the sample's genotype |
| noise\_fraction | `minor_allele_count` divided by `major_allele_count` |
| contributing\_sites | Number of unique sites that contributed to the `minor_allele_count` |

* `noise_n.tsv` This file is identical to `noise_acgt`, however in this case N is used as the minor\_allele and other base changes are ignored
* `noise_del.tsv` This file is identical to `noise_acgt`, however in this case deletions are used as the minor\_allele and other base changes are ignored

## Calculation Details

A single value for the noise level of the sample over the regions listed in the `bed_file` is calculated in the following manner:

$$
\begin{aligned}
&total\_depth_{i} = \sum_{n\ in\ {A, C, G, T}}count(n)\ at\ position\ i\\
\\
&genotype_{i} = max\{count(A), count(C), count(G), count(T)\}\ at\ position\ i\\
\\
&alt\_count_i = \sum_{n\ in\ {A,C,G,T}}{^{0\ if\ n\ =\ genotype_i}_{count(n)\ o.w.}}\\
\\
&noise = 100 \cdot \frac{\sum_j{alt\_count_j}}{\sum_j{total\_depth_j}}\\
\\
&where\ j = positions\ for\ which\ \frac{alt\_count^n_j}{total\_depth_j} < threshold\ for\ n\ in\ {\{A,C,G,T\}}\\
\end{aligned}\\
$$

Essentially, this means for every position which does not have any alt allele which exceeds the threshold, the noise level is the total count of alt bases as such positions, divided by the total number of bases at such positions. 

