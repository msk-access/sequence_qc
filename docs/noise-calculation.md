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
* `noise.html` - HTML report with summary of
  * Top noisy positions with highest alt allele frequencies
  * Histogram of positions from `bed_file` with each count of masked "N" bases

## Calculation Details

For the overall noise level of the sample, a single valued is calculated over the regions listed in the `bed_file` in the following manner:

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

Essentially, this means for every position which does not have any alt allele which exceeds the threshold, the noise level is the total count of alt bases as such positions, divided by the total number of bases at such positions. For the noise-by-substitution calculation, only the specified alternate allele contributes to the numerator and denominator of the noise fraction. 

## Plots Output

The module will output an HTML report with noise metrics, here is an example report file:

{% file src=".gitbook/assets/donor6-t\_noise.html" caption="sequence\_qc Noise.html Report \(v0.1.19\)" %}

### Noise By Substitution:

![](.gitbook/assets/screen-shot-2020-09-25-at-2.59.33-pm.png)

* Uses the previously-defined calculation to express the noise level of this sample, for each of the 12 possible substitution types
* Expected noise level is on the order of 10e-6 for ACCESS duplex and simplex samples
* C&gt;T and G&gt;A noise levels are usually the highest

### Top noisy positions:

![](.gitbook/assets/screen-shot-2020-09-25-at-2.59.43-pm.png)

* The positions from the `bed_file` are sorted by those with the highest noise fraction, and the top positions' noise fractions are plotted 
* Violin plot represents all positions from the bed file, and is expected to have most positions on the low end, with some outliers closer to the supplied threshold

### Fragment Size distribution for noisy positions:

![](.gitbook/assets/screen-shot-2021-05-13-at-12.17.11-pm.png)

* A histogram of fragment sizes is plotted for reads that contain a "noisy" position \(as defined previously\)
* Substitution types can be plotted individually by clicking the legend

### N Counts Histogram:

![](.gitbook/assets/screen-shot-2020-09-25-at-2.59.54-pm%20%281%29.png)

* Each position is counted for "N" or no-calls, and the number of positions with each N count is plotted as a histogram
* ACCESS samples are expected to have a peak below 10 N's, although duplex and simplex samples will have a larger number of N bases than the original uncollapsed or "standard" bam files

