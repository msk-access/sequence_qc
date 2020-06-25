import logging
import pysamstats

import pandas as pd

from pysam import AlignmentFile
from pybedtools import BedTool

FORMAT = '%(asctime)-15s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger("sequence_qc")
logger.setLevel(logging.DEBUG)

EPSILON = 1e-9
OUTPUT_PILEUP_NAME = 'pileup.tsv'
OUTPUT_NOISE_FILENAME = 'noise_positions.tsv'

output_columns = [
    'chrom',
    'pos',
    'ref',
    'A',
    'C',
    'G',
    'T',
    'insertions',
    'deletions',
    'N',
]


def calculate_noise(ref_fasta: str, bam_path: str, bed_file_path: str, noise_threshold: float,
                    noise_output_filename: str = OUTPUT_NOISE_FILENAME, truncate: bool = True,
                    min_mapping_quality: int = 1, min_base_quality: int = 1, output_prefix: str = '',):
    """
    Create file of noise across specified regions in `bed_file` using pybedtools and pysamstats

    :param ref_fasta: string - path to reference fastq
    :param bam_path: string - path to bam
    :param bed_file_path: string - path to bed file
    :param output_prefix: string - prefix for output files
    :param noise_threshold: float - threshold past which to exclude positions from noise calculation
    :param noise_output_filename: string - filename to give output pileup
    :param truncate: int - 0 or 1, whether to exclude reads that only partially overlap the bedfile
    :param min_mapping_quality: int - exclude reads with mapping qualities less than this threshold
    :param min_base_quality: int - exclude bases with less than this base quality
    :return:
    """
    bed_file = BedTool(bed_file_path)
    bam = AlignmentFile(bam_path)
    pileup_df_all = pd.DataFrame()

    # Build data frame of all positions in bed file
    for region in bed_file.intervals:
        chrom = region.chrom.replace('chr', '')
        start = region.start
        stop = region.stop

        pileup = pysamstats.load_pileup('variation', bam, chrom=chrom, start=start, end=stop, fafile=ref_fasta,
                                        truncate=truncate, max_depth=30000, min_baseq=min_base_quality,
                                        min_mapq=min_mapping_quality)

        pileup_df_all = pd.concat([pileup_df_all, pd.DataFrame(pileup)])

    # Convert bytes objects to strings so output tsv is formatted correctly
    for field in ['chrom', 'ref']:
        pileup_df_all.loc[:, field] = pileup_df_all[field].apply(lambda s: s.decode('utf-8'))

    # Save the complete pileup
    pileup_df_all[output_columns].to_csv(output_prefix + OUTPUT_PILEUP_NAME, sep='\t', index=False)

    # Determine per-position genotype and alt count
    pileup_df_all = _calculate_alt_and_geno(pileup_df_all)

    # Include columns for ins / dels / N
    pileup_df_all = _include_indels_and_n_noise(pileup_df_all)

    # Filter to only positions below noise threshold
    thresh_boolv = pileup_df_all.apply(_apply_threshold, axis=1, thresh=noise_threshold)
    below_thresh_positions = pileup_df_all[thresh_boolv]

    # Filter again to positions with noise
    noisy_boolv = (below_thresh_positions['mismatches'] > 0) | \
                  (below_thresh_positions['insertions'] > 0) | \
                  (below_thresh_positions['deletions'] > 0) | \
                  (below_thresh_positions['N'] > 0)

    noisy_positions = below_thresh_positions[noisy_boolv]
    noisy_positions = noisy_positions.sort_values('alt_count')
    relevant_cols = output_columns + ['alt_count', 'geno_count']
    noisy_positions[relevant_cols].to_csv(output_prefix + noise_output_filename, sep='\t', index=False)

    # Calculate sample noise
    alt_count_total = below_thresh_positions['alt_count'].sum()
    geno_count_total = below_thresh_positions['geno_count'].sum()
    noise = alt_count_total / (alt_count_total + geno_count_total + EPSILON)

    logger.info('Alt count, Geno count, Noise: {} {} {}'.format(alt_count_total, geno_count_total, noise))
    return noise


def _apply_threshold(row: pd.Series, thresh: float) -> bool:
    """
    Returns False if any alt allele crosses `thresh` for the given row of the pileup, True otherwise

    :param row: pandas.Series - row that represents single pileup position
    :param thresh: float - threshold past which alt allele fraction should return false
    """
    base_counts = {'A': row['A'], 'C': row['C'], 'G': row['G'], 'T': row['T']}
    genotype = max(base_counts, key=base_counts.get)
    non_geno_bases = ['A', 'C', 'G', 'T']
    non_geno_bases.remove(genotype)

    if any([row[r] / (row['total_acgt'] + EPSILON) > thresh for r in non_geno_bases]):
        return False

    return True


def _calculate_alt_and_geno(noise_df: pd.DataFrame) -> pd.DataFrame:
    """
    Determine the genotype and alt count for each position in the `noise_df`

    :param noise_df: pd.DataFrame
    :return: pd.DataFrame
    """
    noise_df['total_acgt'] = noise_df['A'] + noise_df['C'] + noise_df['G'] + noise_df['T']
    noise_df['geno_count'] = noise_df[['A', 'C', 'G', 'T']].max(axis=1)
    noise_df['alt_count'] = noise_df['total_acgt'] - noise_df['geno_count']
    return noise_df


def _include_indels_and_n_noise(noise_df: pd.DataFrame) -> pd.DataFrame:
    """
    Add additional columns for noise including insertions / deletions / all indels / N

    :param noise_df: pd.DataFrame from pysamstats.pileup with the following columns:
        [A, C, G, T, insertions, deletions, N]
    :return:
    """
    # 1. Noise including insertions as possible genotype or alt allele
    noise_df['total_acgt_ins'] = noise_df['total_acgt'] + noise_df['insertions']
    noise_df['geno_count_ins'] = noise_df[['A', 'C', 'G', 'T', 'insertions', 'deletions']].max(axis=1)
    noise_df['alt_count_ins'] = noise_df['total_acgt_ins'] - noise_df['geno_count_ins']

    # 2. Noise including deletions as possible genotype or alt allele
    noise_df['total_acgt_del'] = noise_df['total_acgt'] + noise_df['deletions']
    noise_df['geno_count_del'] = noise_df[['A', 'C', 'G', 'T', 'deletions']].max(axis=1)
    noise_df['alt_count_del'] = noise_df['total_acgt_del'] - noise_df['geno_count_del']

    # 3. Noise including insertions or deletions as possible genotype or alt allele
    noise_df['total_acgt_indel'] = noise_df['total_acgt'] + noise_df['insertions'] + noise_df['deletions']
    noise_df['geno_count_indel'] = noise_df[['A', 'C', 'G', 'T', 'insertions', 'deletions']].max(axis=1)
    noise_df['alt_count_indel'] = noise_df['total_acgt_indel'] - noise_df['geno_count_indel']

    # 4. Noise including N bases as alt allele (but won't ever be considered as genotype)
    noise_df['total_acgt_N'] = noise_df['total_acgt'] + noise_df['N']
    noise_df['alt_count_N'] = noise_df['total_acgt_N'] - noise_df['geno_count']

    return noise_df
