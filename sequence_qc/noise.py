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

# Output files
NOISE_ACGT = 'noise_acgt.tsv'
NOISE_DEL = 'noise_del.tsv'
NOISE_ACGT_INDEL = 'noise_acgt_indel.tsv'
NOISE_N = 'noise_n.tsv'
# Headers for output files
ALT_COUNT = 'minor_allele_count'
GENO_COUNT = 'major_allele_count'
SAMPLE_ID = 'sample_id'
NOISE_FRACTION = 'noise_fraction'
CONTRIBUTING_SITES = 'contributing_sites'

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


def calculate_noise(ref_fasta: str, bam_path: str, bed_file_path: str, noise_threshold: float, truncate: bool = True,
                    min_mapping_quality: int = 1, min_base_quality: int = 1, output_prefix: str = '',) -> float:
    """
    Create file of noise across specified regions in `bed_file` using pybedtools and pysamstats

    :param ref_fasta: string - path to reference fastq
    :param bam_path: string - path to bam
    :param bed_file_path: string - path to bed file
    :param output_prefix: string - prefix for output files
    :param noise_threshold: float - threshold past which to exclude positions from noise calculation
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

    # Filter to only positions below noise threshold
    thresh_boolv = pileup_df_all.apply(_apply_threshold, axis=1, thresh=noise_threshold)
    below_thresh_positions = pileup_df_all[thresh_boolv]
    # For noise from N's
    thresh_boolv_n = pileup_df_all.apply(_apply_threshold, axis=1, thresh=noise_threshold, with_n=True)
    below_thresh_positions_n = pileup_df_all[thresh_boolv_n]
    # For noise from Deletions
    thresh_boolv_del = pileup_df_all.apply(_apply_threshold, axis=1, thresh=noise_threshold, with_del=True)
    below_thresh_positions_del = pileup_df_all[thresh_boolv_del]

    # Calculate sample noise and contributing sites for SNV / insertions
    noisy_positions = _create_noisy_positions_file(below_thresh_positions, output_prefix)
    noisy_positions.to_csv(output_prefix + OUTPUT_NOISE_FILENAME, sep='\t', index=False)
    contributing_sites = noisy_positions.shape[0]
    alt_count_total = below_thresh_positions[ALT_COUNT].sum()
    geno_count_total = below_thresh_positions[GENO_COUNT].sum()
    noise = alt_count_total / (alt_count_total + geno_count_total + EPSILON)
    # For N's
    noisy_positions_n = _create_noisy_positions_file(below_thresh_positions, output_prefix, use_n=True)
    contributing_sites_n = noisy_positions_n.shape[0]
    alt_count_total_n = below_thresh_positions_n['N'].sum()
    geno_count_total_n = below_thresh_positions_n[GENO_COUNT].sum()
    noise_n = alt_count_total_n / (alt_count_total_n + geno_count_total_n + EPSILON)
    # For Deletions
    noisy_positions_del = _create_noisy_positions_file(below_thresh_positions, output_prefix, use_del=True)
    contributing_sites_del = noisy_positions_del.shape[0]
    alt_count_total_del = below_thresh_positions_del['deletions'].sum()
    geno_count_total_del = below_thresh_positions_del[GENO_COUNT].sum()
    noise_del = alt_count_total_del / (alt_count_total_del + geno_count_total_del + EPSILON)

    _write_noise_file(NOISE_ACGT, alt_count_total, geno_count_total, noise, contributing_sites, output_prefix)
    _write_noise_file(NOISE_N, alt_count_total_n, geno_count_total_n, noise_n, contributing_sites_n, output_prefix)
    _write_noise_file(NOISE_DEL, alt_count_total_del, geno_count_total_del, noise_del, contributing_sites_del, output_prefix)

    return noise


def _write_noise_file(output_filename: str, alt: int, geno: int, noise: float, contributing_sites, output_prefix: str = 'sample_id') -> None:
    """
    Save sample noise info to file

    :param alt:
    :param geno:
    :param noise:
    :return:
    """
    pd.DataFrame({
        SAMPLE_ID: [output_prefix],
        ALT_COUNT: [alt],
        GENO_COUNT: [geno],
        NOISE_FRACTION: [noise],
        CONTRIBUTING_SITES: [contributing_sites]
    }).to_csv(output_prefix + output_filename, sep='\t', index=False)


def _create_noisy_positions_file(pileup_df: pd.DataFrame, output_prefix: str = '', use_n: bool = False,
                                 use_del: bool = False) -> None:
    """
    Filter to only positions with noise and save to a tsv
    """
    noisy_boolv = (pileup_df[ALT_COUNT] > 0) | (pileup_df['insertions'] > 0)
    if use_del:
        noisy_boolv = (pileup_df['deletions'] > 0)
    if use_n:
        noisy_boolv = (pileup_df['N'] > 0)

    noisy_positions = pileup_df[noisy_boolv]
    noisy_positions = noisy_positions.sort_values(ALT_COUNT)
    return noisy_positions


def _apply_threshold(row: pd.Series, thresh: float, with_n: bool = False, with_del: bool = False) -> bool:
    """
    Returns False if any alt allele crosses `thresh` for the given row of the pileup, True otherwise

    :param row: pandas.Series - row that represents single pileup position
    :param thresh: float - threshold past which alt allele fraction should return false
    """
    base_counts = {'A': row['A'], 'C': row['C'], 'G': row['G'], 'T': row['T']}
    if with_del:
        base_counts['deletions'] = row['deletions']
    if with_n:
        base_counts['N'] = row['N']
    genotype = max(base_counts, key=base_counts.get)

    non_geno_bases = ['A', 'C', 'G', 'T']
    if with_del:
        non_geno_bases.append('deletions')
    if with_n:
        non_geno_bases.append('N')
    non_geno_bases.remove(genotype)

    tot = row['A'] + row['C'] + row['G'] + row['T']
    if with_del:
        tot += row['deletions']

    if any([row[r] / (tot + EPSILON) > thresh for r in non_geno_bases]):
        return False

    return True


def _calculate_alt_and_geno(noise_df: pd.DataFrame) -> pd.DataFrame:
    """
    Determine the genotype and alt count for each position in the `noise_df`

    :param noise_df: pd.DataFrame
    :return: pd.DataFrame
    """
    noise_df['total_acgt'] = noise_df['A'] + noise_df['C'] + noise_df['G'] + noise_df['T']
    noise_df[GENO_COUNT] = noise_df[['A', 'C', 'G', 'T']].max(axis=1)
    noise_df[ALT_COUNT] = noise_df['total_acgt'] - noise_df[GENO_COUNT]
    return noise_df
