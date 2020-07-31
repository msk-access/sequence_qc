import logging
import pysamstats
import pandas as pd
from pysam import AlignmentFile
from pybedtools import BedTool

from sequence_qc import plots


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
N_COUNT = 'n_count'
DEL_COUNT = 'del_count'
TOTAL_BASE_COUNT = 'total_base_count'
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
                    min_mapping_quality: int = 1, min_base_quality: int = 1, sample_id: str = '',) -> float:
    """
    Create file of noise across specified regions in `bed_file` using pybedtools and pysamstats

    :param ref_fasta: string - path to reference fastq
    :param bam_path: string - path to bam
    :param bed_file_path: string - path to bed file
    :param sample_id: string - prefix for output files
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
    pileup_df_all[output_columns].to_csv(sample_id + OUTPUT_PILEUP_NAME, sep='\t', index=False)

    # Continue with calculation
    noise = _calculate_noise_from_pileup(pileup_df_all, sample_id, noise_threshold)
    return noise


def _calculate_noise_from_pileup(pileup: pd.DataFrame, sample_id: str, noise_threshold: float) -> float:
    """
    Use the pileup to determine average noise, and create noise output files

    :param pileup:
    :return:
    """
    pileup_df_all = _calculate_alt_and_geno(pileup)

    # Calculate sample noise and contributing sites for SNV / insertions
    #
    # Filter to only positions below noise threshold
    thresh_boolv = pileup_df_all.apply(_apply_threshold, axis=1, thresh=noise_threshold)
    below_thresh_positions = pileup_df_all[thresh_boolv]
    noisy_boolv = (below_thresh_positions[ALT_COUNT] > 0) | (below_thresh_positions['insertions'] > 0)
    noisy_positions = below_thresh_positions[noisy_boolv]
    noisy_positions = noisy_positions.sort_values(ALT_COUNT, ascending=False)

    noisy_positions.to_csv(sample_id + OUTPUT_NOISE_FILENAME, sep='\t', index=False)
    contributing_sites = noisy_positions.shape[0]
    alt_count_total = below_thresh_positions[ALT_COUNT].sum()
    geno_count_total = below_thresh_positions[GENO_COUNT].sum()
    noise = alt_count_total / (alt_count_total + geno_count_total + EPSILON)

    pd.DataFrame({
        SAMPLE_ID: [sample_id],
        ALT_COUNT: [alt_count_total],
        GENO_COUNT: [geno_count_total],
        NOISE_FRACTION: [noise],
        CONTRIBUTING_SITES: [contributing_sites]
    }).to_csv(sample_id + NOISE_ACGT, sep='\t', index=False)

    # For noise from Deletions
    thresh_lambda = lambda row: (row['deletions'] / (row['total_acgt'] + row['deletions'] + EPSILON)) < noise_threshold
    thresh_boolv_del = pileup_df_all.apply(thresh_lambda, axis=1)
    below_thresh_positions_del = pileup_df_all[thresh_boolv_del]
    noisy_positions_del = below_thresh_positions_del[below_thresh_positions_del['deletions'] > 0]
    contributing_sites_del = noisy_positions_del.shape[0]
    alt_count_total_del = below_thresh_positions_del['deletions'].sum()
    total_count_del = alt_count_total_del + below_thresh_positions_del['total_acgt'].sum()
    noise_del = alt_count_total_del / (total_count_del + EPSILON)

    pd.DataFrame({
        SAMPLE_ID: [sample_id],
        DEL_COUNT: [alt_count_total_del],
        TOTAL_BASE_COUNT: [total_count_del],
        NOISE_FRACTION: [noise_del],
        CONTRIBUTING_SITES: [contributing_sites_del]
    }).to_csv(sample_id + NOISE_DEL, sep='\t', index=False)

    # For N's
    noisy_positions_n = pileup_df_all[pileup_df_all['N'] > 0]
    contributing_sites_n = noisy_positions_n.shape[0]
    total_n = noisy_positions_n['N'].sum()
    total_acgt = pileup_df_all['total_acgt'].sum()
    noise_n = total_n / (total_n + total_acgt + EPSILON)

    # Make plots
    plots.all_plots(pileup_df_all, noisy_positions, sample_id)

    pd.DataFrame({
        SAMPLE_ID: [sample_id],
        N_COUNT: [total_n],
        TOTAL_BASE_COUNT: [total_acgt],
        NOISE_FRACTION: [noise_n],
        CONTRIBUTING_SITES: [contributing_sites_n]
    }).to_csv(sample_id + NOISE_N, sep='\t', index=False)

    return noise


def _apply_threshold(row: pd.Series, thresh: float, with_del: bool = False) -> bool:
    """
    Returns False if any alt allele crosses `thresh` for the given row of the pileup, True otherwise

    :param row: pandas.Series - row that represents single pileup position
    :param thresh: float - threshold past which alt allele fraction should return false
    """
    base_counts = {'A': row['A'], 'C': row['C'], 'G': row['G'], 'T': row['T']}
    genotype = max(base_counts, key=base_counts.get)

    non_geno_bases = ['A', 'C', 'G', 'T']
    non_geno_bases.remove(genotype)

    tot = row['A'] + row['C'] + row['G'] + row['T']

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
    noise_df['noise_acgt'] = noise_df[ALT_COUNT] / noise_df['total_acgt']
    return noise_df
