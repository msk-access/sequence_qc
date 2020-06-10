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

noise_output_columns = [
    'chrom',
    'pos',
    'ref',
    'A',
    'C',
    'G',
    'T',
    'insertions',
    'deletions',
]


def load_bed_file(bed_file_path):
    """
    Use pybedtools to read bed file

    :param bed_file_path:
    :return:
    """
    return BedTool(bed_file_path)


def calculate_noise(
    ref_fasta,
    bam_path,
    bed_file_path,
    noise_threshold,
    noise_output_filename='noise_positions.tsv',
    include_insertions=False,
    include_deletions=False,
    include_N=False,
    truncate=True,
    flag_filter=0,
    min_mapping_quality=1,
    min_base_quality=20):
    """
    Create file of noise across specified regions in `bed_file` using pybedtools and pysamstats

    :param ref_fasta:
    :param bam_path:
    :param bed_file_path:
    :param noise_threshold:
    :param include_insertions:
    :param include_deletions:
    :param include_N:
    :param truncate:
    :param flag_filter:
    :param min_mapping_quality:
    :param min_base_quality:
    :return:
    """
    bed_file = load_bed_file(bed_file_path)
    bam = AlignmentFile(bam_path)
    pileup_df_all = pd.DataFrame()

    # Build data frame of all positions in bed file
    for region in bed_file.intervals:
        chrom = region.chrom.replace('chr', '')
        start = region.start
        stop = region.stop

        pileup = pysamstats.load_pileup(
            'variation',
            bam,
            chrom=chrom,
            start=start,
            end=stop,
            fafile=ref_fasta,
            truncate=truncate,
            max_depth=30000,
            # no_del=not add_indels,
            min_baseq=min_base_quality,
            min_mapq=min_mapping_quality)

        pileup_df_all = pd.concat([pileup_df_all, pd.DataFrame(pileup)])

    # Determine per-position genotype and alt count
    pileup_df_all = _calculate_alt_and_geno(pileup_df_all)

    # Include columns for ins / dels / N
    pileup_df_all = _include_indels_and_n_noise(pileup_df_all)

    # Filter to only noisy positions
    boolv = pileup_df_all.apply(_apply_threshold, axis=1, thresh=noise_threshold)
    noise_positions = pileup_df_all[boolv]
    # Convert bytes objects to strings so output tsv is formatted correctly
    noise_positions.loc[:,'chrom'] = noise_positions['chrom'].apply(lambda s: s.decode('utf-8'))
    noise_positions.loc[:,'ref'] = noise_positions['ref'].apply(lambda s: s.decode('utf-8'))
    noise_positions[noise_output_columns].to_csv(noise_output_filename, sep='\t', index=False)

    # Calculate sample noise
    alt_count_total = noise_positions['alt_count'].sum()
    geno_count_total = noise_positions['geno_count'].sum()
    noise = alt_count_total / (alt_count_total + geno_count_total + EPSILON)

    logger.info('Alt count, Geno count, Noise: {} {} {}'.format(alt_count_total, geno_count_total, noise))
    return noise


def _apply_threshold(row, thresh):
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


def _calculate_alt_and_geno(noise_df):
    """
    Determine the genotype and alt count for each position in the `noise_df`

    :param noise_df: pd.DataFrame
    :return: pd.DataFrame
    """
    total_acgt = lambda row: row['A'] + row['C'] + row['G'] + row['T']
    noise_df['total_acgt'] = noise_df.apply(total_acgt, axis=1)
    geno_lambda = lambda row: max(row['A'], row['C'], row['G'], row['T'])
    noise_df['geno_count'] = noise_df.apply(geno_lambda, axis=1)
    alt_lambda = lambda row: row['total_acgt'] - row['geno_count']
    noise_df['alt_count'] = noise_df.apply(alt_lambda, axis=1)
    return noise_df


def _include_indels_and_n_noise(noise_df):
    """
    Add additional columns for noise including insertions / deletions / all indels / N

    :param noise_df: pd.DataFrame
    :return:
    """
    # 1. Noise including insertions as possible genotype or alt allele
    noise_df['total_acgt_ins'] = noise_df['total_acgt'] + noise_df['insertions']
    geno_lambda = lambda row: max(row['A'], row['C'], row['G'], row['T'], row['insertions'])
    noise_df['geno_count_ins'] = noise_df.apply(geno_lambda, axis=1)
    alt_lambda = lambda row: row['total_acgt_ins'] - row['geno_count_ins']
    noise_df['alt_count_ins'] = noise_df.apply(alt_lambda, axis=1)

    # 2. Noise including deletions as possible genotype or alt allele
    noise_df['total_acgt_del'] = noise_df['total_acgt'] + noise_df['deletions']
    geno_lambda = lambda row: max(row['A'], row['C'], row['G'], row['T'], row['deletions'])
    noise_df['geno_count_del'] = noise_df.apply(geno_lambda, axis=1)
    alt_lambda = lambda row: row['total_acgt_del'] - row['geno_count_del']
    noise_df['alt_count_del'] = noise_df.apply(alt_lambda, axis=1)

    # 3. Noise including insertions and deletions as possible genotype (either ins or dels) or alt allele
    noise_df['total_acgt_indel'] = noise_df['total_acgt'] + noise_df['insertions'] + noise_df['deletions']
    geno_lambda = lambda row: max(row['A'], row['C'], row['G'], row['T'], row['insertions'], row['deletions'])
    noise_df['geno_count_indel'] = noise_df.apply(geno_lambda, axis=1)
    alt_lambda = lambda row: row['total_acgt_indel'] - row['geno_count_indel']
    noise_df['alt_count_indel'] = noise_df.apply(alt_lambda, axis=1)

    # 4. Noise including N bases as alt allele (but won't ever be considered as genotype)
    noise_df['total_acgt_N'] = noise_df['total_acgt'] + noise_df['N']
    alt_lambda = lambda row: row['total_acgt_N'] - row['geno_count']
    noise_df['alt_count_N'] = noise_df.apply(alt_lambda, axis=1)

    return noise_df
