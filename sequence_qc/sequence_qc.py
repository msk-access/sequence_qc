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


def load_bed_file(bed_file_path):
    """
    Use pybedtools to read bed file

    :param bed_file_path:
    :return:
    """
    return BedTool(bed_file_path)


def calculate_noise_pysamstats(
    ref_fasta,
    bam_path,
    bed_file_path,
    noise_threshold,
    include_insertions=False,
    include_deletions=False,
    include_N=False,
    truncate=True,
    flag_filter=0,
    min_mapping_quality=1,
    min_base_quality=20):
    """
    Same calculation as above, using pysamstats

    :return:
    """
    bed_file = load_bed_file(bed_file_path)
    mybam = AlignmentFile(bam_path)

    alt_count = 0
    total_count = EPSILON

    for region in bed_file.intervals:
        # todo: shadows builtin chr()
        chr = region.chrom.replace('chr', '')
        start = region.start
        stop = region.stop

        position_coverage = pysamstats.stat_pileup(
            'variation',
            mybam,
            chrom=chr,
            start=start,
            end=stop,
            fafile=ref_fasta,
            truncate=truncate,
            max_depth=30000,
            # no_del=not add_indels,
            min_baseq=min_base_quality,
            min_mapq=min_mapping_quality)

        for rec in position_coverage:
            logger.debug(rec)

            base_counts = {'A': rec['A'], 'C': rec['C'], 'G': rec['G'], 'T': rec['T']}
            genotype = max(base_counts, key=base_counts.get)
            non_geno_bases = ['A', 'C', 'G', 'T']
            non_geno_bases.remove(genotype)

            if include_insertions:
                non_geno_bases.append('insertions')
            if include_deletions:
                non_geno_bases.append('deletions')
            if include_N:
                non_geno_bases.append('N')

            all_reads = rec['A'] + rec['C'] + rec['G'] + rec['T']

            if include_insertions:
                all_reads += rec['insertions']
            if include_deletions:
                all_reads += rec['deletions']
            if include_N:
                all_reads += rec['N']

            all_reads_div = all_reads + EPSILON

            if all([rec[r] / all_reads_div < noise_threshold for r in non_geno_bases]):
                p_alt_count = sum([rec[r] for r in non_geno_bases])

                if include_insertions:
                    p_alt_count += rec['insertions']
                if include_deletions:
                    p_alt_count += rec['deletions']
                if include_N:
                    p_alt_count += rec['N']

                alt_count += p_alt_count
                total_count += all_reads

    logger.debug("Alt base count: {}, Total base count: {}".format(alt_count, total_count))
    return alt_count / total_count


def calculate_noise_pandas(
    ref_fasta,
    bam_path,
    bed_file_path,
    noise_threshold,
    include_insertions=False,
    include_deletions=False,
    include_N=False,
    truncate=True,
    flag_filter=0,
    min_mapping_quality=1,
    min_base_quality=20):
    """


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
    total_acgt = lambda row: row['A'] + row['C'] + row['G'] + row['T']
    pileup_df_all['total_acgt'] = pileup_df_all.apply(total_acgt, axis=1)
    geno_lambda = lambda row: max(row['A'], row['C'], row['G'], row['T'])
    pileup_df_all['geno_count'] = pileup_df_all.apply(geno_lambda, axis=1)
    alt_lambda = lambda row: row['total_acgt'] - row['geno_count']
    pileup_df_all['alt_count'] = pileup_df_all.apply(alt_lambda, axis=1)

    # Filter to only noisy positions
    noise_frac = pileup_df_all['alt_count'] / (pileup_df_all['alt_count'] + pileup_df_all['geno_count'])
    boolv = noise_frac < noise_threshold
    noisy_positions = pileup_df_all[boolv]

    # Calculate sample noise
    alt_count_total = noisy_positions['alt_count'].sum()
    geno_count_total = noisy_positions['geno_count'].sum()
    noise = alt_count_total / (alt_count_total + geno_count_total)

    logger.info('Alt count, Geno count, Noise: {} {} {}'.format(alt_count_total, geno_count_total, noise))
    return noise
