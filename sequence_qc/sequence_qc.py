import logging
import pysamstats

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
    ignore_overlaps=True,
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
