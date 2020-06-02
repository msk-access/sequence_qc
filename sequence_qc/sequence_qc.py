import logging
import pysamstats

from pysam import AlignmentFile
from pybedtools import BedTool
from pyfaidx import Fasta


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



def calculate_noise(
    ref_fasta,
    bam_path,
    bed_file_path,
    noise_threshold,
    add_indels=False,
    truncate=True,
    ignore_overlaps=True,
    flag_filter=0,
    min_mapping_quality=1,
    min_base_quality=20):
    """
    Create file of noise across specified regions in `bed_file` using pybedtools and pysam

    :param ref_fasta: string - Path to reference fasta which will be used to compare to
        `bam_path`.
    :param bam_path: string - Path to a BAM file which will be interrogated for noise
        level.
    :param bed_file_path: string - Path to a bed-format file with positions
        over which to calculate the noise level.
    :param noise_threshold: float - Upper limit to consider a position "noisy".
        This should be a fraction between 0 and 1, which will determine whether a
        position in the supplied bed_file with any alt reads has too many to be considered
        just noise.

    :param add_indels: see pysam.Pileup
    :param truncate: see pysam.Pileup
    :param ignore_overlaps: see pysam.Pileup
    :param flag_filter: see pysam.Pileup

    :return:
    """
    # todo: why do we need to use check_seq?
    ref = Fasta(ref_fasta)
    af = AlignmentFile(bam_path, check_sq=False)
    bed_file = load_bed_file(bed_file_path)

    alt_count = 0
    # Use small number to avoid DivisionByZero
    total_count = EPSILON

    for region in bed_file.intervals:
        # todo: shadows builtin chr()
        chr = region.chrom.replace('chr', '')
        start = region.start
        stop = region.stop

        logger.debug("Processing region {}, start: {}, end: {}".format(chr, start, stop))
        # todo: why do we need to use "nofilter" to get any pileups...?
        # todo: why doesn't min_mapq do anything?
        pileup = af.pileup(chr, start, stop, stepper="nofilter", truncate=truncate, ignore_overlaps=ignore_overlaps, flag_filter=flag_filter)

        for p in pileup:
            logger.debug("Position: {}".format(p.pos))

            refbase = ref[region.chrom][p.pos:p.pos + 1]
            # todo: need this?
            refbase = str(refbase)
            refbase_lower = refbase.lower()
            logger.debug("Ref Base: {}".format(refbase))

            bases = p.get_query_sequences(add_indels=add_indels)
            # Even with using the add_indels=False, deletions will still be returned as empty strings,
            # thus we must remove them from the bases manually
            # todo: confirm this is correct pysam behavior
            if not add_indels:
                bases = list(filter(''.__ne__, bases))
            logger.debug("Pileup: {}".format(bases))

            # Apply mapping quality filter
            mapping_qualities = p.get_mapping_qualities()
            base_qualities = p.get_query_qualities()
            bases = [b for i, b in enumerate(bases) if mapping_qualities[i] > min_mapping_quality]

            # Need to filter base qualities as well
            base_qualities = [b for i, b in enumerate(base_qualities) if mapping_qualities[i] > min_mapping_quality]
            # Apply base quality filter
            bases = [b for i, b in enumerate(bases) if base_qualities[i] > min_base_quality]

            # todo: instead of comparing to both upper and lowercase, try to use samtools "." and "," formatting
            mismatches = list(filter((refbase).__ne__, bases))
            mismatches = list(filter((refbase_lower).__ne__, mismatches))
            mismatches = [m.upper() for m in mismatches]
            mismatches_count = len(mismatches)
            # Use small number to avoid DivisionByZero
            total_base_count = len(bases) + EPSILON
            logger.debug("Mismatches: {}".format(str(mismatches_count)))

            mismatches_A = list(filter(lambda x: x == 'A', mismatches))
            mismatches_C = list(filter(lambda x: x == 'C', mismatches))
            mismatches_G = list(filter(lambda x: x == 'G', mismatches))
            mismatches_T = list(filter(lambda x: x == 'T', mismatches))
            mismatches_D = list(filter(lambda x: x == '', mismatches))
            mismatches_N = list(filter(lambda x: x == 'N', mismatches))
            mismatches_all = [len(mismatches_A), len(mismatches_C), len(mismatches_G), len(mismatches_T), len(mismatches_D), len(mismatches_N)]
            logger.debug("Mismatches A, C, G, T, D, N: {}".format(mismatches_all))

            if all([((m / total_base_count) < noise_threshold) for m in mismatches_all]):
                alt_count += mismatches_count
                total_count += total_base_count

    logger.debug("Alt count: {}".format(alt_count))
    logger.debug("Total base count: {}".format(total_count))
    return alt_count / total_count


def calculate_noise_pysamstats(
    ref_fasta,
    bam_path,
    bed_file_path,
    noise_threshold,
    add_indels=False,
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
            no_del=not add_indels,
            min_baseq=min_base_quality,
            min_mapq=min_mapping_quality)

        for rec in position_coverage:
            print(rec)
            nonref_bases = {'A', 'C', 'G', 'T'}
            nonref_bases.remove(rec['ref'])

            if all([rec[r] / (rec['reads_all'] + EPSILON) < noise_threshold for r in nonref_bases]):
                alt_count += rec['mismatches']
                total_count += rec['reads']

    logging.debug("Alt base count: {}, Total base count: {}".format(alt_count, total_count))
    return alt_count / total_count
