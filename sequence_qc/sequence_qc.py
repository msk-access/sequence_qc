import logging

from pysam import AlignmentFile
from pybedtools import BedTool
from pyfaidx import Fasta


FORMAT = '%(asctime)-15s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger("sequence_qc")
logger.setLevel(logging.DEBUG)

def calculate_noise(
    ref_fasta,
    bam_path,
    bed_file_path,
    noise_threshold,
    add_indels=False,
    truncate=True,
    ignore_overlaps=True,
    flag_filter=0):
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
    bed_file = BedTool(bed_file_path)

    alt_count = 0
    # Use small number to avoid DivisionByZero
    total_count = 1e-9

    for region in bed_file.intervals:
        # todo: shadows builtin chr()
        chr = region.chrom.replace('chr', '')
        start = region.start
        stop = region.stop

        logger.debug("Processing region {}, start: {}, end: {}".format(chr, start, stop))
        # todo: why do we need to use "nofilter" to get any pileups...?
        pileup = af.pileup(chr, start, stop, stepper="nofilter", truncate=truncate, ignore_overlaps=ignore_overlaps, flag_filter=flag_filter)

        for p in pileup:
            logger.debug("Position: {}".format(p.pos))

            refbase = ref[region.chrom][p.pos:p.pos + 1]
            # todo: need this?
            refbase = str(refbase)
            refbase_lower = refbase.lower()
            logger.debug("Ref Base: {}".format(refbase))

            bases = p.get_query_sequences(add_indels=add_indels)
            logger.debug("Pileup: {}".format(''.join(bases)))

            # todo: instead of comparing to both upper and lowercase, try to use samtools "." and "," formatting
            mismatches = list(filter((refbase).__ne__, bases))
            mismatches = list(filter((refbase_lower).__ne__, mismatches))
            mismatches = [m.upper() for m in mismatches]
            mismatches_count = len(mismatches)
            # Use small number to avoid DivisionByZero
            total_base_count = len(bases) + 1e-9
            logger.debug("Mismatches: {}".format(str(mismatches_count)))

            mismatches_A = list(filter(lambda x: x == 'A', mismatches))
            mismatches_C = list(filter(lambda x: x == 'C', mismatches))
            mismatches_G = list(filter(lambda x: x == 'G', mismatches))
            mismatches_T = list(filter(lambda x: x == 'T', mismatches))
            mismatches_all = [len(mismatches_A), len(mismatches_C), len(mismatches_G), len(mismatches_T)]
            logger.debug("Mismatches A, C, G, T: {}".format(mismatches_all))

            if all([((m / total_base_count) < noise_threshold) for m in mismatches_all]):
                alt_count += mismatches_count
                total_count += total_base_count

    return alt_count / total_count
